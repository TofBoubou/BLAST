#include "blast/io/output/output_writer.hpp"
#include "blast/io/output/hdf5_writer.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <iostream>

namespace blast::io::output {

// WriterFactory implementation
auto WriterFactory::create_writer(OutputFormat format) 
    -> std::expected<std::unique_ptr<FormatWriter>, UnsupportedFormatError> {
    
    switch (format) {
        case OutputFormat::HDF5:
            return std::make_unique<HDF5Writer>();
        
        default:
            return std::unexpected(UnsupportedFormatError(format));
    }
}

auto WriterFactory::get_available_formats() noexcept -> std::vector<OutputFormat> {
    return {
        OutputFormat::HDF5
        // Add more as implemented
    };
}

// OutputWriter implementation
OutputWriter::OutputWriter(OutputConfig config) : config_(std::move(config)) {
    initialize_writers();
}

auto OutputWriter::initialize_writers() -> void {
    writers_.clear();
    
    // Create primary writer
    if (auto writer = WriterFactory::create_writer(config_.primary_format)) {
        writers_.push_back(std::move(writer.value()));
    }
    
    // Create additional format writers
    for (auto format : config_.additional_formats) {
        if (auto writer = WriterFactory::create_writer(format)) {
            writers_.push_back(std::move(writer.value()));
        }
    }
}

auto OutputWriter::write_solution(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& sim_config,
    const thermophysics::MixtureInterface& mixture,
    const std::string& case_name,
    ProgressCallback progress
) -> std::expected<std::vector<std::filesystem::path>, OutputError> {
    
    // Convert solution to output dataset
    auto dataset_result = convert_solution(solution, sim_config, mixture);
    if (!dataset_result) {
        return std::unexpected(dataset_result.error());
    }
    auto dataset = std::move(dataset_result.value());
    
    // Compute derived quantities
    if (config_.save_metadata) {
        if (auto wall_result = compute_wall_properties(dataset, mixture); !wall_result) {
            return std::unexpected(wall_result.error());
        }
        
        if (auto integrated_result = compute_integrated_quantities(dataset); !integrated_result) {
            return std::unexpected(integrated_result.error());
        }
    }
    
    // Generate file paths
    auto file_paths = generate_file_paths(case_name, dataset.metadata.creation_time);
    
    // Ensure output directory exists
    if (!file_paths.empty()) {
        auto output_dir = file_paths[0].parent_path();
        std::filesystem::create_directories(output_dir);
    }
    
    std::vector<std::filesystem::path> written_files;
    written_files.reserve(writers_.size());
    
    // Write using each format writer
    for (std::size_t i = 0; i < writers_.size() && i < file_paths.size(); ++i) {
        auto& writer = writers_[i];
        const auto& file_path = file_paths[i];
        
        if (progress) {
            progress(static_cast<double>(i) / writers_.size(), 
                    std::format("Writing {}", file_path.filename().string()));
        }
        
        if (auto write_result = writer->write(file_path, dataset, config_, progress); !write_result) {
            return std::unexpected(FileWriteError(file_path, write_result.error().message()));
        }
        
        written_files.push_back(file_path);
    }
    
    if (progress) {
        progress(1.0, "Output complete");
    }
    
    return written_files;
}

auto OutputWriter::convert_solution(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& sim_config,
    const thermophysics::MixtureInterface& mixture
) const -> std::expected<OutputDataset, OutputError> {
    
    OutputDataset dataset;
    
    // Fill metadata
    dataset.metadata.creation_time = std::chrono::system_clock::now();
    dataset.metadata.simulation_config = sim_config;
    dataset.metadata.convergence.converged = solution.converged;
    dataset.metadata.convergence.total_iterations = solution.total_iterations;
    
    // Grid information
    dataset.metadata.grid.n_eta = sim_config.numerical.n_eta;
    dataset.metadata.grid.eta_max = sim_config.numerical.eta_max;
    dataset.metadata.grid.d_eta = sim_config.numerical.eta_max / (sim_config.numerical.n_eta - 1);
    
    // Generate eta coordinates
    dataset.metadata.grid.eta_coordinates.resize(sim_config.numerical.n_eta);
    for (int i = 0; i < sim_config.numerical.n_eta; ++i) {
        dataset.metadata.grid.eta_coordinates[i] = i * dataset.metadata.grid.d_eta;
    }
    
    dataset.metadata.grid.xi_coordinates = solution.xi_solved;
    
    // Mixture information
    dataset.metadata.mixture.species_names.resize(mixture.n_species());
    dataset.metadata.mixture.species_molecular_weights.resize(mixture.n_species());
    for (std::size_t i = 0; i < mixture.n_species(); ++i) {
        dataset.metadata.mixture.species_names[i] = mixture.species_name(i);
        dataset.metadata.mixture.species_molecular_weights[i] = mixture.species_molecular_weight(i);
    }
    
    auto charges = mixture.species_charges();
    dataset.metadata.mixture.species_charges.assign(charges.begin(), charges.end());
    dataset.metadata.mixture.has_electrons = mixture.has_electrons();
    
    // Convert station data
    dataset.stations.reserve(solution.stations.size());
    
    for (std::size_t station_idx = 0; station_idx < solution.stations.size(); ++station_idx) {
        const auto& station_solution = solution.stations[station_idx];
        StationData station_data;
        
        // Basic coordinates
        station_data.xi = solution.xi_solved[station_idx];
        station_data.x_physical = 0.0; // TODO: Compute from xi
        
        station_data.eta = dataset.metadata.grid.eta_coordinates;
        // TODO: Compute physical y coordinates
        station_data.y_physical.resize(station_data.eta.size(), 0.0);
        
        // Primary variables
        station_data.F = station_solution.F;
        station_data.g = station_solution.g;
        station_data.V = station_solution.V;
        
        // Temperature field
        if (station_idx < solution.temperature_fields.size()) {
            station_data.temperature = solution.temperature_fields[station_idx];
        } else {
            station_data.temperature.resize(station_solution.F.size(), 300.0); // Default
            std::cout << "Temperature field size: " << solution.temperature_fields[0].size() << std::endl;
            std::cout << "Erreur dans le traitement de la température" << std::endl;
        }
        
        // Initialize pressure and density with placeholder values
        station_data.pressure.resize(station_solution.F.size(), 1000.0);
        station_data.density.resize(station_solution.F.size(), 0.01);
        
        // Species concentrations
        station_data.species_concentrations = station_solution.c;
        
        // Initialize optional fields based on config
        if (config_.variables.transport_properties) {
            station_data.viscosity.resize(station_solution.F.size(), 0.0);
            station_data.thermal_conductivity.resize(station_solution.F.size(), 0.0);
            station_data.diffusion_coefficients = core::Matrix<double>(mixture.n_species(), mixture.n_species());
            station_data.diffusion_coefficients.setZero();
        }
        
        if (config_.variables.chemical_rates) {
            station_data.production_rates = core::Matrix<double>(mixture.n_species(), station_solution.F.size());
            station_data.production_rates.setZero();
        }
        
        if (config_.variables.diffusion_fluxes) {
            station_data.diffusion_fluxes = core::Matrix<double>(mixture.n_species(), station_solution.F.size());
            station_data.diffusion_fluxes.setZero();
        }
        
        dataset.stations.push_back(std::move(station_data));
    }
    
    return dataset;
}

auto OutputWriter::compute_wall_properties(
    OutputDataset& dataset,
    const thermophysics::MixtureInterface& mixture
) const -> std::expected<void, OutputError> {
    
    dataset.wall.x_positions.reserve(dataset.stations.size());
    dataset.wall.temperatures.reserve(dataset.stations.size());
    dataset.wall.heat_flux.reserve(dataset.stations.size());
    dataset.wall.shear_stress.reserve(dataset.stations.size());
    
    for (const auto& station : dataset.stations) {
        dataset.wall.x_positions.push_back(station.x_physical);
        
        if (!station.temperature.empty()) {
            dataset.wall.temperatures.push_back(station.temperature[0]); // Wall temperature
        }
        
        // TODO: Compute heat flux and shear stress from gradients
        dataset.wall.heat_flux.push_back(0.0);
        dataset.wall.shear_stress.push_back(0.0);
    }
    
    return {};
}

auto OutputWriter::compute_integrated_quantities(
    OutputDataset& dataset
) const -> std::expected<void, OutputError> {
    
    dataset.integrated.displacement_thickness.reserve(dataset.stations.size());
    dataset.integrated.momentum_thickness.reserve(dataset.stations.size());
    dataset.integrated.shape_factor.reserve(dataset.stations.size());
    dataset.integrated.total_enthalpy_thickness.reserve(dataset.stations.size());
    
    for (const auto& station : dataset.stations) {
        // TODO: Implement boundary layer integral calculations
        dataset.integrated.displacement_thickness.push_back(0.0);
        dataset.integrated.momentum_thickness.push_back(0.0);
        dataset.integrated.shape_factor.push_back(0.0);
        dataset.integrated.total_enthalpy_thickness.push_back(0.0);
    }
    
    return {};
}

auto OutputWriter::generate_file_paths(
    const std::string& case_name,
    const std::chrono::system_clock::time_point& timestamp
) const -> std::vector<std::filesystem::path> {
    
    std::vector<std::filesystem::path> paths;
    paths.reserve(writers_.size());
    
    // Generate timestamp string if requested
    std::string timestamp_str;
    if (config_.include_timestamp) {
        auto time_t = std::chrono::system_clock::to_time_t(timestamp);
        auto tm = *std::localtime(&time_t);
        
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
        timestamp_str = "_" + oss.str();
    }
    
    // Create paths for each writer
    for (const auto& writer : writers_) {
        auto filename = case_name + timestamp_str + std::string(writer->get_extension());
        auto full_path = config_.base_directory / filename;
        paths.push_back(full_path);
    }
    
    return paths;
}

auto OutputWriter::validate_config() const -> std::expected<void, OutputError> {
    // Check if base directory is writable
    if (!std::filesystem::exists(config_.base_directory)) {
        std::error_code ec;
        std::filesystem::create_directories(config_.base_directory, ec);
        if (ec) {
            return std::unexpected(OutputError(
                std::format("Cannot create output directory '{}': {}", 
                           config_.base_directory.string(), ec.message())
            ));
        }
    }
    
    // Validate compression level
    if (config_.compression_level < 0.0 || config_.compression_level > 9.0) {
        return std::unexpected(OutputError("Compression level must be between 0.0 and 9.0"));
    }
    
    return {};
}

auto OutputWriter::get_output_info(const std::string& case_name) const 
    -> std::vector<std::pair<OutputFormat, std::filesystem::path>> {
    
    std::vector<std::pair<OutputFormat, std::filesystem::path>> info;
    
    auto timestamp = std::chrono::system_clock::now();
    auto paths = generate_file_paths(case_name, timestamp);
    
    std::vector<OutputFormat> formats = {config_.primary_format};
    formats.insert(formats.end(), config_.additional_formats.begin(), config_.additional_formats.end());
    
    for (std::size_t i = 0; i < std::min(formats.size(), paths.size()); ++i) {
        info.emplace_back(formats[i], paths[i]);
    }
    
    return info;
}

// Convenience functions implementation
namespace convenience {

auto write_hdf5(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& config,
    const thermophysics::MixtureInterface& mixture,
    const std::filesystem::path& output_path,
    bool include_derivatives
) -> std::expected<void, OutputError> {
    
    OutputConfig out_config;
    out_config.primary_format = OutputFormat::HDF5;
    out_config.save_derivatives = include_derivatives;
    out_config.base_directory = output_path.parent_path();
    
    OutputWriter writer(out_config);
    
    auto result = writer.write_solution(solution, config, mixture, output_path.stem().string());
    if (!result) {
        return std::unexpected(result.error());
    }
    
    return {};
}

auto write_vtk(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& config,
    const thermophysics::MixtureInterface& mixture,
    const std::filesystem::path& output_path
) -> std::expected<void, OutputError> {
    
    OutputConfig out_config;
    out_config.primary_format = OutputFormat::HDF5;
    out_config.base_directory = output_path.parent_path();
    
    OutputWriter writer(out_config);
    
    auto result = writer.write_solution(solution, config, mixture, output_path.stem().string());
    if (!result) {
        return std::unexpected(result.error());
    }
    
    return {};
}

auto write_auto(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& config,
    const thermophysics::MixtureInterface& mixture,
    const std::filesystem::path& output_path
) -> std::expected<void, OutputError> {
    
    auto extension = output_path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    
    OutputFormat format;
    if (extension == ".h5" || extension == ".hdf5") {
        format = OutputFormat::HDF5;
    } else {
        return std::unexpected(OutputError(
            std::format("Cannot auto-detect format for extension '{}'", extension)
        ));
    }
    
    OutputConfig out_config;
    out_config.primary_format = format;
    out_config.base_directory = output_path.parent_path();
    
    OutputWriter writer(out_config);
    
    auto result = writer.write_solution(solution, config, mixture, output_path.stem().string());
    if (!result) {
        return std::unexpected(result.error());
    }
    
    return {};
}

} // namespace convenience

} // namespace blast::io::output