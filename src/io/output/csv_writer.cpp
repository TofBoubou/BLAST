#include "blast/io/output/csv_writer.hpp"
#include <algorithm>
#include <sstream>

namespace blast::io::output {

// CSVWriter implementation
auto CSVWriter::write(
    const std::filesystem::path& file_path,
    const OutputDataset& dataset,
    const OutputConfig& config,
    ProgressCallback progress
) const -> std::expected<void, OutputError> {
    
    try {
        if (progress) progress(0.0, "Starting CSV write");
        
        // Create output directory
        auto output_dir = file_path.parent_path();
        std::filesystem::create_directories(output_dir);
        
        if (csv_config_.separate_files_per_station) {
            return write_station_files(file_path, dataset, config, progress);
        } else if (csv_config_.separate_files_per_variable) {
            return write_variable_files(file_path, dataset, config);
        } else {
            return write_single_file(file_path, dataset, config);
        }
        
    } catch (const std::exception& e) {
        return std::unexpected(OutputError(
            std::format("CSV write failed: {}", e.what())
        ));
    }
}

auto CSVWriter::write_station_files(
    const std::filesystem::path& base_path,
    const OutputDataset& dataset,
    const OutputConfig& config,
    ProgressCallback progress
) const -> std::expected<void, OutputError> {
    
    auto base_name = base_path.stem().string();
    auto output_dir = base_path.parent_path() / base_name;
    std::filesystem::create_directories(output_dir);
    
    // Write metadata file
    auto metadata_path = output_dir / "metadata.txt";
    if (auto result = write_metadata_file(metadata_path, dataset.metadata); !result) {
        return std::unexpected(result.error());
    }
    
    // Write station files
    for (std::size_t i = 0; i < dataset.stations.size(); ++i) {
        if (progress) {
            double station_progress = static_cast<double>(i) / dataset.stations.size();
            progress(station_progress, std::format("Writing station {}", i));
        }
        
        auto station_filename = std::format("station_{:04d}.csv", i);
        auto station_path = output_dir / station_filename;
        
        if (auto result = write_station_csv(station_path, dataset.stations[i], config, dataset.metadata, i); !result) {
            return std::unexpected(result.error());
        }
    }
    
    // Write summary files
    ProfileCSVWriter profile_writer(csv_config_);
    
    // Wall properties
    auto wall_path = output_dir / "wall_properties.csv";
    if (auto result = profile_writer.write_wall_properties(wall_path, dataset.wall, dataset.metadata.grid.xi_coordinates); !result) {
        return std::unexpected(result.error());
    }
    
    // Integrated quantities
    auto integrated_path = output_dir / "integrated_quantities.csv";
    if (auto result = profile_writer.write_integrated_quantities(integrated_path, dataset.integrated, dataset.metadata.grid.xi_coordinates); !result) {
        return std::unexpected(result.error());
    }
    
    if (progress) progress(1.0, "CSV write complete");
    
    return {};
}

auto CSVWriter::write_station_csv(
    const std::filesystem::path& file_path,
    const StationData& station,
    const OutputConfig& config,
    const SimulationMetadata& metadata,
    std::size_t station_index
) const -> std::expected<void, OutputError> {
    
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(file_path, "Cannot open file for writing"));
    }
    
    // Set output format
    if (csv_config_.scientific_notation) {
        file << std::scientific;
    } else {
        file << std::fixed;
    }
    file << std::setprecision(csv_config_.precision);
    
    // Write header
    if (csv_config_.include_headers) {
        file << "# Station " << station_index << " - xi = " << station.xi << csv_config_.line_ending;
        file << "# x_physical = " << station.x_physical << " m" << csv_config_.line_ending;
        
        auto headers = create_station_header(config, metadata);
        for (std::size_t i = 0; i < headers.size(); ++i) {
            if (i > 0) file << csv_config_.delimiter;
            file << headers[i];
        }
        file << csv_config_.line_ending;
    }
    
    // Write data rows
    const auto n_eta = station.eta.size();
    
    for (std::size_t i = 0; i < n_eta; ++i) {
        std::vector<std::string> row_data;
        
        // Coordinates
        row_data.push_back(format_value(station.eta[i]));
        if (i < station.y_physical.size()) {
            row_data.push_back(format_value(station.y_physical[i]));
        } else {
            row_data.push_back("0.0");
        }
        
        // Flow variables
        if (config.variables.flow_variables) {
            row_data.push_back(format_value(station.F[i]));
            row_data.push_back(format_value(station.g[i]));
            row_data.push_back(format_value(station.V[i]));
            
            if (i < station.temperature.size()) {
                row_data.push_back(format_value(station.temperature[i]));
            }
            if (i < station.pressure.size()) {
                row_data.push_back(format_value(station.pressure[i]));
            }
            if (i < station.density.size()) {
                row_data.push_back(format_value(station.density[i]));
            }
        }
        
        // Species concentrations
        if (config.variables.species_concentrations) {
            for (std::size_t j = 0; j < station.species_concentrations.rows(); ++j) {
                row_data.push_back(format_value(station.species_concentrations(j, i)));
            }
        }
        
        // Transport properties
        if (config.variables.transport_properties && !station.viscosity.empty()) {
            if (i < station.viscosity.size()) {
                row_data.push_back(format_value(station.viscosity[i]));
            }
            if (i < station.thermal_conductivity.size()) {
                row_data.push_back(format_value(station.thermal_conductivity[i]));
            }
        }
        
        // Write row
        for (std::size_t j = 0; j < row_data.size(); ++j) {
            if (j > 0) file << csv_config_.delimiter;
            file << row_data[j];
        }
        file << csv_config_.line_ending;
    }
    
    return {};
}

auto CSVWriter::create_station_header(
    const OutputConfig& config,
    const SimulationMetadata& metadata
) const -> std::vector<std::string> {
    
    std::vector<std::string> headers;
    
    // Coordinates
    headers.push_back("eta");
    headers.push_back("y_physical");
    
    // Flow variables
    if (config.variables.flow_variables) {
        headers.push_back("F");
        headers.push_back("g");
        headers.push_back("V");
        headers.push_back("Temperature");
        headers.push_back("Pressure");
        headers.push_back("Density");
    }
    
    // Species concentrations
    if (config.variables.species_concentrations) {
        for (std::size_t i = 0; i < metadata.mixture.species_names.size(); ++i) {
            headers.push_back(std::format("c_{}", metadata.mixture.species_names[i]));
        }
    }
    
    // Transport properties
    if (config.variables.transport_properties) {
        headers.push_back("Viscosity");
        headers.push_back("ThermalConductivity");
    }
    
    return headers;
}

auto CSVWriter::format_value(double value) const -> std::string {
    std::ostringstream oss;
    
    if (csv_config_.scientific_notation) {
        oss << std::scientific;
    } else {
        oss << std::fixed;
    }
    oss << std::setprecision(csv_config_.precision) << value;
    
    return oss.str();
}

auto CSVWriter::write_metadata_file(
    const std::filesystem::path& file_path,
    const SimulationMetadata& metadata
) const -> std::expected<void, OutputError> {
    
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(file_path, "Cannot open metadata file"));
    }
    
    file << "# BLAST Boundary Layer Solution Metadata" << csv_config_.line_ending;
    file << "# Generated: " << std::format("{:%Y-%m-%d %H:%M:%S}", metadata.creation_time) << csv_config_.line_ending;
    file << csv_config_.line_ending;
    
    file << "BLAST_Version: " << metadata.blast_version << csv_config_.line_ending;
    file << csv_config_.line_ending;
    
    // Grid information
    file << "# Grid Information" << csv_config_.line_ending;
    file << "n_eta: " << metadata.grid.n_eta << csv_config_.line_ending;
    file << "eta_max: " << metadata.grid.eta_max << csv_config_.line_ending;
    file << "d_eta: " << metadata.grid.d_eta << csv_config_.line_ending;
    file << "n_stations: " << metadata.grid.xi_coordinates.size() << csv_config_.line_ending;
    file << csv_config_.line_ending;
    
    // Mixture information
    file << "# Mixture Information" << csv_config_.line_ending;
    file << "n_species: " << metadata.mixture.species_names.size() << csv_config_.line_ending;
    file << "has_electrons: " << (metadata.mixture.has_electrons ? "true" : "false") << csv_config_.line_ending;
    file << "species_names: ";
    for (std::size_t i = 0; i < metadata.mixture.species_names.size(); ++i) {
        if (i > 0) file << csv_config_.delimiter;
        file << metadata.mixture.species_names[i];
    }
    file << csv_config_.line_ending;
    file << csv_config_.line_ending;
    
    // Convergence information
    file << "# Convergence Information" << csv_config_.line_ending;
    file << "converged: " << (metadata.convergence.converged ? "true" : "false") << csv_config_.line_ending;
    file << "total_iterations: " << metadata.convergence.total_iterations << csv_config_.line_ending;
    file << "final_residual: " << metadata.convergence.final_residual << csv_config_.line_ending;
    
    return {};
}

// ProfileCSVWriter implementation
auto ProfileCSVWriter::write_profiles(
    const std::filesystem::path& file_path,
    const StationData& station,
    const std::vector<std::string>& species_names,
    std::size_t station_index
) const -> std::expected<void, OutputError> {
    
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(file_path, "Cannot open file for writing"));
    }
    
    // Set output format
    if (config_.scientific_notation) {
        file << std::scientific;
    } else {
        file << std::fixed;
    }
    file << std::setprecision(config_.precision);
    
    // Write header
    if (config_.include_headers) {
        file << "# Boundary Layer Profiles - Station " << station_index << config_.line_ending;
        file << "# xi = " << station.xi << ", x_physical = " << station.x_physical << " m" << config_.line_ending;
        
        file << "eta" << config_.delimiter << "y_physical" << config_.delimiter;
        file << "F" << config_.delimiter << "g" << config_.delimiter << "V" << config_.delimiter;
        file << "Temperature" << config_.delimiter << "Pressure" << config_.delimiter << "Density";
        
        for (const auto& species_name : species_names) {
            file << config_.delimiter << "c_" << species_name;
        }
        file << config_.line_ending;
    }
    
    // Write data
    const auto n_eta = station.eta.size();
    for (std::size_t i = 0; i < n_eta; ++i) {
        file << station.eta[i] << config_.delimiter;
        file << (i < station.y_physical.size() ? station.y_physical[i] : 0.0) << config_.delimiter;
        file << station.F[i] << config_.delimiter;
        file << station.g[i] << config_.delimiter;
        file << station.V[i] << config_.delimiter;
        file << (i < station.temperature.size() ? station.temperature[i] : 0.0) << config_.delimiter;
        file << (i < station.pressure.size() ? station.pressure[i] : 0.0) << config_.delimiter;
        file << (i < station.density.size() ? station.density[i] : 0.0);
        
        for (std::size_t j = 0; j < species_names.size(); ++j) {
            file << config_.delimiter;
            if (j < station.species_concentrations.rows()) {
                file << station.species_concentrations(j, i);
            } else {
                file << 0.0;
            }
        }
        file << config_.line_ending;
    }
    
    return {};
}

auto ProfileCSVWriter::write_wall_properties(
    const std::filesystem::path& file_path,
    const OutputDataset::WallData& wall_data,
    const std::vector<double>& xi_coordinates
) const -> std::expected<void, OutputError> {
    
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(file_path, "Cannot open file for writing"));
    }
    
    // Set output format
    if (config_.scientific_notation) {
        file << std::scientific;
    } else {
        file << std::fixed;
    }
    file << std::setprecision(config_.precision);
    
    // Write header
    if (config_.include_headers) {
        file << "# Wall Properties" << config_.line_ending;
        file << "xi" << config_.delimiter << "x_position" << config_.delimiter;
        file << "temperature" << config_.delimiter << "heat_flux" << config_.delimiter;
        file << "shear_stress" << config_.line_ending;
    }
    
    // Write data
    const auto n_points = std::min({xi_coordinates.size(), wall_data.x_positions.size(),
                                   wall_data.temperatures.size(), wall_data.heat_flux.size(),
                                   wall_data.shear_stress.size()});
    
    for (std::size_t i = 0; i < n_points; ++i) {
        file << xi_coordinates[i] << config_.delimiter;
        file << wall_data.x_positions[i] << config_.delimiter;
        file << wall_data.temperatures[i] << config_.delimiter;
        file << wall_data.heat_flux[i] << config_.delimiter;
        file << wall_data.shear_stress[i] << config_.line_ending;
    }
    
    return {};
}

auto ProfileCSVWriter::write_integrated_quantities(
    const std::filesystem::path& file_path,
    const OutputDataset::IntegratedQuantities& quantities,
    const std::vector<double>& xi_coordinates
) const -> std::expected<void, OutputError> {
    
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(file_path, "Cannot open file for writing"));
    }
    
    // Set output format
    if (config_.scientific_notation) {
        file << std::scientific;
    } else {
        file << std::fixed;
    }
    file << std::setprecision(config_.precision);
    
    // Write header
    if (config_.include_headers) {
        file << "# Integrated Boundary Layer Quantities" << config_.line_ending;
        file << "xi" << config_.delimiter << "displacement_thickness" << config_.delimiter;
        file << "momentum_thickness" << config_.delimiter << "shape_factor" << config_.delimiter;
        file << "total_enthalpy_thickness" << config_.line_ending;
    }
    
    // Write data
    const auto n_points = std::min({xi_coordinates.size(), quantities.displacement_thickness.size(),
                                   quantities.momentum_thickness.size(), quantities.shape_factor.size(),
                                   quantities.total_enthalpy_thickness.size()});
    
    for (std::size_t i = 0; i < n_points; ++i) {
        file << xi_coordinates[i] << config_.delimiter;
        file << quantities.displacement_thickness[i] << config_.delimiter;
        file << quantities.momentum_thickness[i] << config_.delimiter;
        file << quantities.shape_factor[i] << config_.delimiter;
        file << quantities.total_enthalpy_thickness[i] << config_.line_ending;
    }
    
    return {};
}

// Convenience functions implementation
namespace csv {

auto write_csv_series(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& config,
    const thermophysics::MixtureInterface& mixture,
    const std::filesystem::path& output_directory
) -> std::expected<std::vector<std::filesystem::path>, OutputError> {
    
    // Create CSV output configuration
    OutputConfig csv_output_config;
    csv_output_config.base_directory = output_directory;
    csv_output_config.primary_format = OutputFormat::CSV;
    csv_output_config.variables.flow_variables = true;
    csv_output_config.variables.species_concentrations = true;
    csv_output_config.save_metadata = true;
    
    // Create CSV writer
    CSVConfig csv_config;
    csv_config.separate_files_per_station = true;
    csv_config.include_headers = true;
    csv_config.precision = 8;
    
    CSVWriter csv_writer(csv_config);
    
    // Convert solution to output dataset (reusing existing conversion logic)
    OutputWriter temp_writer(csv_output_config);
    
    auto write_result = temp_writer.write_solution(solution, config, mixture, "csv_output");
    if (!write_result) {
        return std::unexpected(write_result.error());
    }
    
    return write_result.value();
}

auto write_profiles_only(
    const blast::boundary_layer::solver::SolutionResult& solution,
    const Configuration& config,
    const thermophysics::MixtureInterface& mixture,
    const std::filesystem::path& output_directory
) -> std::expected<std::vector<std::filesystem::path>, OutputError> {
    
    std::filesystem::create_directories(output_directory);
    
    CSVConfig csv_config;
    csv_config.include_headers = true;
    csv_config.precision = 8;
    
    ProfileCSVWriter profile_writer(csv_config);
    
    // Get species names
    std::vector<std::string> species_names;
    for (std::size_t i = 0; i < mixture.n_species(); ++i) {
        species_names.emplace_back(mixture.species_name(i));
    }
    
    std::vector<std::filesystem::path> written_files;
    
    // Write profile for each station
    for (std::size_t i = 0; i < solution.stations.size(); ++i) {
        auto filename = std::format("profiles_station_{:04d}.csv", i);
        auto file_path = output_directory / filename;
        
        // Convert SolutionState to StationData (simplified)
        StationData station_data;
        station_data.xi = solution.xi_solved[i];
        station_data.x_physical = 0.0; // TODO: Compute from xi
        station_data.F = solution.stations[i].F;
        station_data.g = solution.stations[i].g;
        station_data.V = solution.stations[i].V;
        station_data.species_concentrations = solution.stations[i].c;
        
        // Generate eta coordinates
        const double eta_max = config.numerical.eta_max;
        const auto n_eta = solution.stations[i].F.size();
        station_data.eta.resize(n_eta);
        for (std::size_t j = 0; j < n_eta; ++j) {
            station_data.eta[j] = static_cast<double>(j) * eta_max / (n_eta - 1);
        }
        
        // Initialize other fields with defaults
        station_data.y_physical.resize(n_eta, 0.0);
        station_data.temperature.resize(n_eta, 300.0);
        station_data.pressure.resize(n_eta, 1000.0);
        station_data.density.resize(n_eta, 0.01);
        
        auto write_result = profile_writer.write_profiles(file_path, station_data, species_names, i);
        if (!write_result) {
            return std::unexpected(write_result.error());
        }
        
        written_files.push_back(file_path);
    }
    
    return written_files;
}

auto validate_csv_file(
    const std::filesystem::path& file_path,
    char expected_delimiter
) -> std::expected<void, OutputError> {
    
    if (!std::filesystem::exists(file_path)) {
        return std::unexpected(OutputError(
            std::format("CSV file does not exist: {}", file_path.string())
        ));
    }
    
    std::ifstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(OutputError(
            std::format("Cannot open CSV file: {}", file_path.string())
        ));
    }
    
    std::string first_line;
    std::getline(file, first_line);
    
    // Skip comment lines
    while (first_line.empty() || first_line[0] == '#') {
        if (!std::getline(file, first_line)) {
            return std::unexpected(OutputError("CSV file contains no data"));
        }
    }
    
    // Check if delimiter exists
    if (first_line.find(expected_delimiter) == std::string::npos) {
        return std::unexpected(OutputError(
            std::format("Expected delimiter '{}' not found in CSV file", expected_delimiter)
        ));
    }
    
    return {};
}

} // namespace csv

} // namespace blast::io::output