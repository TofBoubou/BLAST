#include "blast/io/output/vtk_writer.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <zlib.h>
#include <array>

namespace blast::io::output {

// VTKXMLWriter implementation
auto VTKXMLWriter::write(
    const std::filesystem::path& file_path,
    const OutputDataset& dataset,
    const OutputConfig& config,
    ProgressCallback progress
) const -> std::expected<void, OutputError> {
    
    try {
        if (progress) progress(0.0, "Starting VTK XML write");
        
        if (vtk_config_.write_time_series && dataset.stations.size() > 1) {
            return write_time_series(file_path, dataset, config);
        } else {
            return write_structured_grid(file_path, dataset, config, progress);
        }
        
    } catch (const std::exception& e) {
        return std::unexpected(OutputError(
            std::format("VTK XML write failed: {}", e.what())
        ));
    }
}

auto VTKXMLWriter::write_structured_grid(
    const std::filesystem::path& file_path,
    const OutputDataset& dataset,
    const OutputConfig& config,
    ProgressCallback progress
) const -> std::expected<void, OutputError> {
    
    if (dataset.stations.empty()) {
        return std::unexpected(OutputError("No station data to write"));
    }
    
    std::ofstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(file_path, "Cannot open file for writing"));
    }
    
    const auto& stations = dataset.stations;
    const auto n_stations = stations.size();
    const auto n_eta = stations[0].eta.size();
    
    if (progress) progress(0.1, "Writing VTK header");
    
    // Write XML header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"";
    if (vtk_config_.compress_binary && vtk_config_.data_mode == VTKConfig::DataMode::Binary) {
        file << " compressor=\"vtkZLibDataCompressor\"";
    }
    file << ">\n";
    
    // Write structured grid
    file << "  <StructuredGrid WholeExtent=\"0 " << (n_stations - 1) 
         << " 0 " << (n_eta - 1) << " 0 0\">\n";
    file << "    <Piece Extent=\"0 " << (n_stations - 1) 
         << " 0 " << (n_eta - 1) << " 0 0\">\n";
    
    if (progress) progress(0.3, "Writing grid points");
    
    // Write points (coordinates)
    auto points_result = generate_grid_points(stations);
    if (!points_result) {
        return std::unexpected(points_result.error());
    }
    file << points_result.value();
    
    if (progress) progress(0.6, "Writing point data");
    
    // Write point data (variables)
    auto point_data_result = generate_point_data(stations, config);
    if (!point_data_result) {
        return std::unexpected(point_data_result.error());
    }
    file << point_data_result.value();
    
    if (progress) progress(0.9, "Finalizing VTK file");
    
    // Close tags
    file << "    </Piece>\n";
    file << "  </StructuredGrid>\n";
    file << "</VTKFile>\n";
    
    if (progress) progress(1.0, "VTK XML write complete");
    
    return {};
}

auto VTKXMLWriter::write_time_series(
    const std::filesystem::path& base_path,
    const OutputDataset& dataset,
    const OutputConfig& config
) const -> std::expected<void, OutputError> {
    
    std::vector<std::filesystem::path> vts_files;
    vts_files.reserve(dataset.stations.size());
    
    // Write individual time steps
    for (std::size_t i = 0; i < dataset.stations.size(); ++i) {
        OutputDataset single_station_dataset = dataset;
        single_station_dataset.stations = {dataset.stations[i]};
        
        auto step_name = std::format("{}_{:04d}.vts", base_path.stem().string(), i);
        auto step_path = base_path.parent_path() / step_name;
        
        auto write_result = write_structured_grid(step_path, single_station_dataset, config, nullptr);
        if (!write_result) {
            return std::unexpected(write_result.error());
        }
        
        vts_files.push_back(step_path);
    }
    
    // Create ParaView collection file
    ParaViewCollectionWriter collection_writer(base_path.parent_path() / (base_path.stem().string() + ".pvd"));
    
    for (std::size_t i = 0; i < vts_files.size(); ++i) {
        collection_writer.add_time_step(
            dataset.stations[i].xi,  // Use xi as time
            vts_files[i],
            "boundary_layer"
        );
    }
    
    return collection_writer.write_collection();
}

auto VTKXMLWriter::generate_grid_points(
    const std::vector<StationData>& stations
) const -> std::expected<std::string, OutputError> {
    
    const auto n_stations = stations.size();
    const auto n_eta = stations[0].eta.size();
    const auto total_points = n_stations * n_eta;
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    
    oss << "      <Points>\n";
    
    if (vtk_config_.data_mode == VTKConfig::DataMode::ASCII) {
        oss << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        
        for (std::size_t station_idx = 0; station_idx < n_stations; ++station_idx) {
            const auto& station = stations[station_idx];
            for (std::size_t eta_idx = 0; eta_idx < n_eta; ++eta_idx) {
                // VTK coordinates: (x, y, z)
                // BLAST coordinates: x_physical (streamwise), y_physical (normal), 0 (spanwise)
                double x = station.x_physical;
                double y = (eta_idx < station.y_physical.size()) ? station.y_physical[eta_idx] : 0.0;
                double z = 0.0;
                
                oss << "          " << x << " " << y << " " << z << "\n";
            }
        }
        
        oss << "        </DataArray>\n";
    } else {
        // Binary format
        oss << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">\n";
        
        std::vector<double> coordinates;
        coordinates.reserve(total_points * 3);
        
        for (std::size_t station_idx = 0; station_idx < n_stations; ++station_idx) {
            const auto& station = stations[station_idx];
            for (std::size_t eta_idx = 0; eta_idx < n_eta; ++eta_idx) {
                coordinates.push_back(station.x_physical);
                coordinates.push_back((eta_idx < station.y_physical.size()) ? station.y_physical[eta_idx] : 0.0);
                coordinates.push_back(0.0);
            }
        }
        
        auto encoded_result = encode_binary_data(coordinates);
        if (!encoded_result) {
            return std::unexpected(encoded_result.error());
        }
        
        oss << "          " << encoded_result.value() << "\n";
        oss << "        </DataArray>\n";
    }
    
    oss << "      </Points>\n";
    return oss.str();
}

auto VTKXMLWriter::generate_point_data(
    const std::vector<StationData>& stations,
    const OutputConfig& config
) const -> std::expected<std::string, OutputError> {
    
    const auto n_stations = stations.size();
    const auto n_eta = stations[0].eta.size();
    const auto total_points = n_stations * n_eta;
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    
    oss << "      <PointData>\n";
    
    // Flow variables
    if (config.variables.flow_variables) {
        // F field
        std::vector<double> F_data;
        F_data.reserve(total_points);
        for (const auto& station : stations) {
            F_data.insert(F_data.end(), station.F.begin(), station.F.end());
        }
        oss << generate_field_data("F", F_data);
        
        // g field
        std::vector<double> g_data;
        g_data.reserve(total_points);
        for (const auto& station : stations) {
            g_data.insert(g_data.end(), station.g.begin(), station.g.end());
        }
        oss << generate_field_data("g", g_data);
        
        // V field
        std::vector<double> V_data;
        V_data.reserve(total_points);
        for (const auto& station : stations) {
            V_data.insert(V_data.end(), station.V.begin(), station.V.end());
        }
        oss << generate_field_data("V", V_data);
        
        // Temperature
        std::vector<double> T_data;
        T_data.reserve(total_points);
        for (const auto& station : stations) {
            T_data.insert(T_data.end(), station.temperature.begin(), station.temperature.end());
        }
        oss << generate_field_data("Temperature", T_data);
        
        // Pressure
        std::vector<double> P_data;
        P_data.reserve(total_points);
        for (const auto& station : stations) {
            P_data.insert(P_data.end(), station.pressure.begin(), station.pressure.end());
        }
        oss << generate_field_data("Pressure", P_data);
        
        // Density
        std::vector<double> rho_data;
        rho_data.reserve(total_points);
        for (const auto& station : stations) {
            rho_data.insert(rho_data.end(), station.density.begin(), station.density.end());
        }
        oss << generate_field_data("Density", rho_data);
    }
    
    // Species concentrations
    if (config.variables.species_concentrations && !stations.empty()) {
        const auto n_species = stations[0].species_concentrations.rows();
        
        for (std::size_t species_idx = 0; species_idx < n_species; ++species_idx) {
            std::vector<double> species_data;
            species_data.reserve(total_points);
            
            for (const auto& station : stations) {
                for (std::size_t eta_idx = 0; eta_idx < n_eta; ++eta_idx) {
                    species_data.push_back(station.species_concentrations(species_idx, eta_idx));
                }
            }
            
            auto species_name = std::format("Species_{}", species_idx);
            oss << generate_field_data(species_name, species_data);
        }
    }
    
    // Velocity vector for visualization
    std::vector<std::array<double, 3>> velocity_vectors;
    velocity_vectors.reserve(total_points);
    
    for (const auto& station : stations) {
        for (std::size_t eta_idx = 0; eta_idx < n_eta; ++eta_idx) {
            // Convert F and V to physical velocity components
            double u = (eta_idx < station.F.size()) ? station.F[eta_idx] : 0.0;  // Streamwise
            double v = (eta_idx < station.V.size()) ? station.V[eta_idx] : 0.0;  // Normal
            double w = 0.0;  // Spanwise
            
            velocity_vectors.push_back({u, v, w});
        }
    }
    oss << generate_vector_data("Velocity", velocity_vectors);
    
    oss << "      </PointData>\n";
    return oss.str();
}

auto VTKXMLWriter::generate_field_data(
    const std::string& name,
    const std::vector<double>& data,
    const std::string& type
) const -> std::string {
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    
    if (vtk_config_.data_mode == VTKConfig::DataMode::ASCII) {
        oss << "        <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">\n";
        oss << "          ";
        
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (i > 0 && i % 6 == 0) oss << "\n          ";
            oss << data[i] << " ";
        }
        oss << "\n        </DataArray>\n";
    } else {
        oss << "        <DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"binary\">\n";
        
        auto encoded_result = encode_binary_data(data);
        if (encoded_result) {
            oss << "          " << encoded_result.value() << "\n";
        }
        
        oss << "        </DataArray>\n";
    }
    
    return oss.str();
}

auto VTKXMLWriter::generate_vector_data(
    const std::string& name,
    const std::vector<std::array<double, 3>>& data
) const -> std::string {
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    
    if (vtk_config_.data_mode == VTKConfig::DataMode::ASCII) {
        oss << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        
        for (const auto& vec : data) {
            oss << "          " << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
        }
        
        oss << "        </DataArray>\n";
    } else {
        oss << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"binary\">\n";
        
        std::vector<double> flat_data;
        flat_data.reserve(data.size() * 3);
        for (const auto& vec : data) {
            flat_data.insert(flat_data.end(), vec.begin(), vec.end());
        }
        
        auto encoded_result = encode_binary_data(flat_data);
        if (encoded_result) {
            oss << "          " << encoded_result.value() << "\n";
        }
        
        oss << "        </DataArray>\n";
    }
    
    return oss.str();
}

auto VTKXMLWriter::encode_binary_data(
    const std::vector<double>& data
) const -> std::expected<std::string, OutputError> {
    
    if (data.empty()) {
        return std::string();
    }
    
    try {
        // Convert to bytes
        const auto byte_size = data.size() * sizeof(double);
        const auto* bytes = reinterpret_cast<const uint8_t*>(data.data());
        std::vector<uint8_t> byte_data(bytes, bytes + byte_size);
        
        // Compress if enabled
        if (vtk_config_.compress_binary) {
            auto compressed_result = compress_binary_data(byte_data);
            if (!compressed_result) {
                return std::unexpected(compressed_result.error());
            }
            byte_data = std::move(compressed_result.value());
        }
        
        // Base64 encode
        const std::string base64_chars = 
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        
        std::string encoded;
        int val = 0, valb = -6;
        for (uint8_t c : byte_data) {
            val = (val << 8) + c;
            valb += 8;
            while (valb >= 0) {
                encoded.push_back(base64_chars[(val >> valb) & 0x3F]);
                valb -= 6;
            }
        }
        if (valb > -6) {
            encoded.push_back(base64_chars[((val << 8) >> (valb + 8)) & 0x3F]);
        }
        while (encoded.size() % 4) {
            encoded.push_back('=');
        }
        
        return encoded;
        
    } catch (const std::exception& e) {
        return std::unexpected(OutputError(
            std::format("Binary encoding failed: {}", e.what())
        ));
    }
}

auto VTKXMLWriter::compress_binary_data(
    const std::vector<uint8_t>& data
) const -> std::expected<std::vector<uint8_t>, OutputError> {
    
    try {
        uLongf compressed_size = compressBound(data.size());
        std::vector<uint8_t> compressed_data(compressed_size);
        
        int result = compress(compressed_data.data(), &compressed_size, 
                             data.data(), data.size());
        
        if (result != Z_OK) {
            return std::unexpected(OutputError(
                std::format("Compression failed with code: {}", result)
            ));
        }
        
        compressed_data.resize(compressed_size);
        return compressed_data;
        
    } catch (const std::exception& e) {
        return std::unexpected(OutputError(
            std::format("Compression error: {}", e.what())
        ));
    }
}

// ParaView Collection Writer implementation
auto ParaViewCollectionWriter::add_time_step(
    double time,
    const std::filesystem::path& file_path,
    const std::string& part_name
) -> void {
    time_steps_.emplace_back(TimeStep{time, file_path, part_name});
}

auto ParaViewCollectionWriter::write_collection() const 
    -> std::expected<void, OutputError> {
    
    std::ofstream file(collection_path_);
    if (!file.is_open()) {
        return std::unexpected(FileWriteError(collection_path_, "Cannot open collection file"));
    }
    
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <Collection>\n";
    
    for (const auto& step : time_steps_) {
        auto relative_path = std::filesystem::relative(step.file_path, collection_path_.parent_path());
        file << "    <DataSet timestep=\"" << std::scientific << step.time 
             << "\" part=\"0\" file=\"" << relative_path.string() << "\"/>\n";
    }
    
    file << "  </Collection>\n";
    file << "</VTKFile>\n";
    
    return {};
}

// VTK utility functions
namespace vtk {

auto transform_coordinates(
    const std::vector<double>& eta,
    const std::vector<double>& xi,
    const std::vector<double>& y_physical,
    const std::vector<double>& x_physical
) -> std::vector<std::array<double, 3>> {
    
    std::vector<std::array<double, 3>> coords;
    const auto n_points = std::min({eta.size(), y_physical.size()});
    coords.reserve(n_points);
    
    for (std::size_t i = 0; i < n_points; ++i) {
        coords.push_back({
            (i < x_physical.size()) ? x_physical[i] : 0.0,  // x (streamwise)
            y_physical[i],                                   // y (normal)
            0.0                                             // z (spanwise)
        });
    }
    
    return coords;
}

auto create_structured_topology(
    int n_eta,
    int n_xi
) -> std::vector<std::vector<int>> {
    
    std::vector<std::vector<int>> cells;
    cells.reserve((n_eta - 1) * (n_xi - 1));
    
    for (int i = 0; i < n_xi - 1; ++i) {
        for (int j = 0; j < n_eta - 1; ++j) {
            // Quad cell indices (VTK ordering)
            std::vector<int> cell = {
                i * n_eta + j,           // (i, j)
                (i + 1) * n_eta + j,     // (i+1, j)
                (i + 1) * n_eta + j + 1, // (i+1, j+1)
                i * n_eta + j + 1        // (i, j+1)
            };
            cells.push_back(cell);
        }
    }
    
    return cells;
}

auto validate_vtk_file(const std::filesystem::path& file_path)
    -> std::expected<void, OutputError> {
    
    if (!std::filesystem::exists(file_path)) {
        return std::unexpected(OutputError(
            std::format("VTK file does not exist: {}", file_path.string())
        ));
    }
    
    std::ifstream file(file_path);
    if (!file.is_open()) {
        return std::unexpected(OutputError(
            std::format("Cannot open VTK file: {}", file_path.string())
        ));
    }
    
    std::string first_line;
    std::getline(file, first_line);
    
    if (first_line.find("<?xml") != std::string::npos || 
        first_line.find("# vtk DataFile") != std::string::npos) {
        return {};
    }
    
    return std::unexpected(OutputError("Invalid VTK file format"));
}

auto get_vtk_type_string(const std::type_info& type) -> std::string {
    if (type == typeid(float)) return "Float32";
    if (type == typeid(double)) return "Float64";
    if (type == typeid(int)) return "Int32";
    if (type == typeid(long)) return "Int64";
    if (type == typeid(unsigned int)) return "UInt32";
    if (type == typeid(unsigned long)) return "UInt64";
    return "Float64"; // Default
}

auto compute_cell_data(
    const std::vector<double>& point_data,
    int n_eta,
    int n_xi
) -> std::vector<double> {
    
    std::vector<double> cell_data;
    const int n_cells = (n_eta - 1) * (n_xi - 1);
    cell_data.reserve(n_cells);
    
    for (int i = 0; i < n_xi - 1; ++i) {
        for (int j = 0; j < n_eta - 1; ++j) {
            // Average of corner values
            double avg = 0.0;
            avg += point_data[i * n_eta + j];
            avg += point_data[(i + 1) * n_eta + j];
            avg += point_data[(i + 1) * n_eta + j + 1];
            avg += point_data[i * n_eta + j + 1];
            avg /= 4.0;
            
            cell_data.push_back(avg);
        }
    }
    
    return cell_data;
}

} // namespace vtk

} // namespace blast::io::output