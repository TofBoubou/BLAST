#include "blast/io/output/hdf5_writer.hpp"
#include <cstring>
#include <iomanip>
#include <sstream>
#include <vector>

namespace blast::io::output {

// HDF5Writer implementation
auto HDF5Writer::write(const std::filesystem::path& file_path, const OutputDataset& dataset, const OutputConfig& config,
                       ProgressCallback progress) const -> std::expected<void, OutputError> {

  try {
    if (progress)
      progress(0.0, "Creating HDF5 file");

    // Create HDF5 file
    auto file_result = create_file(file_path);
    if (!file_result) {
      return std::unexpected(file_result.error());
    }
    auto file = std::move(file_result.value());

    if (progress)
      progress(0.1, "Writing metadata");

    // Write metadata
    if (config.save_metadata) {
      if (auto meta_result = write_metadata(file, dataset.metadata); !meta_result) {
        return std::unexpected(meta_result.error());
      }
    }

    if (progress)
      progress(0.3, "Writing station data");

    // Write station data
    if (auto stations_result = write_stations(file, dataset.stations, config, progress); !stations_result) {
      return std::unexpected(stations_result.error());
    }

    if (progress)
      progress(0.8, "Writing wall and integrated data");

    if (progress)
      progress(1.0, "HDF5 write complete");

    return {};

  } catch (const std::exception& e) {
    return std::unexpected(OutputError(std::format("HDF5 write failed: {}", e.what())));
  }
}

auto HDF5Writer::create_file(const std::filesystem::path& file_path) const -> std::expected<FileHandle, OutputError> {

  // Create file access property list
  auto fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (fapl < 0) {
    return std::unexpected(OutputError("Failed to create file access property list"));
  }

  // Set close degree (for proper cleanup)
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

  // Create the file
  auto file_id = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);

  if (file_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create HDF5 file: {}", file_path.string())));
  }

  try {
    return FileHandle(file_id);
  } catch (const std::exception& e) {
    H5Fclose(file_id);
    return std::unexpected(OutputError(e.what()));
  }
}

auto HDF5Writer::write_metadata(FileHandle& file,
                                const SimulationMetadata& metadata) const -> std::expected<void, OutputError> {

  // Create metadata group
  auto metadata_group_result = create_group(file, "metadata");
  if (!metadata_group_result) {
    return std::unexpected(metadata_group_result.error());
  }
  auto metadata_group = std::move(metadata_group_result.value());

  // Write basic metadata
  if (auto result = write_string(metadata_group, "blast_version", metadata.blast_version); !result) {
    return std::unexpected(result.error());
  }

  // Write creation time
  auto time_t = std::chrono::system_clock::to_time_t(metadata.creation_time);
  auto tm = *std::gmtime(&time_t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
  if (auto result = write_string(metadata_group, "creation_time", oss.str()); !result) {
    return std::unexpected(result.error());
  }

  // Write grid information
  auto grid_group_result = create_group(metadata_group, "grid");
  if (!grid_group_result) {
    return std::unexpected(grid_group_result.error());
  }
  auto grid_group = std::move(grid_group_result.value());

  if (auto result = write_scalar(grid_group, "n_eta", metadata.grid.n_eta); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_scalar(grid_group, "eta_max", metadata.grid.eta_max); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_scalar(grid_group, "d_eta", metadata.grid.d_eta); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_vector(grid_group, "xi_coordinates", metadata.grid.xi_coordinates); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_vector(grid_group, "eta_coordinates", metadata.grid.eta_coordinates); !result) {
    return std::unexpected(result.error());
  }

  // Write mixture information
  auto mixture_group_result = create_group(metadata_group, "mixture");
  if (!mixture_group_result) {
    return std::unexpected(mixture_group_result.error());
  }
  auto mixture_group = std::move(mixture_group_result.value());

  if (auto result = write_string_array(mixture_group, "species_names", metadata.mixture.species_names); !result) {
    return std::unexpected(result.error());
  }
  if (auto result =
          write_vector(mixture_group, "species_molecular_weights", metadata.mixture.species_molecular_weights);
      !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_vector(mixture_group, "species_charges", metadata.mixture.species_charges); !result) {
    return std::unexpected(result.error());
  }

  // Write convergence information
  auto conv_group_result = create_group(metadata_group, "convergence");
  if (!conv_group_result) {
    return std::unexpected(conv_group_result.error());
  }
  auto conv_group = std::move(conv_group_result.value());

  if (auto result = write_scalar(conv_group, "converged", metadata.convergence.converged ? 1.0 : 0.0); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_scalar(conv_group, "total_iterations", metadata.convergence.total_iterations); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_scalar(conv_group, "final_residual", metadata.convergence.final_residual); !result) {
    return std::unexpected(result.error());
  }

  return {};
}

auto HDF5Writer::write_stations(FileHandle& file, const std::vector<StationData>& stations, const OutputConfig& config,
                                ProgressCallback progress) const -> std::expected<void, OutputError> {

  // Create stations group
  auto stations_group_result = create_group(file, "stations");
  if (!stations_group_result) {
    return std::unexpected(stations_group_result.error());
  }
  auto stations_group = std::move(stations_group_result.value());

  // Write each station
  for (std::size_t i = 0; i < stations.size(); ++i) {
    if (progress) {
      double station_progress = 0.3 + 0.5 * (static_cast<double>(i) / stations.size());
      progress(station_progress, std::format("Writing station {}", i));
    }

    // Create station group
    auto station_name = std::format("station_{:03d}", i);
    auto station_group_result = create_group(stations_group, station_name);
    if (!station_group_result) {
      return std::unexpected(station_group_result.error());
    }
    auto station_group = std::move(station_group_result.value());

    // Write station data
    if (auto result = write_station_data(station_group, stations[i], config); !result) {
      return std::unexpected(result.error());
    }
  }

  return {};
}

auto HDF5Writer::write_station_data(GroupHandle& station_group, const StationData& station,
                                    const OutputConfig& config) const -> std::expected<void, OutputError> {

  // Coordinates
  if (auto result = write_scalar(station_group, "xi", station.xi, "", "Streamwise coordinate"); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_scalar(station_group, "x_physical", station.x_physical, "m", "Physical x coordinate");
      !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_vector(station_group, "eta", station.eta, "", "Similarity coordinate"); !result) {
    return std::unexpected(result.error());
  }
  if (auto result = write_vector(station_group, "y_physical", station.y_physical, "m", "Physical y coordinate");
      !result) {
    return std::unexpected(result.error());
  }

  // Flow variables
  if (config.variables.flow_variables) {
    if (auto result = write_vector(station_group, "F", station.F, "", "Dimensionless stream function"); !result) {
      return std::unexpected(result.error());
    }
    if (auto result = write_vector(station_group, "g", station.g, "", "Dimensionless enthalpy"); !result) {
      return std::unexpected(result.error());
    }
    if (auto result = write_vector(station_group, "V", station.V, "", "Velocity field"); !result) {
      return std::unexpected(result.error());
    }
    if (auto result = write_vector(station_group, "temperature", station.temperature, "K", "Temperature"); !result) {
      return std::unexpected(result.error());
    }
  }

  // Species concentrations
  if (config.variables.species_concentrations) {
    if (auto result = write_matrix(station_group, "species_concentrations", station.species_concentrations, "",
                                   "Species mass fractions");
        !result) {
      return std::unexpected(result.error());
    }
  }

  return {};
}

// Utility function implementations
auto HDF5Writer::create_group(hid_t parent, const std::string& name) const -> std::expected<GroupHandle, OutputError> {

  auto group_id = H5Gcreate2(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (group_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create group '{}'", name)));
  }

  try {
    return GroupHandle(group_id);
  } catch (const std::exception& e) {
    H5Gclose(group_id);
    return std::unexpected(OutputError(e.what()));
  }
}

auto HDF5Writer::write_vector(hid_t parent, const std::string& name, const std::vector<double>& data,
                              const std::string& units,
                              const std::string& description) const -> std::expected<void, OutputError> {

  if (data.empty()) {
    return {}; // Skip empty datasets
  }

  // Create dataspace
  hsize_t dims = data.size();
  auto space_id = H5Screate_simple(1, &dims, nullptr);
  if (space_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create dataspace for '{}'", name)));
  }
  DataspaceHandle space(space_id);

  // Create dataset properties - use chunked properties if compression is
  // enabled
  auto prop_result =
      (hdf5_config_.compression_level > 0) ? create_chunked_properties(data.size()) : create_dataset_properties();

  if (!prop_result) {
    return std::unexpected(prop_result.error());
  }
  auto props = std::move(prop_result.value());

  auto dataset_id = H5Dcreate2(parent, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, props, H5P_DEFAULT);
  if (dataset_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create dataset '{}'", name)));
  }
  DatasetHandle dataset(dataset_id);

  // Write data
  auto status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  if (status < 0) {
    return std::unexpected(OutputError(std::format("Failed to write data for '{}'", name)));
  }

  // Write attributes
  if (!units.empty()) {
    if (auto result = write_string(dataset, "units", units); !result) {
      return std::unexpected(result.error());
    }
  }
  if (!description.empty()) {
    if (auto result = write_string(dataset, "description", description); !result) {
      return std::unexpected(result.error());
    }
  }

  return {};
}

auto HDF5Writer::write_matrix(hid_t parent, const std::string& name, const core::Matrix<double>& data,
                              const std::string& units,
                              const std::string& description) const -> std::expected<void, OutputError> {

  if (data.rows() == 0 || data.cols() == 0) {
    return {}; // Skip empty matrices
  }

  // Create dataspace
  hsize_t dims[2] = {static_cast<hsize_t>(data.rows()), static_cast<hsize_t>(data.cols())};
  auto space_id = H5Screate_simple(2, dims, nullptr);
  if (space_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create dataspace for matrix '{}'", name)));
  }
  DataspaceHandle space(space_id);

  // Create dataset properties - use chunked properties if compression is
  // enabled
  auto prop_result = (hdf5_config_.compression_level > 0) ? create_chunked_properties_2d(data.rows(), data.cols())
                                                          : create_dataset_properties();

  if (!prop_result) {
    return std::unexpected(prop_result.error());
  }
  auto props = std::move(prop_result.value());

  auto dataset_id = H5Dcreate2(parent, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, props, H5P_DEFAULT);
  if (dataset_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create dataset '{}'", name)));
  }
  DatasetHandle dataset(dataset_id);

  // Write data (Eigen matrices are column-major, HDF5 expects row-major)
  auto& eigen_data = data.eigen();
  std::vector<double> row_major_data(data.rows() * data.cols());
  for (std::size_t i = 0; i < data.rows(); ++i) {
    for (std::size_t j = 0; j < data.cols(); ++j) {
      row_major_data[i * data.cols() + j] = eigen_data(i, j);
    }
  }

  auto status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, row_major_data.data());
  if (status < 0) {
    return std::unexpected(OutputError(std::format("Failed to write matrix data for '{}'", name)));
  }

  // Write attributes
  if (!units.empty()) {
    if (auto result = write_string(dataset, "units", units); !result) {
      return std::unexpected(result.error());
    }
  }
  if (!description.empty()) {
    if (auto result = write_string(dataset, "description", description); !result) {
      return std::unexpected(result.error());
    }
  }

  return {};
}

auto HDF5Writer::create_chunked_properties_2d(std::size_t rows,
                                              std::size_t cols) const -> std::expected<PropertyHandle, OutputError> {

  auto plist_id = H5Pcreate(H5P_DATASET_CREATE);
  if (plist_id < 0) {
    return std::unexpected(OutputError("Failed to create dataset property list"));
  }
  PropertyHandle props(plist_id);

  // Set chunking for 2D data - chunk sizes must be <= data sizes
  hsize_t chunk_dims[2];
  chunk_dims[0] = std::min(rows, hdf5_config_.chunk_size);
  chunk_dims[1] = std::min(cols, hdf5_config_.chunk_size);

  // Ensure minimum chunk size of 1 in each dimension
  if (chunk_dims[0] == 0)
    chunk_dims[0] = 1;
  if (chunk_dims[1] == 0)
    chunk_dims[1] = 1;

  if (H5Pset_chunk(props, 2, chunk_dims) < 0) {
    return std::unexpected(OutputError("Failed to set 2D chunking"));
  }

  // Set compression if enabled
  if (hdf5_config_.compression_level > 0) {
    if (hdf5_config_.use_shuffle_filter) {
      H5Pset_shuffle(props);
    }

    H5Pset_deflate(props, hdf5_config_.compression_level);

    if (hdf5_config_.use_fletcher32) {
      H5Pset_fletcher32(props);
    }
  }

  return props;
}

auto HDF5Writer::write_scalar(hid_t parent, const std::string& name, double value, const std::string& units,
                              const std::string& description) const -> std::expected<void, OutputError> {

  // Create scalar dataspace
  auto space_id = H5Screate(H5S_SCALAR);
  if (space_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create scalar dataspace for '{}'", name)));
  }
  DataspaceHandle space(space_id);

  // Create dataset
  auto dataset_id = H5Dcreate2(parent, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create scalar dataset '{}'", name)));
  }
  DatasetHandle dataset(dataset_id);

  // Write data
  auto status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
  if (status < 0) {
    return std::unexpected(OutputError(std::format("Failed to write scalar value for '{}'", name)));
  }

  // Write attributes
  if (!units.empty()) {
    if (auto result = write_string(dataset, "units", units); !result) {
      return std::unexpected(result.error());
    }
  }
  if (!description.empty()) {
    if (auto result = write_string(dataset, "description", description); !result) {
      return std::unexpected(result.error());
    }
  }

  return {};
}

auto HDF5Writer::write_string(hid_t parent, const std::string& name,
                              const std::string& value) const -> std::expected<void, OutputError> {

  // Create string type
  auto str_type = H5Tcopy(H5T_C_S1);
  if (str_type < 0) {
    return std::unexpected(OutputError("Failed to create string type"));
  }
  TypeHandle string_type(str_type);

  H5Tset_size(string_type, value.length());
  H5Tset_strpad(string_type, H5T_STR_NULLTERM);

  // Create scalar dataspace
  auto space_id = H5Screate(H5S_SCALAR);
  if (space_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create dataspace for string '{}'", name)));
  }
  DataspaceHandle space(space_id);

  // Check if this is an attribute or dataset
  H5I_type_t obj_type = H5Iget_type(parent);
  if (obj_type == H5I_DATASET || obj_type == H5I_GROUP) {
    // Write as attribute
    auto attr_id = H5Acreate2(parent, name.c_str(), string_type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
      return std::unexpected(OutputError(std::format("Failed to create string attribute '{}'", name)));
    }

    auto status = H5Awrite(attr_id, string_type, value.c_str());
    H5Aclose(attr_id);

    if (status < 0) {
      return std::unexpected(OutputError(std::format("Failed to write string attribute '{}'", name)));
    }
  }

  return {};
}

auto HDF5Writer::write_string_array(hid_t parent, const std::string& name,
                                    const std::vector<std::string>& values) const -> std::expected<void, OutputError> {

  if (values.empty()) {
    return {}; // Skip empty arrays
  }

  // Find maximum string length
  std::size_t max_len = 0;
  for (const auto& str : values) {
    max_len = std::max(max_len, str.length());
  }
  ++max_len; // For null terminator

  // Create string type
  auto str_type = H5Tcopy(H5T_C_S1);
  if (str_type < 0) {
    return std::unexpected(OutputError("Failed to create string type"));
  }
  TypeHandle string_type(str_type);

  H5Tset_size(string_type, max_len);
  H5Tset_strpad(string_type, H5T_STR_NULLTERM);

  // Create dataspace
  hsize_t dims = values.size();
  auto space_id = H5Screate_simple(1, &dims, nullptr);
  if (space_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create dataspace for string array '{}'", name)));
  }
  DataspaceHandle space(space_id);

  // Create dataset
  auto dataset_id = H5Dcreate2(parent, name.c_str(), string_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id < 0) {
    return std::unexpected(OutputError(std::format("Failed to create string array dataset '{}'", name)));
  }
  DatasetHandle dataset(dataset_id);

  // Prepare data buffer
  std::vector<char> buffer(values.size() * max_len, '\0');
  for (std::size_t i = 0; i < values.size(); ++i) {
    std::strncpy(&buffer[i * max_len], values[i].c_str(), max_len - 1);
  }

  // Write data
  auto status = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
  if (status < 0) {
    return std::unexpected(OutputError(std::format("Failed to write string array data for '{}'", name)));
  }

  return {};
}

auto HDF5Writer::create_dataset_properties() const -> std::expected<PropertyHandle, OutputError> {

  auto plist_id = H5Pcreate(H5P_DATASET_CREATE);
  if (plist_id < 0) {
    return std::unexpected(OutputError("Failed to create dataset property list"));
  }
  PropertyHandle props(plist_id);

  // Set compression if enabled - but no chunking here since we don't know the
  // size
  if (hdf5_config_.compression_level > 0) {
    if (hdf5_config_.use_shuffle_filter) {
      H5Pset_shuffle(props);
    }

    H5Pset_deflate(props, hdf5_config_.compression_level);

    if (hdf5_config_.use_fletcher32) {
      H5Pset_fletcher32(props);
    }
  }

  return props;
}

auto HDF5Writer::create_chunked_properties(std::size_t size) const -> std::expected<PropertyHandle, OutputError> {

  auto plist_id = H5Pcreate(H5P_DATASET_CREATE);
  if (plist_id < 0) {
    return std::unexpected(OutputError("Failed to create dataset property list"));
  }
  PropertyHandle props(plist_id);

  // Set chunking - chunk size must be <= data size
  hsize_t chunk_size = std::min(size, hdf5_config_.chunk_size);
  if (chunk_size == 0) {
    chunk_size = 1; // Minimum chunk size
  }

  if (H5Pset_chunk(props, 1, &chunk_size) < 0) {
    return std::unexpected(OutputError("Failed to set chunking"));
  }

  // Set compression if enabled
  if (hdf5_config_.compression_level > 0) {
    if (hdf5_config_.use_shuffle_filter) {
      H5Pset_shuffle(props);
    }

    H5Pset_deflate(props, hdf5_config_.compression_level);

    if (hdf5_config_.use_fletcher32) {
      H5Pset_fletcher32(props);
    }
  }

  return props;
}

auto HDF5Writer::get_hdf5_error() const -> std::string {
  // Get the current error stack
  std::string error_msg = "HDF5 error occurred";

  // TODO: Implement detailed HDF5 error extraction
  // This would involve walking the HDF5 error stack

  return error_msg;
}

// HDF5 convenience functions
namespace hdf5 {

auto initialize() -> std::expected<void, OutputError> {
  if (H5open() < 0) {
    return std::unexpected(OutputError("Failed to initialize HDF5 library"));
  }
  return {};
}

auto finalize() -> void {
  H5close();
}

auto check_version() -> std::expected<std::string, OutputError> {
  unsigned majnum, minnum, relnum;
  if (H5get_libversion(&majnum, &minnum, &relnum) < 0) {
    return std::unexpected(OutputError("Failed to get HDF5 version"));
  }

  return std::format("{}.{}.{}", majnum, minnum, relnum);
}

auto validate_file(const std::filesystem::path& file_path) -> std::expected<void, OutputError> {

  if (!std::filesystem::exists(file_path)) {
    return std::unexpected(OutputError(std::format("File does not exist: {}", file_path.string())));
  }

  // Check if it's a valid HDF5 file
  auto result = H5Fis_hdf5(file_path.c_str());
  if (result <= 0) {
    return std::unexpected(OutputError(std::format("Not a valid HDF5 file: {}", file_path.string())));
  }

  return {};
}

} // namespace hdf5

} // namespace blast::io::output