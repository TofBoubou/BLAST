#pragma once
#include "output_writer.hpp"
#include <hdf5.h>
#include <map>
#include <memory>
#include <string_view>

namespace blast::io::output {

// HDF5-specific configuration
struct HDF5Config {
  int compression_level = 6;      // 0-9, higher = better compression
  bool use_chunking = true;       // Enable chunked layout
  bool use_shuffle_filter = true; // Reorder bytes for better compression
  bool use_fletcher32 = false;    // Checksum filter
  std::size_t chunk_size = 1024;  // Chunk size for datasets
  bool store_as_double = true;    // Use double precision by default
};

// RAII wrapper for HDF5 handles
template <typename HandleType, auto CloseFunc> class HDF5Handle {
private:
  HandleType handle_;

public:
  explicit HDF5Handle(HandleType handle) : handle_(handle) {
    if (handle_ < 0) {
      throw OutputError("Invalid HDF5 handle");
    }
  }

  ~HDF5Handle() {
    if (handle_ >= 0) {
      CloseFunc(handle_);
    }
  }

  // Move semantics only
  HDF5Handle(HDF5Handle&& other) noexcept : handle_(other.handle_) { other.handle_ = -1; }

  HDF5Handle& operator=(HDF5Handle&& other) noexcept {
    if (this != &other) {
      if (handle_ >= 0) {
        CloseFunc(handle_);
      }
      handle_ = other.handle_;
      other.handle_ = -1;
    }
    return *this;
  }

  // Delete copy operations
  HDF5Handle(const HDF5Handle&) = delete;
  HDF5Handle& operator=(const HDF5Handle&) = delete;

  [[nodiscard]] auto get() const noexcept -> HandleType { return handle_; }
  [[nodiscard]] auto valid() const noexcept -> bool { return handle_ >= 0; }

  // Implicit conversion for C API
  operator HandleType() const noexcept { return handle_; }
};

// Type aliases for HDF5 handles
using FileHandle = HDF5Handle<hid_t, H5Fclose>;
using GroupHandle = HDF5Handle<hid_t, H5Gclose>;
using DatasetHandle = HDF5Handle<hid_t, H5Dclose>;
using DataspaceHandle = HDF5Handle<hid_t, H5Sclose>;
using PropertyHandle = HDF5Handle<hid_t, H5Pclose>;
using TypeHandle = HDF5Handle<hid_t, H5Tclose>;

// HDF5 implementation of FormatWriter
class HDF5Writer : public FormatWriter {
private:
  HDF5Config hdf5_config_;

  // Core writing functions
  [[nodiscard]] auto
  create_file(const std::filesystem::path& file_path) const -> std::expected<FileHandle, OutputError>;

  [[nodiscard]] auto write_metadata(FileHandle& file,
                                    const SimulationMetadata& metadata) const -> std::expected<void, OutputError>;

  [[nodiscard]] auto write_stations(FileHandle& file, const std::vector<StationData>& stations,
                                    const OutputConfig& config,
                                    ProgressCallback progress) const -> std::expected<void, OutputError>;

  [[nodiscard]] auto write_station_data(GroupHandle& station_group, const StationData& station,
                                        const OutputConfig& config) const -> std::expected<void, OutputError>;

  // Utility functions for HDF5 operations
  [[nodiscard]] auto create_group(hid_t parent,
                                  const std::string& name) const -> std::expected<GroupHandle, OutputError>;

  [[nodiscard]] auto write_vector(hid_t parent, const std::string& name, const std::vector<double>& data,
                                  const std::string& units = "",
                                  const std::string& description = "") const -> std::expected<void, OutputError>;

  [[nodiscard]] auto write_matrix(hid_t parent, const std::string& name, const core::Matrix<double>& data,
                                  const std::string& units = "",
                                  const std::string& description = "") const -> std::expected<void, OutputError>;

  [[nodiscard]] auto create_chunked_properties_2d(std::size_t rows,
                                                  std::size_t cols) const -> std::expected<PropertyHandle, OutputError>;

  [[nodiscard]] auto write_scalar(hid_t parent, const std::string& name, double value, const std::string& units = "",
                                  const std::string& description = "") const -> std::expected<void, OutputError>;

  [[nodiscard]] auto write_string(hid_t parent, const std::string& name,
                                  const std::string& value) const -> std::expected<void, OutputError>;

  [[nodiscard]] auto
  write_string_array(hid_t parent, const std::string& name,
                     const std::vector<std::string>& values) const -> std::expected<void, OutputError>;

  // Property list creation
  [[nodiscard]] auto create_dataset_properties() const -> std::expected<PropertyHandle, OutputError>;

  [[nodiscard]] auto create_chunked_properties(std::size_t size) const -> std::expected<PropertyHandle, OutputError>;

  // Error handling
  [[nodiscard]] auto get_hdf5_error() const -> std::string;

public:
  explicit HDF5Writer(HDF5Config config = {}) : hdf5_config_(config) {}

  [[nodiscard]] auto write(const std::filesystem::path& file_path, const OutputDataset& dataset,
                           const OutputConfig& config,
                           ProgressCallback progress = nullptr) const -> std::expected<void, OutputError> override;

  [[nodiscard]] auto get_extension() const noexcept -> std::string_view override { return ".h5"; }

  [[nodiscard]] auto supports_compression() const noexcept -> bool override { return true; }

  [[nodiscard]] auto supports_metadata() const noexcept -> bool override { return true; }

  // HDF5-specific configuration
  auto set_hdf5_config(HDF5Config config) noexcept -> void { hdf5_config_ = config; }

  [[nodiscard]] auto get_hdf5_config() const noexcept -> const HDF5Config& { return hdf5_config_; }
};

// HDF5 reader for post-processing analysis
class HDF5Reader {
private:
  FileHandle file_;

public:
  explicit HDF5Reader(const std::filesystem::path& file_path);

  // Read complete dataset
  [[nodiscard]] auto read_dataset() const -> std::expected<OutputDataset, OutputError>;

  // Read specific components
  [[nodiscard]] auto read_metadata() const -> std::expected<SimulationMetadata, OutputError>;

  [[nodiscard]] auto read_station(int station_index) const -> std::expected<StationData, OutputError>;

  [[nodiscard]] auto read_station_variable(int station_index, const std::string& variable_name) const
      -> std::expected<std::vector<double>, OutputError>;

  [[nodiscard]] auto get_station_count() const -> std::expected<int, OutputError>;

  [[nodiscard]] auto
  get_variable_names(int station_index = 0) const -> std::expected<std::vector<std::string>, OutputError>;

  // Utility functions
  [[nodiscard]] auto is_valid_file(const std::filesystem::path& file_path) const noexcept -> bool;

  [[nodiscard]] auto get_file_info() const -> std::expected<std::map<std::string, std::string>, OutputError>;
};

// Convenience functions for HDF5
namespace hdf5 {

// Initialize HDF5 library (call once at program start)
auto initialize() -> std::expected<void, OutputError>;

// Cleanup HDF5 library (call at program end)
auto finalize() -> void;

// Check HDF5 version compatibility
[[nodiscard]] auto check_version() -> std::expected<std::string, OutputError>;

// Validate HDF5 file
[[nodiscard]] auto validate_file(const std::filesystem::path& file_path) -> std::expected<void, OutputError>;

// Get dataset information without loading data
[[nodiscard]] auto get_dataset_info(const std::filesystem::path& file_path)
    -> std::expected<std::map<std::string, std::variant<int, double, std::string>>, OutputError>;

} // namespace hdf5

} // namespace blast::io::output