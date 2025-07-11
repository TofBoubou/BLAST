#pragma once
#include "output_writer.hpp"
#include <sstream>

namespace blast::io::output {

// VTK-specific configuration
struct VTKConfig {
    enum class Format { XML, Legacy };
    enum class DataMode { ASCII, Binary, Appended };
    
    Format format = Format::XML;
    DataMode data_mode = DataMode::Binary;
    bool compress_binary = true;
    int precision = 6;              // For ASCII output
    bool write_time_series = true;  // Create .pvd collection file
    bool include_ghost_cells = false;
    
    // For structured grids
    bool write_as_structured = true;  // vs unstructured
};

// VTK XML writer for modern ParaView compatibility
class VTKXMLWriter : public FormatWriter {
private:
    VTKConfig vtk_config_;
    
    // Core writing functions
    [[nodiscard]] auto write_structured_grid(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config,
        ProgressCallback progress
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_time_series(
        const std::filesystem::path& base_path,
        const OutputDataset& dataset,
        const OutputConfig& config
    ) const -> std::expected<void, OutputError>;
    
    // XML generation helpers
    [[nodiscard]] auto generate_xml_header() const -> std::string;
    
    [[nodiscard]] auto generate_grid_points(
        const std::vector<StationData>& stations
    ) const -> std::expected<std::string, OutputError>;
    
    [[nodiscard]] auto generate_point_data(
        const std::vector<StationData>& stations,
        const OutputConfig& config
    ) const -> std::expected<std::string, OutputError>;
    
    [[nodiscard]] auto generate_field_data(
        const std::string& name,
        const std::vector<double>& data,
        const std::string& type = "Float64"
    ) const -> std::string;
    
    [[nodiscard]] auto generate_vector_data(
        const std::string& name,
        const std::vector<std::array<double, 3>>& data
    ) const -> std::string;
    
    // Binary data encoding
    [[nodiscard]] auto encode_binary_data(
        const std::vector<double>& data
    ) const -> std::expected<std::string, OutputError>;
    
    [[nodiscard]] auto compress_binary_data(
        const std::vector<uint8_t>& data
    ) const -> std::expected<std::vector<uint8_t>, OutputError>;
    
    // Coordinate transformation
    [[nodiscard]] auto compute_physical_coordinates(
        const std::vector<StationData>& stations
    ) const -> std::expected<std::vector<std::array<double, 3>>, OutputError>;

public:
    explicit VTKXMLWriter(VTKConfig config = {}) : vtk_config_(config) {}
    
    [[nodiscard]] auto write(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config,
        ProgressCallback progress = nullptr
    ) const -> std::expected<void, OutputError> override;
    
    [[nodiscard]] auto get_extension() const noexcept -> std::string_view override {
        return vtk_config_.format == VTKConfig::Format::XML ? ".vts" : ".vtk";
    }
    
    [[nodiscard]] auto supports_compression() const noexcept -> bool override {
        return vtk_config_.format == VTKConfig::Format::XML;
    }
    
    [[nodiscard]] auto supports_metadata() const noexcept -> bool override {
        return true;
    }
    
    // VTK-specific configuration
    auto set_vtk_config(VTKConfig config) noexcept -> void {
        vtk_config_ = config;
    }
    
    [[nodiscard]] auto get_vtk_config() const noexcept -> const VTKConfig& {
        return vtk_config_;
    }
};

// VTK Legacy writer for compatibility with older tools
class VTKLegacyWriter : public FormatWriter {
private:
    VTKConfig vtk_config_;
    
    [[nodiscard]] auto write_legacy_format(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_header(std::ostream& os) const -> void;
    
    [[nodiscard]] auto write_coordinates(
        std::ostream& os,
        const std::vector<StationData>& stations
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_scalar_field(
        std::ostream& os,
        const std::string& name,
        const std::vector<double>& data
    ) const -> void;
    
    [[nodiscard]] auto write_vector_field(
        std::ostream& os,
        const std::string& name,
        const std::vector<std::array<double, 3>>& data
    ) const -> void;

public:
    explicit VTKLegacyWriter(VTKConfig config = {}) : vtk_config_(config) {
        vtk_config_.format = VTKConfig::Format::Legacy;
    }
    
    [[nodiscard]] auto write(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config,
        ProgressCallback progress = nullptr
    ) const -> std::expected<void, OutputError> override;
    
    [[nodiscard]] auto get_extension() const noexcept -> std::string_view override {
        return ".vtk";
    }
    
    [[nodiscard]] auto supports_compression() const noexcept -> bool override {
        return false;
    }
    
    [[nodiscard]] auto supports_metadata() const noexcept -> bool override {
        return false;  // Limited metadata support in legacy format
    }
};

// ParaView collection file writer (.pvd)
class ParaViewCollectionWriter {
private:
    struct TimeStep {
        double time;
        std::filesystem::path file_path;
        std::string part_name;
    };
    
    std::vector<TimeStep> time_steps_;
    std::filesystem::path collection_path_;

public:
    explicit ParaViewCollectionWriter(const std::filesystem::path& collection_path)
        : collection_path_(collection_path) {}
    
    auto add_time_step(
        double time,
        const std::filesystem::path& file_path,
        const std::string& part_name = "boundary_layer"
    ) -> void;
    
    [[nodiscard]] auto write_collection() const 
        -> std::expected<void, OutputError>;
    
    auto clear() -> void { time_steps_.clear(); }
    
    [[nodiscard]] auto get_time_step_count() const noexcept -> std::size_t {
        return time_steps_.size();
    }
};

// Utility functions for VTK
namespace vtk {
    
    // Convert BLAST coordinate system to VTK/ParaView conventions
    [[nodiscard]] auto transform_coordinates(
        const std::vector<double>& eta,
        const std::vector<double>& xi,
        const std::vector<double>& y_physical,
        const std::vector<double>& x_physical
    ) -> std::vector<std::array<double, 3>>;
    
    // Create structured grid topology
    [[nodiscard]] auto create_structured_topology(
        int n_eta,
        int n_xi
    ) -> std::vector<std::vector<int>>;
    
    // Validate VTK file
    [[nodiscard]] auto validate_vtk_file(const std::filesystem::path& file_path)
        -> std::expected<void, OutputError>;
    
    // Get VTK data type string
    [[nodiscard]] auto get_vtk_type_string(const std::type_info& type) -> std::string;
    
    // Compute cell-centered data from point data
    [[nodiscard]] auto compute_cell_data(
        const std::vector<double>& point_data,
        int n_eta,
        int n_xi
    ) -> std::vector<double>;
    
    // Generate streamlines for visualization
    [[nodiscard]] auto generate_streamlines(
        const std::vector<StationData>& stations,
        int n_streamlines = 10
    ) -> std::expected<std::vector<std::vector<std::array<double, 3>>>, OutputError>;

} // namespace vtk

} // namespace blast::io::output