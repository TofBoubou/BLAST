#pragma once
#include "output_writer.hpp"
#include <fstream>
#include <iomanip>
#include <map>

namespace blast::io::output {

// CSV-specific configuration
struct CSVConfig {
    char delimiter = ',';
    int precision = 8;
    bool include_headers = true;
    bool scientific_notation = false;
    std::string line_ending = "\n";
    
    // File organization
    bool separate_files_per_station = true;
    bool separate_files_per_variable = false;
    bool single_combined_file = false;
};

// CSV writer for simple data analysis
class CSVWriter : public FormatWriter {
private:
    CSVConfig csv_config_;
    
    // Core writing functions
    [[nodiscard]] auto write_single_file(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_station_files(
        const std::filesystem::path& base_path,
        const OutputDataset& dataset,
        const OutputConfig& config,
        ProgressCallback progress
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_variable_files(
        const std::filesystem::path& base_path,
        const OutputDataset& dataset,
        const OutputConfig& config
    ) const -> std::expected<void, OutputError>;
    
    // Station data writing
    [[nodiscard]] auto write_station_csv(
        const std::filesystem::path& file_path,
        const StationData& station,
        const OutputConfig& config,
        const SimulationMetadata& metadata,
        std::size_t station_index
    ) const -> std::expected<void, OutputError>;
    
    // Utility functions
    [[nodiscard]] auto format_value(double value) const -> std::string;
    
    [[nodiscard]] auto create_station_header(
        const OutputConfig& config,
        const SimulationMetadata& metadata
    ) const -> std::vector<std::string>;
    
    [[nodiscard]] auto write_metadata_file(
        const std::filesystem::path& file_path,
        const SimulationMetadata& metadata
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_variable_file(
        const std::filesystem::path& file_path,
        const std::vector<StationData>& stations,
        const std::string& variable_name,
        const OutputConfig& config
    ) const -> std::expected<void, OutputError>;
    
    [[nodiscard]] auto write_species_file(
        const std::filesystem::path& file_path,
        const std::vector<StationData>& stations,
        const SimulationMetadata& metadata,
        const OutputConfig& config
    ) const -> std::expected<void, OutputError>;

public:
    explicit CSVWriter(CSVConfig config = {}) : csv_config_(config) {}
    
    [[nodiscard]] auto write(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config,
        ProgressCallback progress = nullptr
    ) const -> std::expected<void, OutputError> override;
    
    [[nodiscard]] auto get_extension() const noexcept -> std::string_view override {
        return ".csv";
    }
    
    [[nodiscard]] auto supports_compression() const noexcept -> bool override {
        return false;
    }
    
    [[nodiscard]] auto supports_metadata() const noexcept -> bool override {
        return false;  // Limited metadata support via separate files
    }
    
    // CSV-specific configuration
    auto set_csv_config(CSVConfig config) noexcept -> void {
        csv_config_ = config;
    }
    
    [[nodiscard]] auto get_csv_config() const noexcept -> const CSVConfig& {
        return csv_config_;
    }
};

// Specialized CSV writers for specific use cases
class ProfileCSVWriter {
private:
    CSVConfig config_;
    
public:
    explicit ProfileCSVWriter(CSVConfig config = {}) : config_(config) {}
    
    // Write boundary layer profiles for a single station
    [[nodiscard]] auto write_profiles(
        const std::filesystem::path& file_path,
        const StationData& station,
        const std::vector<std::string>& species_names,
        std::size_t station_index = 0
    ) const -> std::expected<void, OutputError>;
    
    // Write wall properties along the body
    [[nodiscard]] auto write_wall_properties(
        const std::filesystem::path& file_path,
        const OutputDataset::WallData& wall_data,
        const std::vector<double>& xi_coordinates
    ) const -> std::expected<void, OutputError>;
    
    // Write integrated quantities
    [[nodiscard]] auto write_integrated_quantities(
        const std::filesystem::path& file_path,
        const OutputDataset::IntegratedQuantities& quantities,
        const std::vector<double>& xi_coordinates
    ) const -> std::expected<void, OutputError>;
};

// CSV reader for post-processing
class CSVReader {
private:
    CSVConfig config_;
    
    [[nodiscard]] auto parse_line(
        const std::string& line,
        std::vector<std::string>& tokens
    ) const -> void;
    
    [[nodiscard]] auto parse_header(
        const std::string& header_line
    ) const -> std::vector<std::string>;

public:
    explicit CSVReader(CSVConfig config = {}) : config_(config) {}
    
    // Read station data from CSV
    [[nodiscard]] auto read_station_data(
        const std::filesystem::path& file_path
    ) const -> std::expected<std::map<std::string, std::vector<double>>, OutputError>;
    
    // Read time series data
    [[nodiscard]] auto read_time_series(
        const std::filesystem::path& file_path,
        const std::string& variable_name
    ) const -> std::expected<std::pair<std::vector<double>, std::vector<double>>, OutputError>;
    
    // Get available variables from CSV file
    [[nodiscard]] auto get_variable_names(
        const std::filesystem::path& file_path
    ) const -> std::expected<std::vector<std::string>, OutputError>;
};

// Convenience functions for CSV operations
namespace csv {
    
    // Write complete solution as CSV series
    [[nodiscard]] auto write_csv_series(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& output_directory
    ) -> std::expected<std::vector<std::filesystem::path>, OutputError>;
    
    // Write only profiles (most common use case)
    [[nodiscard]] auto write_profiles_only(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& output_directory
    ) -> std::expected<std::vector<std::filesystem::path>, OutputError>;
    
    // Write wall properties
    [[nodiscard]] auto write_wall_data(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& file_path
    ) -> std::expected<void, OutputError>;
    
    // Validate CSV file format
    [[nodiscard]] auto validate_csv_file(
        const std::filesystem::path& file_path,
        char expected_delimiter = ','
    ) -> std::expected<void, OutputError>;
    
    // Convert CSV to other formats
    [[nodiscard]] auto csv_to_hdf5(
        const std::filesystem::path& csv_directory,
        const std::filesystem::path& hdf5_file
    ) -> std::expected<void, OutputError>;

} // namespace csv

} // namespace blast::io::output