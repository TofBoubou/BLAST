#pragma once
#include "output_types.hpp"
#include "../../boundary_layer/solver/boundary_layer_solver.hpp"
#include "../../thermophysics/mixture_interface.hpp"
#include <memory>
#include <expected>

namespace blast::io::output {

// Abstract base class for format-specific writers
class FormatWriter {
public:
    virtual ~FormatWriter() = default;
    
    [[nodiscard]] virtual auto write(
        const std::filesystem::path& file_path,
        const OutputDataset& dataset,
        const OutputConfig& config,
        ProgressCallback progress = nullptr
    ) const -> std::expected<void, OutputError> = 0;
    
    [[nodiscard]] virtual auto get_extension() const noexcept -> std::string_view = 0;
    [[nodiscard]] virtual auto supports_compression() const noexcept -> bool = 0;
    [[nodiscard]] virtual auto supports_metadata() const noexcept -> bool = 0;
};

// Factory for creating format-specific writers
class WriterFactory {
public:
    [[nodiscard]] static auto create_writer(OutputFormat format) 
        -> std::expected<std::unique_ptr<FormatWriter>, UnsupportedFormatError>;
    
    [[nodiscard]] static auto get_available_formats() noexcept 
        -> std::vector<OutputFormat>;
};

// Main output writer class
class OutputWriter {
private:
    OutputConfig config_;
    std::vector<std::unique_ptr<FormatWriter>> writers_;
    
    // Convert solution result to output dataset
    [[nodiscard]] auto convert_solution(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& sim_config,
        const thermophysics::MixtureInterface& mixture
    ) const -> std::expected<OutputDataset, OutputError>;
    
    // Compute derived quantities
    auto compute_wall_properties(
        OutputDataset& dataset,
        const thermophysics::MixtureInterface& mixture
    ) const -> std::expected<void, OutputError>;
    
    auto compute_integrated_quantities(
        OutputDataset& dataset
    ) const -> std::expected<void, OutputError>;
    
    // Generate output file paths
    [[nodiscard]] auto generate_file_paths(
        const std::string& case_name,
        const std::chrono::system_clock::time_point& timestamp
    ) const -> std::vector<std::filesystem::path>;

public:
    explicit OutputWriter(OutputConfig config = {});
    
    // Main interface - write complete solution
    [[nodiscard]] auto write_solution(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& sim_config,
        const thermophysics::MixtureInterface& mixture,
        const std::string& case_name = "simulation",
        ProgressCallback progress = nullptr
    ) -> std::expected<std::vector<std::filesystem::path>, OutputError>;
    
    // Streaming interface for large simulations
    class StreamWriter {
    private:
        const OutputWriter& parent_;
        OutputDataset dataset_;
        std::vector<std::filesystem::path> output_paths_;
        bool finalized_ = false;
        
    public:
        explicit StreamWriter(const OutputWriter& parent, 
                             const Configuration& config,
                             const thermophysics::MixtureInterface& mixture,
                             const std::string& case_name);
        
        [[nodiscard]] auto add_station(
            const blast::boundary_layer::equations::SolutionState& station,
            double xi,
            double x_physical,
            const std::vector<double>& temperature
        ) -> std::expected<void, OutputError>;
        
        [[nodiscard]] auto finalize() -> std::expected<std::vector<std::filesystem::path>, OutputError>;
        
        ~StreamWriter();
    };
    
    [[nodiscard]] auto create_stream_writer(
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::string& case_name = "simulation"
    ) -> std::expected<std::unique_ptr<StreamWriter>, OutputError>;
    
    // Utility functions
    [[nodiscard]] auto get_config() const noexcept -> const OutputConfig& { 
        return config_; 
    }
    
    auto set_config(OutputConfig config) -> void {
        config_ = std::move(config);
        initialize_writers();
    }
    
    // Validation
    [[nodiscard]] auto validate_config() const -> std::expected<void, OutputError>;
    
    // Information
    [[nodiscard]] auto get_output_info(const std::string& case_name = "simulation") const 
        -> std::vector<std::pair<OutputFormat, std::filesystem::path>>;

private:
    auto initialize_writers() -> void;
};

// Convenience functions for common use cases
namespace convenience {

    // Write solution with default HDF5 format
    [[nodiscard]] auto write_hdf5(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& output_path,
        bool include_derivatives = true
    ) -> std::expected<void, OutputError>;
    
    // Write solution in VTK format for visualization
    [[nodiscard]] auto write_vtk(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& output_path
    ) -> std::expected<void, OutputError>;
    
    // Write CSV files for basic analysis
    [[nodiscard]] auto write_csv_series(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& output_directory
    ) -> std::expected<std::vector<std::filesystem::path>, OutputError>;
    
    // Quick write with auto-detection of format from extension
    [[nodiscard]] auto write_auto(
        const blast::boundary_layer::solver::SolutionResult& solution,
        const Configuration& config,
        const thermophysics::MixtureInterface& mixture,
        const std::filesystem::path& output_path
    ) -> std::expected<void, OutputError>;

} // namespace convenience

} // namespace blast::io::output