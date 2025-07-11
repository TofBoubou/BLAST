#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include "../config_types.hpp"
#include "../../boundary_layer/solver/boundary_layer_solver.hpp"
#include <filesystem>
#include <chrono>
#include <variant>
#include <expected>

namespace blast::io::output {

// Output format enumeration
enum class OutputFormat {
    HDF5
};

// Output configuration
struct OutputConfig {
    std::filesystem::path base_directory = "BLAST_outputs";
    std::string case_name = "simulation";
    OutputFormat primary_format = OutputFormat::HDF5;
    std::vector<OutputFormat> additional_formats;
    bool save_metadata = true;
    bool save_derivatives = true;
    bool save_coefficients = false;  // Can be large
    bool compress_data = true;
    double compression_level = 6.0;  // 0-9 for HDF5
    
    // Time stamping
    bool include_timestamp = true;
    
    // Variable selection
    struct VariableSelection {
        bool flow_variables = true;      // F, g, V, T, P, rho
        bool species_concentrations = true;
        bool transport_properties = false;
        bool chemical_rates = false;
        bool diffusion_fluxes = false;
    } variables;
};

// Metadata container
struct SimulationMetadata {
    std::string blast_version = "1.0.0";
    std::chrono::system_clock::time_point creation_time;
    Configuration simulation_config;
    
    struct GridInfo {
        int n_eta;
        double eta_max;
        double d_eta;
        std::vector<double> xi_coordinates;
        std::vector<double> eta_coordinates;
    } grid;
    
    struct MixtureInfo {
        std::vector<std::string> species_names;
        std::vector<double> species_molecular_weights;
        std::vector<double> species_charges;
        bool has_electrons;
    } mixture;
    
    struct ConvergenceInfo {
        bool converged;
        int total_iterations;
        double final_residual;
    } convergence;
};

// Field data for a single station
struct StationData {
    // Coordinates
    double xi;
    double x_physical;
    std::vector<double> eta;
    std::vector<double> y_physical;
    
    // Primary variables
    std::vector<double> F;           // Dimensionless stream function
    std::vector<double> g;           // Dimensionless enthalpy
    std::vector<double> V;           // Velocity field
    std::vector<double> temperature;
    std::vector<double> pressure;
    std::vector<double> density;
    
    // Species data
    core::Matrix<double> species_concentrations;  // [n_species x n_eta]
    
    // Transport properties (optional)
    std::vector<double> viscosity;
    std::vector<double> thermal_conductivity;
    core::Matrix<double> diffusion_coefficients;
    
    // Chemical data (optional)
    core::Matrix<double> production_rates;
    core::Matrix<double> diffusion_fluxes;
    
    // Derivatives (optional)
    struct Derivatives {
        std::vector<double> dF_deta;
        std::vector<double> dg_deta;
        std::vector<double> dT_deta;
        core::Matrix<double> dc_deta;
        
        // Xi derivatives
        std::vector<double> dF_dxi;
        std::vector<double> dg_dxi;
        core::Matrix<double> dc_dxi;
    } derivatives;
};

// Complete output dataset
struct OutputDataset {
    SimulationMetadata metadata;
    std::vector<StationData> stations;
    
    // Wall properties
    struct WallData {
        std::vector<double> x_positions;
        std::vector<double> temperatures;
        std::vector<double> heat_flux;
        std::vector<double> shear_stress;
    } wall;
    
    // Integrated quantities
    struct IntegratedQuantities {
        std::vector<double> displacement_thickness;
        std::vector<double> momentum_thickness;
        std::vector<double> shape_factor;
        std::vector<double> total_enthalpy_thickness;
    } integrated;
};

// Progress callback for large outputs
using ProgressCallback = std::function<void(double progress, const std::string& stage)>;

// Output error types
class OutputError : public core::BlastException {
public:
    explicit OutputError(std::string_view message,
                        std::source_location location = std::source_location::current())
        : BlastException(std::format("Output Error: {}", message), location) {}
};

class FileWriteError : public OutputError {
private:
    std::filesystem::path file_path_;
    
public:
    explicit FileWriteError(const std::filesystem::path& path, std::string_view message,
                           std::source_location location = std::source_location::current())
        : OutputError(std::format("File '{}': {}", path.string(), message), location)
        , file_path_(path) {}
    
    [[nodiscard]] auto file_path() const noexcept -> const std::filesystem::path& { 
        return file_path_; 
    }
};

class UnsupportedFormatError : public OutputError {
public:
    explicit UnsupportedFormatError(OutputFormat format,
                                   std::source_location location = std::source_location::current())
        : OutputError(std::format("Unsupported output format: {}", static_cast<int>(format)), location) {}
};

} // namespace blast::io::output