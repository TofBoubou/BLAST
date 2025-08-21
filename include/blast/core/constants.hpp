#pragma once

#include <cstddef>

namespace blast::constants {

// ================================================================================================
// FUNDAMENTAL PHYSICAL CONSTANTS
// ================================================================================================

namespace physical {
/// Universal gas constant [J/(kmol·K)]
inline constexpr double universal_gas_constant = 8.31446261815324;

/// Boltzmann constant [J/K]
inline constexpr double boltzmann_constant = 1.380649e-23;

/// Avogadro's number [mol^-1]
inline constexpr double avogadro_number = 6.02214076e23;

/// Mathematical constant π
inline constexpr double pi = 3.14159265358979323846;

/// Stefan-Boltzmann constant [W/(m²·K⁴)]
inline constexpr double stefan_boltzmann = 5.670374419e-8;
}  // namespace physical

// ================================================================================================
// NUMERICAL TOLERANCES
// ================================================================================================

namespace tolerance {
/// Standard convergence tolerance for iterative methods
inline constexpr double standard = 1e-6;

/// High precision tolerance for critical computations
inline constexpr double high_precision = 1e-9;

/// Ultra high precision for numerical stability
inline constexpr double ultra_precision = 1e-12;

/// Mass fraction sum tolerance
inline constexpr double mass_fraction_sum = 1e-6;

/// Coordinate transform tolerance
inline constexpr double coordinate_transform = 1e-6;
}  // namespace tolerance

// ================================================================================================
// ITERATION LIMITS
// ================================================================================================

namespace iteration_limits {
/// Default maximum iterations for boundary layer solver
inline constexpr int boundary_layer_max = 1000000;

/// Maximum iterations for enthalpy-temperature solver
inline constexpr int enthalpy_temperature_max = 500000;

/// Maximum iterations for Stefan problems
inline constexpr int stefan_max = 50000;
}  // namespace iteration_limits

// ================================================================================================
// DEFAULT GRID PARAMETERS
// ================================================================================================

namespace grid {
/// Default number of eta points in boundary layer grid
inline constexpr int default_n_eta = 101;

/// Default maximum eta value for boundary layer
inline constexpr double default_eta_max = 8.0;

/// Minimum required points for numerical differentiation
inline constexpr std::size_t min_differentiation_points = 5;

/// Default temperature discretization points for abacus
inline constexpr int default_temperature_points = 100;
}  // namespace grid

// ================================================================================================
// NUMERICAL METHOD COEFFICIENTS
// ================================================================================================

namespace numerical_methods {
/// Simpson's 5-point integration coefficients
namespace simpson_5pt {
inline constexpr double coeff_1 = 17.0;
inline constexpr double coeff_2 = 42.0;
inline constexpr double coeff_3 = -16.0;
inline constexpr double coeff_4 = 6.0;
inline constexpr double coeff_5 = -1.0;
inline constexpr double divisor = 48.0;
}  // namespace simpson_5pt

/// Simpson's 3-point integration coefficients
namespace simpson_3pt {
inline constexpr double coeff_1 = 1.0;
inline constexpr double coeff_2 = 4.0;
inline constexpr double coeff_3 = 1.0;
inline constexpr double divisor = 3.0;
}  // namespace simpson_3pt

/// Trapezoidal rule coefficient
inline constexpr double trapezoidal_coeff = 0.5;

/// Tridiagonal solver and collocation method constants
namespace tridiagonal_solver {
/// Collocation point positions
inline constexpr double left_point = -1.0;
inline constexpr double center_point = 0.0;
inline constexpr double right_point = 1.0;

/// Common collocation coefficients
inline constexpr double half_coeff = 0.5;
inline constexpr double one_half_coeff = 1.5;
inline constexpr double minus_half_coeff = -0.5;
inline constexpr double minus_one_half_coeff = -1.5;
inline constexpr double two_coeff = 2.0;
inline constexpr double minus_two_coeff = -2.0;
inline constexpr double six_coeff = 6.0;
inline constexpr double minus_six_coeff = -6.0;
inline constexpr double ten_coeff = 10.0;
inline constexpr double minus_ten_coeff = -10.0;

/// Numerical tolerance for matrix operations
inline constexpr double diagonal_tolerance = 1e-14;
inline constexpr double solution_normalization_tolerance = 1e-9;
}  // namespace tridiagonal_solver
}  // namespace numerical_methods

// ================================================================================================
// DEFAULT PHYSICAL CONDITIONS
// ================================================================================================

namespace defaults {
/// Default ambient temperature [K]
inline constexpr double ambient_temperature = 300.0;

/// Default wall temperature for stable conditions [K]
inline constexpr double stable_wall_temperature = 3100.0;

/// Default edge temperature for stable conditions [K]
inline constexpr double stable_edge_temperature = 3100.0;

/// Default pressure for stable conditions [Pa]
inline constexpr double stable_pressure = 7000.0;

/// Default emissivity for non-radiating surfaces
inline constexpr double default_emissivity = 0.0;

/// Maximum valid emissivity value
inline constexpr double max_emissivity = 1.0;

/// Default abacus temperature range
namespace abacus {
inline constexpr double temperature_min = 300.0;
inline constexpr double temperature_max = 2000.0;
}  // namespace abacus
}  // namespace defaults

// ================================================================================================
// FILE I/O AND FORMATTING CONSTANTS
// ================================================================================================

namespace io {
/// HDF5 default compression level (0-9, higher = better compression)
inline constexpr int default_hdf5_compression = 6;

/// Default HDF5 chunk size for datasets
inline constexpr std::size_t default_hdf5_chunk_size = 1024;

/// Bytes to KB conversion factor
inline constexpr double bytes_to_kb = 1024.0;

/// Invalid HDF5 handle value
inline constexpr int invalid_hdf5_handle = -1;

/// Number of dimensions for temperature range attributes
inline constexpr std::size_t temperature_range_dims = 2;

/// Default BLAST version string
inline constexpr const char* default_blast_version = "1.0.0";
}  // namespace io

// ================================================================================================
// ARRAY AND INDEXING CONSTANTS
// ================================================================================================

namespace indexing {
/// First array index
inline constexpr std::size_t first = 0;

/// Second array index  
inline constexpr std::size_t second = 1;

/// Third array index
inline constexpr std::size_t third = 2;

/// Step size for progress reporting
inline constexpr int progress_report_step = 10;

/// Minimum required array size for range validation
inline constexpr std::size_t min_range_size = 2;

/// Number of variables for single temperature mode
inline constexpr std::size_t single_temp_vars = 2;

/// Maximum number of profile points to print
inline constexpr std::size_t max_profile_print_points = 11;
}  // namespace indexing

// ================================================================================================
// HERMITE INTERPOLATION COEFFICIENTS
// ================================================================================================

namespace hermite {
/// Hermite basis function coefficients for cubic interpolation
namespace basis_functions {
inline constexpr double h00_coeff_t3 = 2.0;
inline constexpr double h00_coeff_t2 = -3.0;
inline constexpr double h00_constant = 1.0;

inline constexpr double h10_coeff_t3 = 1.0;
inline constexpr double h10_coeff_t2 = -2.0;

inline constexpr double h01_coeff_t3 = -2.0;
inline constexpr double h01_coeff_t2 = 3.0;

inline constexpr double h11_coeff_t3 = 1.0;
inline constexpr double h11_coeff_t2 = -1.0;
}  // namespace basis_functions
}  // namespace hermite

// ================================================================================================
// STRING PROCESSING CONSTANTS
// ================================================================================================

namespace string_processing {
/// Length of comma-space separator for options
inline constexpr std::size_t option_separator_length = 2;

/// Format precision for floating point display
inline constexpr int float_precision_2 = 2;
inline constexpr int float_precision_3 = 3;
inline constexpr int float_precision_4 = 4;

/// Field widths for tabular output
inline constexpr int narrow_field_width = 6;
inline constexpr int medium_field_width = 8;
inline constexpr int wide_field_width = 10;
inline constexpr int separator_width = 36;
}  // namespace string_processing

// ================================================================================================
// UNIT CONVERSION FACTORS
// ================================================================================================

namespace conversion {
/// Progress percentage conversion factor
inline constexpr double to_percentage = 100.0;
}  // namespace conversion

}  // namespace blast::constants