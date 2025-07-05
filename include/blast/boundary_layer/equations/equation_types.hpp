#pragma once
#include "../../core/containers.hpp"
#include "../../core/exceptions.hpp"
#include <vector>
#include <span>
#include <concepts>

namespace blast::boundary_layer::equations {

// Physical quantity concept for type safety
template<typename T>
concept PhysicalQuantity = std::floating_point<T> && requires(T t) {
    { t >= T{0} } -> std::convertible_to<bool>;
};

// Numeric range concept
template<typename Range>
concept NumericRange = std::ranges::range<Range> && 
                       std::floating_point<std::ranges::range_value_t<Range>>;

// Complete solution state for boundary layer equations
struct SolutionState {
    std::vector<double> V;          // Continuity (velocity-like variable)
    std::vector<double> F;          // Momentum (dimensionless stream function)
    std::vector<double> g;          // Energy (dimensionless enthalpy)
    core::Matrix<double> c;         // Species mass fractions [n_species x n_eta]
    
    // Constructor for initialization
    SolutionState(std::size_t n_eta, std::size_t n_species) 
        : V(n_eta), F(n_eta), g(n_eta), c(n_species, n_eta) {}
    
    // Default constructor
    SolutionState() = default;
};

// Geometry factors for different body types and stations
struct GeometryFactors {
    double J_fact;      // Diffusion flux factor
    double W_fact;      // Chemical production factor  
    double bc_fact;     // Boundary condition factor
    
    constexpr GeometryFactors(double j, double w, double bc) noexcept
        : J_fact(j), W_fact(w), bc_fact(bc) {}
};

// Error type for equation solving operations
class EquationError : public core::BlastException {
public:
    explicit EquationError(std::string_view message,
                          std::source_location location = std::source_location::current())
        : BlastException(std::format("Equation Error: {}", message), location) {}
        
    template<typename... Args>
    explicit EquationError(std::string_view format_str, 
                          std::source_location location,
                          Args&&... args) 
        : BlastException(std::format("Equation Error: {}", 
                        std::vformat(format_str, std::make_format_args(args...))), location) {}
};

} // namespace blast::boundary_layer::equations