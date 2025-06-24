
#pragma once
#include "../../../core/containers.hpp"
#include "../../../core/exceptions.hpp"
#include "../../../io/config_types.hpp"
#include "coordinate_transform.hpp"
#include <vector>
#include <span>
#include <expected>
#include <concepts>
#include <ranges>

namespace blast::boundary_layer::grid {

template<typename T>
concept NumericRange = std::ranges::range<T> && // concept is like en exigeance. template<"Concept" T> then a function 
                      std::floating_point<std::ranges::range_value_t<T>>; // contain floats?

template<typename T>
concept GridConfigType = requires(const T& config) {
    config.n_eta;
    config.eta_max;
};

// Compile-time constants
namespace constants {
    constexpr double default_tolerance = 1e-6;
    constexpr size_t default_grid_reserve = 100;
    constexpr double step_reduction_factor = 0.5;
}


class BoundaryLayerGrid {
private:
    std::vector<double> xi_;
    std::vector<double> xi_output_;
    std::vector<double> eta_;
    double d_eta_;
    int n_eta_;
    double eta_max_;

public:
    template<GridConfigType NumericalConfig> // NumericalConfig is restrained by the type GridConfigType
    [[nodiscard]] static constexpr auto create_stagnation_grid(
        NumericalConfig&& numerical_config, // different from io::NumericalConfig
        const io::OuterEdgeConfig& edge_config
    ) -> std::expected<BoundaryLayerGrid, GridError>; // some kind of constructior, it's a factory
    
    template<GridConfigType NumericalConfig>
    [[nodiscard]] static constexpr auto create_downstream_grid(
        NumericalConfig&& numerical_config,
        const io::OuterEdgeConfig& edge_config,
        const io::OutputConfig& output_config
    ) -> std::expected<BoundaryLayerGrid, GridError>;
    
    // Modern accessors with noexcept
    [[nodiscard]] constexpr auto xi_coordinates() const noexcept -> std::span<const double> { return xi_; } // std::span is a light view, so no copy
    [[nodiscard]] constexpr auto xi_output_coordinates() const noexcept -> std::span<const double> { return xi_output_; }
    [[nodiscard]] constexpr auto eta_coordinates() const noexcept -> std::span<const double> { return eta_; }
    
    [[nodiscard]] constexpr auto n_eta() const noexcept -> int { return n_eta_; }
    [[nodiscard]] constexpr auto eta_max() const noexcept -> double { return eta_max_; }
    [[nodiscard]] constexpr auto d_eta() const noexcept -> double { return d_eta_; }
    
    // Grid operations
    [[nodiscard]] auto find_xi_interval(double xi_target) const noexcept
        -> std::expected<std::pair<int, int>, GridError>;
    
    template<NumericRange XGrid>
    [[nodiscard]] auto interpolate_x_from_xi(double xi_target, XGrid&& x_grid) const 
        -> std::expected<double, GridError>;

private:
    explicit constexpr BoundaryLayerGrid(int n_eta, double eta_max) noexcept;
    
    constexpr auto generate_eta_distribution() noexcept -> void;
    auto generate_xi_distribution(const io::OuterEdgeConfig& edge_config) 
        -> std::expected<void, GridError>;
    auto generate_xi_output_distribution(const io::OutputConfig& output_config,
                                        std::span<const double> x_edge) 
        -> std::expected<void, GridError>;
};

// Utility functions with constexpr where possible
[[nodiscard]] constexpr auto compute_xi_step_size(
    double current_xi, double d_xi, int iterations,
    const io::NumericalConfig::StepControl& step_control
) noexcept -> double;

[[nodiscard]] constexpr auto reduce_step_size(double current_xi, double& d_xi) noexcept -> double;

} // namespace blast::boundary_layer::grid