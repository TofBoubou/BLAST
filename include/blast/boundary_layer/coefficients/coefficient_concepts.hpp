#pragma once
#include <concepts>
#include <span>
#include <expected>

namespace blast::boundary_layer::coefficients {

// Physical quantity constraints
template<typename T>
concept Temperature = std::floating_point<T> && requires(T t) {
    { t > T{0} } -> std::convertible_to<bool>;
};

template<typename T>
concept Pressure = std::floating_point<T> && requires(T t) {
    { t > T{0} } -> std::convertible_to<bool>;
};

template<typename T>
concept Density = std::floating_point<T> && requires(T t) {
    { t > T{0} } -> std::convertible_to<bool>;
};

template<typename T>
concept MassFraction = std::floating_point<T> && requires(T t) {
    { t >= T{0} && t <= T{1} } -> std::convertible_to<bool>;
};

// Container concepts
template<typename T>
concept SpeciesContainer = requires(T t) {
    typename T::value_type;
    requires MassFraction<typename T::value_type>;
    { t.size() } -> std::convertible_to<std::size_t>;
};

// Coefficient error type
class CoefficientError : public std::runtime_error {
public:
    explicit CoefficientError(std::string_view message)
        : std::runtime_error(std::string(message)) {}
};

} // namespace blast::boundary_layer::coefficients