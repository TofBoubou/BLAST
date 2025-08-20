#pragma once

#include "blast/boundary_layer/solver/solver_errors.hpp"
#include <expected>
#include <format>
#include <functional>
#include <type_traits>

namespace blast::boundary_layer::solver {

/**
 * @brief Utilities for working with std::expected to reduce boilerplate
 * 
 * This provides monadic operations and helpers to eliminate the repetitive
 * pattern of checking std::expected results found throughout the codebase.
 */
namespace expected_utils {

/**
 * @brief BLAST_TRY macro for simplified error propagation
 * 
 * Usage: auto value = BLAST_TRY(some_expected_result);
 * Replaces:
 *   auto result = some_expected_result;
 *   if (!result) return std::unexpected(result.error());
 *   auto value = result.value();
 */
#define BLAST_TRY(expr) \
    ({ \
        auto __tmp = (expr); \
        if (!__tmp) return std::unexpected(__tmp.error()); \
        if constexpr (!std::is_void_v<decltype(__tmp.value())>) { \
            std::move(__tmp.value()); \
        } \
    })

/**
 * @brief BLAST_TRY_VOID macro for void expected results
 * 
 * Usage: BLAST_TRY_VOID(some_void_expected_result);
 * Replaces:
 *   auto result = some_void_expected_result;
 *   if (!result) return std::unexpected(result.error());
 */
#define BLAST_TRY_VOID(expr) \
    do { \
        auto __tmp = (expr); \
        if (!__tmp) return std::unexpected(__tmp.error()); \
    } while(0)

/**
 * @brief BLAST_TRY_WITH_CONTEXT macro for error propagation with additional context
 * 
 * Usage: auto value = BLAST_TRY_WITH_CONTEXT(some_expected_result, "Failed to compute X");
 */
#define BLAST_TRY_WITH_CONTEXT(expr, context) \
    ({ \
        auto __tmp = (expr); \
        if (!__tmp) { \
            return std::unexpected(NumericError( \
                std::format("{}: {}", (context), __tmp.error().message()))); \
        } \
        std::move(__tmp.value()); \
    })

/**
 * @brief Chain multiple expected operations with automatic error propagation
 * 
 * @tparam T Initial value type
 * @tparam F Function type
 * @param initial Initial expected value
 * @param func Function to apply if initial has value
 * @return Result of chaining or first error encountered
 */
template<typename T, typename F>
[[nodiscard]] constexpr auto and_then(std::expected<T, SolverError> initial, F&& func)
    -> decltype(func(std::move(initial.value()))) {
    
    if (!initial) {
        return std::unexpected(initial.error());
    }
    return func(std::move(initial.value()));
}

/**
 * @brief Transform the value of an expected if it has one
 * 
 * @tparam T Initial value type
 * @tparam F Function type
 * @param initial Initial expected value
 * @param func Transformation function
 * @return Transformed expected or original error
 */
template<typename T, typename F>
[[nodiscard]] constexpr auto transform(std::expected<T, SolverError> initial, F&& func)
    -> std::expected<decltype(func(std::move(initial.value()))), SolverError> {
    
    if (!initial) {
        return std::unexpected(initial.error());
    }
    return func(std::move(initial.value()));
}

/**
 * @brief Transform an error to add context
 * 
 * @tparam T Value type
 * @param initial Initial expected value
 * @param context Context to add to error
 * @return Same expected or error with added context
 */
template<typename T>
[[nodiscard]] constexpr auto transform_error(std::expected<T, SolverError> initial, 
                                            const std::string& context)
    -> std::expected<T, SolverError> {
    
    if (!initial) {
        return std::unexpected(NumericError(
            std::format("{}: {}", context, initial.error().message())));
    }
    return initial;
}

/**
 * @brief Collect multiple expected values into a vector
 * 
 * @tparam T Value type
 * @param expecteds Vector of expected values
 * @return Vector of values or first error encountered
 */
template<typename T>
[[nodiscard]] auto collect_all(std::vector<std::expected<T, SolverError>>&& expecteds)
    -> std::expected<std::vector<T>, SolverError> {
    
    std::vector<T> results;
    results.reserve(expecteds.size());
    
    for (auto& expected : expecteds) {
        if (!expected) {
            return std::unexpected(expected.error());
        }
        results.push_back(std::move(expected.value()));
    }
    
    return results;
}

/**
 * @brief Execute a function with error handling and context
 * 
 * @tparam F Function type
 * @param func Function that may throw or return expected
 * @param context Error context
 * @return Expected result with error handling
 */
template<typename F>
[[nodiscard]] auto with_context(F&& func, const std::string& context)
    -> std::expected<decltype(func()), SolverError> {
    
    try {
        if constexpr (requires { func().has_value(); }) {
            // Function returns std::expected
            auto result = func();
            return transform_error(std::move(result), context);
        } else {
            // Function returns value directly
            return func();
        }
    } catch (const std::exception& e) {
        return std::unexpected(NumericError(
            std::format("{}: {}", context, e.what())));
    }
}

} // namespace expected_utils

// Note: Macros BLAST_TRY, BLAST_TRY_VOID and BLAST_TRY_WITH_CONTEXT are available globally

} // namespace blast::boundary_layer::solver