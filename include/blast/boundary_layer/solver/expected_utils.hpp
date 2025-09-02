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
 * @brief Assign-or-return helper (portable)
 *
 * Usage:
 *   BLAST_TRY_ASSIGN(value, some_expected_result);
 * Expands to:
 *   auto __tmp = some_expected_result;
 *   if (!__tmp) return std::unexpected(__tmp.error());
 *   value = std::move(__tmp.value());
 */
#define BLAST_TRY_ASSIGN(lhs, expr)                                                           \
    do {                                                                                      \
        auto blast_try_tmp = (expr);                                                          \
        if (!blast_try_tmp)                                                                   \
            return std::unexpected(blast_try_tmp.error());                                    \
        lhs = std::move(blast_try_tmp.value());                                               \
    } while (0)

/**
 * @brief Deprecated alias for non-portable BLAST_TRY.
 *
 * The old GNU statement-expression version allowed: auto v = BLAST_TRY(expr);
 * That construct is not portable (fails on MSVC). Use BLAST_TRY_ASSIGN instead:
 *   Type v; BLAST_TRY_ASSIGN(v, expr);
 */
#define BLAST_TRY(...) static_assert(false, "Use BLAST_TRY_ASSIGN(var, expr) instead of BLAST_TRY(expr)")

/**
 * @brief BLAST_TRY_VOID macro for void expected results
 * 
 * Usage: BLAST_TRY_VOID(some_void_expected_result);
 * Replaces:
 *   auto result = some_void_expected_result;
 *   if (!result) return std::unexpected(result.error());
 */
#define BLAST_TRY_VOID(expr)                                                                  \
    do {                                                                                      \
        auto blast_try_tmp_void = (expr);                                                     \
        if (!blast_try_tmp_void)                                                              \
            return std::unexpected(blast_try_tmp_void.error());                               \
    } while (0)

/**
 * @brief BLAST_TRY_WITH_CONTEXT macro for error propagation with additional context
 * 
 * Usage: auto value = BLAST_TRY_WITH_CONTEXT(some_expected_result, "Failed to compute X");
 */
/**
 * @brief Assign-or-return helper with added error context (portable)
 *
 * Usage:
 *   BLAST_TRY_ASSIGN_CTX(value, expr, "message");
 */
#define BLAST_TRY_ASSIGN_CTX(lhs, expr, context)                                              \
    do {                                                                                      \
        auto blast_try_tmp_ctx = (expr);                                                      \
        if (!blast_try_tmp_ctx) {                                                             \
            return std::unexpected(NumericError(                                              \
                std::format("{}: {}", (context), blast_try_tmp_ctx.error().message())));     \
        }                                                                                     \
        lhs = std::move(blast_try_tmp_ctx.value());                                           \
    } while (0)

/**
 * @brief Void-or-return helper with added error context (portable)
 */
#define BLAST_TRY_VOID_CTX(expr, context)                                                     \
    do {                                                                                      \
        auto blast_try_tmp_vctx = (expr);                                                     \
        if (!blast_try_tmp_vctx) {                                                            \
            return std::unexpected(NumericError(                                              \
                std::format("{}: {}", (context), blast_try_tmp_vctx.error().message())));    \
        }                                                                                     \
    } while (0)

// Deprecated alias (non-portable before). Force a clear error if used.
#define BLAST_TRY_WITH_CONTEXT(...) static_assert(false, "Use BLAST_TRY_ASSIGN_CTX(lhs, expr, ctx) or BLAST_TRY_VOID_CTX(expr, ctx)")

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
