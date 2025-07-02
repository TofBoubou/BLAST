#pragma once
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <span>
#include <memory>
#include <concepts>


namespace blast::core{

template<typename T>
using Vector = std::vector<T>;

template<typename T, std::size_t N>
using Array = std::array<T, N>;

template<typename T>
using Span = std::span<T>;

template<typename T>
using UniquePtr = std::unique_ptr<T>;

template<typename T>
using SharedPtr = std::shared_ptr<T>;

template<typename Scalar = double>
using MathVector = Eigen::Vector<Scalar, Eigen::Dynamic>;

template<typename Scalar = double, int Size = Eigen::Dynamic>
using FixedMathVector = Eigen::Vector<Scalar, Size>;

template<typename Scalar = double>
using MathMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

using PhysicsVector = MathVector<double>;        
using SpeciesVector = MathVector<double>;        
using TemperatureField = MathVector<double>;       

template<typename Scalar = double>
class Matrix{
private:
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> data_;

public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols) : data_(rows, cols) {}

    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;

    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default;

    [[nodiscard]] auto rows() const noexcept -> std::size_t {return data_.rows();}
    [[nodiscard]] auto cols() const noexcept -> std::size_t {return data_.cols();}

    auto operator()(std::size_t i, std::size_t j) -> Scalar& {return data_(i,j);}
    [[nodiscard]] auto operator()(std::size_t i, std::size_t j) const -> const Scalar& {return data_(i,j);}

    [[nodiscard]] auto eigen() -> auto& {return data_;}
    [[nodiscard]] auto eigen() const -> auto& {return data_;}
    void setZero() { data_.setZero(); }
};
}