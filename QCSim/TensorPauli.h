#pragma once

#include <Eigen/Core>
#include <complex>
#include <stdexcept>
#include <string>

namespace QC
{
namespace TensorNetworks
{
// Optional observable environment payload is bounded per call. No state-dependent environment survives a public batch
// query.
constexpr size_t ObservableCacheBytes = 4 * 1024 * 1024;

inline char CanonicalPauli(char p)
{
    switch (p)
    {
    case 'I':
    case 'i':
        return 'I';
    case 'X':
    case 'x':
        return 'X';
    case 'Y':
    case 'y':
        return 'Y';
    case 'Z':
    case 'z':
        return 'Z';
    default:
        throw std::invalid_argument("Invalid Pauli character");
    }
}

inline void ValidatePauliString(const std::string &pauli, size_t sites)
{
    if (pauli.size() != sites)
        throw std::invalid_argument("Pauli string length must match the number of qubits");
    for (char p : pauli)
        CanonicalPauli(p);
}

inline Eigen::Matrix2cd PauliOperator(char p)
{
    Eigen::Matrix2cd result = Eigen::Matrix2cd::Zero();
    switch (CanonicalPauli(p))
    {
    case 'I':
        result.setIdentity();
        break;
    case 'X':
        result(0, 1) = result(1, 0) = 1.;
        break;
    case 'Y':
        result(0, 1) = std::complex<double>(0., -1.);
        result(1, 0) = std::complex<double>(0., 1.);
        break;
    case 'Z':
        result(0, 0) = 1.;
        result(1, 1) = -1.;
        break;
    }
    return result;
}
} // namespace TensorNetworks
} // namespace QC
