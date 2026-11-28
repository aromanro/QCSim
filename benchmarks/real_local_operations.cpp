#include "DensityMatrix.h"

#include <chrono>
#include <iomanip>
#include <random>

using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;
using Gate = QC::Gates::AppliedGate<>;
using Density = QC::DensityMatrix<>;
using State = QC::QubitRegister<>;
using Targets = std::array<size_t, 3>;
static std::vector<double> batchTimes;

static Matrix unitary(int dimension, bool complex, unsigned seed)
{
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal;
    Eigen::MatrixXd raw(dimension, dimension);
    for (int c = 0; c < dimension; ++c)
        for (int r = 0; r < dimension; ++r)
            raw(r, c) = normal(generator);
    Eigen::MatrixXd orthogonal = raw.householderQr().householderQ();
    Matrix result = orthogonal.cast<std::complex<double>>();
    if (complex)
        for (int r = 0; r < dimension; ++r)
            result.row(r) *= std::polar(1., 0.3 + 0.7 * r);
    return result;
}

static std::vector<Matrix> channel(int dimension, bool complex)
{
    Eigen::VectorXd p(dimension);
    for (int i = 0; i < dimension; ++i)
        p[i] = .18 + .64 * i / (dimension - 1);
    const auto v = unitary(dimension, complex, 33);
    return {unitary(dimension, complex, 31) * p.cwiseSqrt().asDiagonal() * v,
            unitary(dimension, complex, 32) * (Eigen::VectorXd::Ones(dimension) - p).cwiseSqrt().asDiagonal() * v};
}

static Vector initial(int qubits)
{
    Vector psi(1 << qubits);
    std::mt19937 generator(801);
    std::normal_distribution<double> normal;
    for (auto &amplitude : psi)
    {
        const auto real = normal(generator);
        amplitude = {real, normal(generator)};
    }
    return psi.normalized();
}

template <class Function> static double timeOperation(Function &&operation, int batches, int repetitions)
{
    operation(0);
    batchTimes.clear();
    for (int b = 0; b < batches; ++b)
    {
        const auto start = std::chrono::steady_clock::now();
        for (int k = 0; k < repetitions; ++k)
            operation(1 + k + b * repetitions);
        const auto elapsed = std::chrono::steady_clock::now() - start;
        batchTimes.push_back(std::chrono::duration<double, std::micro>(elapsed).count() / repetitions);
    }
    auto times = batchTimes;
    std::sort(times.begin(), times.end());
    return times[times.size() / 2];
}

#ifndef NATIVE_THREE_KRAUS
static void applyTrajectoryReference(State &state, const std::vector<Matrix> &ops, double draw, const Targets &q)
{
    const Matrix reduced = state.LocalDensityMatrix(q[0], q[1], q[2]);
    std::array<double, 2> weights;
    for (int k = 0; k < 2; ++k)
        weights[k] = (ops[k] * reduced * ops[k].adjoint()).trace().real();
    const auto k = draw * (weights[0] + weights[1]) < weights[0] ? 0 : 1;
    state.ApplyGate(Gate(Matrix(ops[k] / std::sqrt(weights[k])), q[0], q[1], q[2]));
}
#endif

int main(int argc, char **argv)
{
    const std::string filter = argc > 1 ? argv[1] : "all";
    const int batches = argc > 2 ? std::stoi(argv[2]) : 5;
    const int repetitions = argc > 3 ? std::stoi(argv[3]) : 8;
    const int densityQubits = argc > 4 ? std::stoi(argv[4]) : 9;
    const int stateQubits = argc > 5 ? std::stoi(argv[5]) : 18;
    const int threads = argc > 6 ? std::stoi(argv[6]) : 1;
    if (threads != 1 && threads != 4)
        throw std::invalid_argument("Select one or four threads");
#if defined(_WIN32) && defined(PIN_BENCH)
    const DWORD_PTR affinity = threads == 1 ? 0x10 : 0x550; // P-core threads 4, 6, 8 and 10 on the measurement host
    if (!SetProcessAffinityMask(GetCurrentProcess(), affinity))
        throw std::runtime_error("Cannot set benchmark processor affinity");
#endif
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    State::SetParallelThreads(threads);
    Eigen::setNbThreads(1);
    std::cout << "operation,width,coefficient,placement,qubits,threads,us,checksum,batch_us\n";
    for (const std::string op :
         {"density-gate", "density-channel", "state-channel", "state-local-norm", "density-sparse", "state-sparse"})
    {
        if (filter != "all" && filter != "real" && op != filter)
            continue;
        const bool density = op.compare(0, 7, "density") == 0;
        const bool sparse = op == "density-sparse" || op == "state-sparse";
        const int qubits = density ? densityQubits : stateQubits;
        const Vector input = initial(qubits);
        for (int width = op == "density-channel" ? 1 : 2; width <= (sparse ? 2 : 3); ++width)
            for (bool complex : {false, true})
                for (bool separated : {false, true})
                {
                    if (sparse && complex)
                        continue;
                    if (filter == "real" && (complex || sparse || width < 2 || op == "state-local-norm"))
                        continue;
                    const Targets q =
                        separated ? Targets{size_t(qubits - 2), 2, size_t(qubits / 2 + 1)} : Targets{0, 1, 2};
                    const Matrix u = unitary(1 << width, complex, 30);
                    auto ops = channel(1 << width, complex);
                    if (sparse)
                    {
                        const int dimension = 1 << width;
                        ops[0] = Matrix::Identity(dimension, dimension) * std::sqrt(.6);
                        ops[1] = Matrix::Zero(dimension, dimension);
                        for (int j = 0; j < dimension; ++j)
                            ops[1](j ^ (dimension - 1), j) = std::sqrt(.4);
                    }
                    double us, checksum;
                    if (density)
                    {
                        Density simulator(qubits, 1);
                        simulator.SetMultithreading(threads > 1);
                        simulator.setFromStatevector(input);
                        if (op == "density-gate")
                        {
                            const Gate gate(u, q[0], q[1], q[2]);
                            us = timeOperation([&](int) { simulator.ApplyGate(gate); }, batches, repetitions);
                        }
                        else if (width < 3)
                            us = timeOperation([&](int) { simulator.ApplyChannel(ops, q[0], q[1]); }, batches,
                                               repetitions);
                        else
                        {
#ifdef NATIVE_THREE_KRAUS
                            us = timeOperation([&](int) { simulator.ApplyChannel(ops, q[0], q[1], q[2]); }, batches,
                                               repetitions);
#else
                            Density branch(qubits, 1);
                            branch.SetMultithreading(threads > 1);
                            Matrix sum(input.size(), input.size());
                            const Gate k0(ops[0], q[0], q[1], q[2]), k1(ops[1], q[0], q[1], q[2]);
                            us = timeOperation(
                                [&](int) {
                                    branch.setDensityMatrix(simulator.getDensityMatrix());
                                    simulator.ApplyGate(k0);
                                    branch.ApplyGate(k1);
                                    sum = simulator.getDensityMatrix() + branch.getDensityMatrix();
                                    simulator.setDensityMatrix(sum);
                                },
                                batches, repetitions);
#endif
                        }
                        const auto &rho = simulator.getDensityMatrix();
                        checksum = rho.trace().real() + rho(1, 6).real() + rho(10, 3).imag();
                    }
                    else
                    {
                        State simulator(qubits, 1);
                        simulator.SetMultithreading(threads > 1);
                        simulator.setRegisterStorage(input);
                        double norm = 0.;
                        if (op == "state-channel" || sparse)
                            us = timeOperation(
                                [&](int index) {
                                    const double draw = index % 2 ? .87 : .13;
                                    if (width == 2)
                                        simulator.ApplyKrausOperator(ops, draw, q[0], q[1]);
                                    else
#ifdef NATIVE_THREE_KRAUS
                                        simulator.ApplyKrausOperator(ops, draw, q[0], q[1], q[2]);
#else
                                        applyTrajectoryReference(simulator, ops, draw, q);
#endif
                                },
                                batches, repetitions);
                        else
                            us = timeOperation(
                                [&](int) {
                                    norm = width == 2 ? simulator.ApplyLocalMapNorm(u, q[0], q[1])
                                                      : simulator.ApplyLocalMapNorm(u, q[0], q[1], q[2]);
                                },
                                batches, repetitions);
                        checksum = norm + simulator.getRegisterStorage().squaredNorm() +
                                   simulator.getBasisStateAmplitude(1).real() +
                                   simulator.getBasisStateAmplitude(6).imag();
                    }
                    std::cout << op << ',' << width << ',' << (complex ? "complex" : "real") << ','
                              << (separated ? "spread" : "low") << ',' << qubits << ',' << threads << ','
                              << std::setprecision(10) << us << ',' << std::setprecision(16) << checksum << ',';
                    for (size_t i = 0; i < batchTimes.size(); ++i)
                        std::cout << (i ? "|" : "") << std::setprecision(10) << batchTimes[i];
                    std::cout << std::endl;
                }
    }
}
