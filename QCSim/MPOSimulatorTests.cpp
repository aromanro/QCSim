#include <iostream>
#include <iterator>
#include <memory>
#include <vector>
#include <unordered_map>
#include <utility>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <set>

#include "Tests.h"

#include "QubitRegister.h"
#include "MPOSimulator.h"
#include "DensityMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>


#define NR_QUBITS_LIMIT_MPO 7

template<class Callable>
static bool MPO_ExpectInvalidArgument(Callable&& callable, const char* description)
{
	try
	{
		callable();
	}
	catch (const std::invalid_argument&)
	{
		return true;
	}
	catch (const std::exception& ex)
	{
		std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
		return false;
	}

	std::cout << description << " did not throw std::invalid_argument" << std::endl;
	return false;
}

template<class Callable>
static bool MPO_ExpectRuntimeError(Callable&& callable, const char* description)
{
	try
	{
		callable();
	}
	catch (const std::runtime_error&)
	{
		return true;
	}
	catch (const std::exception& ex)
	{
		std::cout << description << " threw the wrong exception: " << ex.what() << std::endl;
		return false;
	}

	std::cout << description << " did not throw std::runtime_error" << std::endl;
	return false;
}

class MPOWrongDeclaredArityGate final : public QC::Gates::AppliedGate<>
{
public:
	explicit MPOWrongDeclaredArityGate(const Eigen::MatrixXcd& op)
		: QC::Gates::AppliedGate<>(op)
	{
	}

	size_t getQubitsNumber() const override
	{
		return 1;
	}
};

// these mirror the helpers from MPSSimulatorTests.cpp, kept local to avoid coupling the two test files
static void FillOneQubitGatesMPO(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::HadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::HyGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SDGGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::TGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::TDGGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PhaseShiftGate<>>(0.38 * M_PI));
	gates.emplace_back(std::make_shared<QC::Gates::PauliXGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PauliYGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::PauliZGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::SquareRootNOTDagGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::RxGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::RyGate<>>(M_PI / 7));
	gates.emplace_back(std::make_shared<QC::Gates::RzGate<>>(M_PI / 5));
	gates.emplace_back(std::make_shared<QC::Gates::UGate<>>(M_PI / 3, M_PI / 5));
}

static void FillTwoQubitGatesMPO(std::vector<std::shared_ptr<QC::Gates::QuantumGateWithOp<>>>& gates)
{
	gates.emplace_back(std::make_shared<QC::Gates::SwapGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::iSwapGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::CNOTGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledYGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledZGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledHadamardGate<>>());
	gates.emplace_back(std::make_shared<QC::Gates::ControlledPhaseShiftGate<>>(M_PI / 3));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledUGate<>>(M_PI / 3, M_PI / 7));
	gates.emplace_back(std::make_shared<QC::Gates::ControlledRxGate<>>(M_PI / 5));
}

// builds the reference density matrix rho = |psi><psi| from the statevector simulator
static Eigen::MatrixXcd ReferenceDensityMatrix(const QC::QubitRegister<>& reg)
{
	return reg.getDensityMatrix();
}

static bool CompareDensityMatrices(const Eigen::MatrixXcd& rhoRef, const Eigen::MatrixXcd& rhoMPO, int nrQubits, double err = 1E-3)
{
	if (rhoRef.rows() != rhoMPO.rows() || rhoRef.cols() != rhoMPO.cols())
	{
		std::cout << "Density matrix dimensions mismatch for the MPO simulator for " << nrQubits << " qubits" << std::endl;
		return false;
	}

	for (Eigen::Index r = 0; r < rhoRef.rows(); ++r)
		for (Eigen::Index c = 0; c < rhoRef.cols(); ++c)
			if (!approxEqual(rhoRef(r, c), rhoMPO(r, c), err))
			{
				std::cout << "Density matrix element (" << r << ", " << c << ") simulation test failed for the MPO simulator for " << nrQubits << " qubits" << std::endl;
				std::cout << "Reference: " << rhoRef(r, c) << " vs MPO: " << rhoMPO(r, c) << std::endl;

				return false;
			}

	return true;
}

// ... unchanged file omitted for brevity ...

static bool TraceRestoreAndHermitizeTestMPO()
{
	std::cout << "\nMPO simulator optional trace restore and hermitization" << std::endl;

	{
		QC::TensorNetworks::MPOSimulator defaults(2);
		using KrausCompletenessCheck = QC::TensorNetworks::MPOSimulatorInterface::KrausCompletenessCheck;
		if (defaults.getKrausCompletenessCheck() != KrausCompletenessCheck::Ignore ||
			defaults.getRestoreTraceAfterTruncation() || defaults.getHermitizeAfterTruncation())
		{
			std::cout << "Optional MPO patches are not off by default" << std::endl;
			return false;
		}
	}

	QC::Gates::CNOTGate<> cnot;

	QC::TensorNetworks::MPOSimulator drifted(2);
	drifted.setToMixtureOfBasisStates(
		std::vector<std::pair<QC::TensorNetworks::MPOSimulatorInterface::IndexType, double>>{
			{0, 0.1}, {1, 0.2}, {2, 0.3}, {3, 0.4}
		});
	drifted.setLimitBondDimension(1);
	drifted.Trim();
	if (!std::isfinite(drifted.Trace().real()) || !std::isfinite(drifted.Trace().imag()))
	{
		std::cout << "Truncation produced a non-finite trace" << std::endl;
		return false;
	}
	if (approxEqual(drifted.Trace(), std::complex<double>(1., 0.), 1E-10))
	{
		std::cout << "Truncation of a four-term mixture at chi=1 did not drift the trace" << std::endl;
		return false;
	}

	QC::TensorNetworks::MPOSimulator restored(2);
	restored.setToMixtureOfBasisStates(
		std::vector<std::pair<QC::TensorNetworks::MPOSimulatorInterface::IndexType, double>>{
			{0, 0.1}, {1, 0.2}, {2, 0.3}, {3, 0.4}
		});
	restored.setRestoreTraceAfterTruncation(true);
	restored.setLimitBondDimension(1);
	restored.Trim();
	if (!approxEqual(restored.Trace(), std::complex<double>(1., 0.), 1E-10))
	{
		std::cout << "Trace restore after Trim did not return Tr(rho) to 1: " << restored.Trace() << std::endl;
		return false;
	}

	QC::TensorNetworks::MPOSimulator restoredGate(2);
	restoredGate.setToMixtureOfBasisStates(
		std::vector<std::pair<QC::TensorNetworks::MPOSimulatorInterface::IndexType, double>>{
			{0, 0.1}, {1, 0.2}, {2, 0.3}, {3, 0.4}
		});
	restoredGate.setRestoreTraceAfterTruncation(true);
	restoredGate.setLimitBondDimension(1);
	restoredGate.ApplyGate(cnot, 1, 0);
	if (!approxEqual(restoredGate.Trace(), std::complex<double>(1., 0.), 1E-10))
	{
		std::cout << "Trace restore after a truncating two-qubit gate did not return Tr(rho) to 1: "
			<< restoredGate.Trace() << std::endl;
		return false;
	}

	// ... rest unchanged ...
}

static bool DiagnosticsTestMPO()
{
	std::cout << "\nMPO simulator diagnostics (trace, Tr(rho^2), Hermiticity)" << std::endl;

	QC::TensorNetworks::MPOSimulator pure(1);
	if (!approxEqual(pure.Trace(), std::complex<double>(1., 0.), 1E-12) ||
		!approxEqual(pure.TraceOfSquare(), std::complex<double>(1., 0.), 1E-12) ||
		!approxEqual(pure.Purity(), 1., 1E-12) ||
		!pure.IsHermitian(1E-12))
	{
		std::cout << "Computational-basis projector failed the MPO diagnostics" << std::endl;
		return false;
	}

	QC::TensorNetworks::MPOSimulator mixed(1);
	mixed.setToMixtureOfBasisStates(
		std::vector<std::pair<QC::TensorNetworks::MPOSimulatorInterface::IndexType, double>>{
			{0, 0.5}, {1, 0.5}
		});
	if (!approxEqual(mixed.TraceOfSquare(), std::complex<double>(0.5, 0.), 1E-12) ||
		!approxEqual(mixed.Purity(), 0.5, 1E-12))
	{
		std::cout << "Maximally mixed qubit failed Tr(rho^2) / purity diagnostics" << std::endl;
		return false;
	}

	// ... rest unchanged ...
}
