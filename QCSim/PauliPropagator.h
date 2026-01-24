#pragma once

#include "PauliStringXZCoeff.h"

#include <unordered_set>
#include <string>
#include <random>

namespace QC
{
	struct PauliStringXZHash
	{
		std::size_t operator()(const PauliStringXZ& p) const noexcept
		{
			std::size_t h1 = std::hash<std::vector<bool>>{}(p.X);
			std::size_t h2 = std::hash<std::vector<bool>>{}(p.Z);
			return h1 ^ (h2 << 1);
		}
	};

	enum class OperationType : unsigned char
	{
		X,
		Y,
		Z,
		H,
		K,
		S,
		SDG,
		SX,
		SXDG,
		CX,
		CY,
		CZ,
		SWAP,
		ISWAP,
		ISWAPDG
	};

	class Operator {
	public:
		Operator() : type(OperationType::X), qubits(1, 0) {}

		Operator(OperationType type, int q1 = 0, int q2 = 0, int q3 = 0)
			: type(type), qubits(GetNrQubitsForType(type))
		{
			qubits[0] = q1;
			if (GetNrQubits() > 1)
				qubits[1] = q2;
			if (GetNrQubits() > 2)
				qubits[2] = q3;
		}

		int GetNrQubits() const
		{
			return static_cast<int>(qubits.size());
		}

		int GetQubit(size_t index) const
		{
			return qubits[index];
		}

		OperationType GetType() const
		{
			return type;
		}

	private:
		static int GetNrQubitsForType(OperationType type)
		{
			if (static_cast<int>(type) >= static_cast<int>(OperationType::CX))
				return 2;
			
			return 1;
		}

		OperationType type;
		std::vector<int> qubits;
	};

	// the Pauli propagation execution starts with the operator and applies the gates in reverse order on it

	// expectation value is <psi|P|psi> = <0|U^t P U|0>
	// U is the circuit operator, P is the Pauli string here
	// so as U is decomposed in individual gates, U = G_n ... G_2 G_1
	// transforming P is done as P' = G_1^t G_2^t ... G_n^t P G_n ... G_2 G_1
	// so one needs to apply the gates in reverse order

	// more, as the Pauli string classes are originally designed for the generators of the clifford simulator
	// ApplyGate does P' = G P G^t
	// so here we need to apply the ^t version of the gates (in many cases G^t = G, that is, many are hermitian, only for the others care must be taken)
	// see the Apply... implementation for details

	class PauliPropagator {
	public:
		PauliPropagator()
		{
			std::random_device rd;
			rng.seed(rd());
		}

		double ExpectationValue(const std::string& pauliString)
		{
			PauliStringXZWithCoefficient pauli;
			pauli.Resize(nrQubits);

			for (int i = 0; i < pauliString.size(); ++i)
			{
				char c = pauliString[i];
				switch (c)
				{
				case 'X':
				case 'x':
					pauli.X[i] = true;
					break;
				case 'Y':
				case 'y':
					// Y = iXZ
					pauli.X[i] = true;
					pauli.Z[i] = true;
					break;
				case 'Z':
				case 'z':
					pauli.Z[i] = true;
					break;
				}
			}

			return ExpectationValue(pauli);
		}

		double ExpectationValue(const PauliStringXZWithCoefficient& pauliString)
		{
			pauliStringsIn.clear();
			pauliStringsOut.clear();

			PauliStringXZWithCoefficient pstr = pauliString;
			pstr.Resize(nrQubits);

			pauliStringsIn.insert(std::move(pstr));

			Execute();

			return ExpectationValue();
		}

		double ExpectationValue(const std::unordered_set<PauliStringXZWithCoefficient, PauliStringXZHash>& pauliStrings)
		{
			pauliStringsIn.clear();
			pauliStringsOut.clear();

			for (const auto& pstr : pauliStrings)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.Resize(nrQubits);
				pauliStringsIn.insert(std::move(pstrNew));
			}

			Execute();

			return ExpectationValue();
		}

		double Probability0(int qubit)
		{
			PauliStringXZWithCoefficient pstr(nrQubits);
			pstr.Z[qubit] = true;

			return 0.5 * (1.0 + ExpectationValue(pstr));
		}

		double Probability1(int qubit)
		{
			PauliStringXZWithCoefficient pstr(nrQubits);
			pstr.Z[qubit] = true;

			return 0.5 * (1.0 - ExpectationValue(pstr));
		}

		std::vector<bool> Sample(const std::vector<int>& qubits)
		{
			std::vector<bool> results;
			results.reserve(qubits.size());

			std::unordered_set<PauliStringXZWithCoefficient, PauliStringXZHash> pauliStrings;
			PauliStringXZWithCoefficient pauliStr(nrQubits);
			pauliStrings.insert(pauliStr);
			
			// NOTE: the scaling of the coefficients with 0.5 each time is not done, so the expectation values are not normalized
			// the division when computing expecZ below cancels this out

			double expectation = 1.0;
			for (int q = 0; q < qubits.size(); ++q)
			{
				const int qubit = qubits[q];

				std::unordered_set<PauliStringXZWithCoefficient, PauliStringXZHash> newPauliStrings;

				for (const auto& ps : pauliStrings)
				{
					pauliStr = ps;
					pauliStr.Z[qubit] = true;
					newPauliStrings.insert(std::move(pauliStr));
				}

				const double newExpectation = ExpectationValue(newPauliStrings);
				const double expecZ = newExpectation / expectation; // conditional expectation value of Z on qubit given previous measurement result

				// p1 = = <P1> = <(I - Z)/2> = 0.5 * (1 - <Z>)
				const double p1 = 0.5 * (1.0 - expecZ);
				const double prob = uniformZeroOne(rng);
				const bool measuredOne = prob < p1;
				results.push_back(measuredOne);

				if (q == qubits.size() - 1) // no need to go further if this is the last qubit measured
					break;

				// update expectation - see the note above
				if (measuredOne)
					expectation -= newExpectation; // I - Z
				else
					expectation += newExpectation; // I + Z

				for (auto& ps : newPauliStrings)
				{
					if (measuredOne) // also adjust the sign for individual pauli strings if -Z is needed
						ps.Coefficient *= -1.0;

					pauliStrings.insert(ps);
				}
			}

			return results;
		}

		int GetNrQubits() const
		{
			return nrQubits;
		}

		void SetNrQubits(int nQubits)
		{
			nrQubits = nQubits;
		}

		void SaveState()
		{
			pos = operations.size();
		}

		void RestoreState()
		{
			operations.resize(pos);
		}

		void ApplyX(int qubit)
		{
			operations.emplace_back(OperationType::X, static_cast<int>(qubit));
		}

		void ApplyY(int qubit)
		{
			operations.emplace_back(OperationType::Y, static_cast<int>(qubit));
		}

		void ApplyZ(int qubit)
		{
			operations.emplace_back(OperationType::Z, static_cast<int>(qubit));
		}

		void ApplyH(int qubit)
		{
			operations.emplace_back(OperationType::H, static_cast<int>(qubit));
		}

		void ApplyK(int qubit)
		{
			operations.emplace_back(OperationType::K, static_cast<int>(qubit));
		}

		void ApplyS(int qubit)
		{
			operations.emplace_back(OperationType::SDG, static_cast<int>(qubit));
		}

		void ApplySDG(int qubit)
		{
			operations.emplace_back(OperationType::S, static_cast<int>(qubit));
		}

		void ApplySX(int qubit)
		{
			operations.emplace_back(OperationType::SXDG, static_cast<int>(qubit));
		}

		void ApplySXDG(int qubit)
		{
			operations.emplace_back(OperationType::SX, static_cast<int>(qubit));
		}

		void ApplyCX(int control, int target)
		{
			operations.emplace_back(OperationType::CX, static_cast<int>(control), static_cast<int>(target));
		}

		void ApplyCY(int control, int target)
		{
			operations.emplace_back(OperationType::CY, static_cast<int>(control), static_cast<int>(target));
		}

		void ApplyCZ(int control, int target)
		{
			operations.emplace_back(OperationType::CZ, static_cast<int>(control), static_cast<int>(target));
		}

		void ApplySWAP(int qubit1, int qubit2)
		{
			operations.emplace_back(OperationType::SWAP, static_cast<int>(qubit1), static_cast<int>(qubit2));
		}

		void ApplyISWAP(int qubit1, int qubit2)
		{
			operations.emplace_back(OperationType::ISWAPDG, static_cast<int>(qubit1), static_cast<int>(qubit2));
		}

		void ApplyISWAPDG(int qubit1, int qubit2)
		{
			operations.emplace_back(OperationType::ISWAP, static_cast<int>(qubit1), static_cast<int>(qubit2));
		}

	private:
		inline void ExecuteX(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyX(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteY(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyY(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteZ(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyZ(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteH(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyH(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteK(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyK(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteS(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyS(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteSDG(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplySdag(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteSX(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplySx(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteSXDG(int qubit)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplySxDag(static_cast<size_t>(qubit));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteCX(int control, int target)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyCX(static_cast<size_t>(control), static_cast<size_t>(target));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteCY(int control, int target)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyCY(static_cast<size_t>(control), static_cast<size_t>(target));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteCZ(int control, int target)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyCZ(static_cast<size_t>(control), static_cast<size_t>(target));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteSWAP(int qubit1, int qubit2)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplySwap(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteISWAP(int control, int target)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyISwap(static_cast<size_t>(control), static_cast<size_t>(target));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteISWAPDG(int control, int target)
		{
			for (auto& pstr : pauliStringsIn)
			{
				PauliStringXZWithCoefficient pstrNew = pstr;
				pstrNew.ApplyISwapDag(static_cast<size_t>(control), static_cast<size_t>(target));
				Insert(std::move(pstrNew));
			}
		}

		inline void ExecuteTwoQubitOp(OperationType opType, int control, int target)
		{
			switch (opType)
			{
			case OperationType::CX:
				ExecuteCX(control, target);
				break;
			case OperationType::CY:
				ExecuteCY(control, target);
				break;
			case OperationType::CZ:
				ExecuteCZ(control, target);
				break;
			case OperationType::SWAP:
				ExecuteSWAP(control, target);
				break;
			case OperationType::ISWAP:
				ExecuteISWAP(control, target);
				break;
			case OperationType::ISWAPDG:
				ExecuteISWAPDG(control, target);
				break;
			default:
				break;
			}
		}

		void Execute()
		{
			for (int i = static_cast<int>(operations.size()) - 1; i >= 0; --i)
			{
				pauliStringsOut.clear();
				const Operator& op = operations[i];

				switch (op.GetType())
				{
				case OperationType::X:
					ExecuteX(op.GetQubit(0));
					break;
				case OperationType::Y:
					ExecuteY(op.GetQubit(0));
					break;
				case OperationType::Z:
					ExecuteZ(op.GetQubit(0));
					break;
				case OperationType::H:
					ExecuteH(op.GetQubit(0));
					break;
				case OperationType::K:
					ExecuteK(op.GetQubit(0));
					break;
				case OperationType::S:
					ExecuteS(op.GetQubit(0));
					break;
				case OperationType::SDG:
					ExecuteSDG(op.GetQubit(0));
					break;
				case OperationType::SX:
					ExecuteSX(op.GetQubit(0));
					break;
				case OperationType::SXDG:
					ExecuteSXDG(op.GetQubit(0));
					break;
				default:
					ExecuteTwoQubitOp(op.GetType(), op.GetQubit(0), op.GetQubit(1));
					break;
				}

				pauliStringsIn = std::move(pauliStringsOut);
			}

			pauliStringsOut = std::move(pauliStringsIn);
		}

		double ExpectationValue()
		{
			// this is for the future, when for example it will support non-clifford gates (maybe rotation gates)
			// such gates will expand the number of pauli strings
			// so it won't end up with a single pauli string at the end
			// TODO: can be parallelized if there are many pauli strings
			double expValue = 0.0;
			for (const auto& pstr : pauliStringsOut)
				expValue += pstr.ExpectationValue();

			return expValue;
		}

		void Insert(PauliStringXZWithCoefficient&& pstrNew)
		{
			auto it = pauliStringsOut.find(pstrNew);
			if (it != pauliStringsOut.end())
			{
				PauliStringXZWithCoefficient pstrExisting = *it;
				pauliStringsOut.erase(it);
				pstrNew.Coefficient += pstrExisting.Coefficient;
			}

			pauliStringsOut.insert(std::move(pstrNew));
		}

		int nrQubits = 0;
		std::unordered_set<PauliStringXZWithCoefficient, PauliStringXZHash> pauliStringsIn;
		std::unordered_set<PauliStringXZWithCoefficient, PauliStringXZHash> pauliStringsOut;
		std::vector<Operator> operations;
		size_t pos = 0;

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne{ 0., 1. };
	};

}
