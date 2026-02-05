#pragma once

#include "PauliPropOper.h"

#include <random>
#include <unordered_set>

namespace QC
{
	using PauliStringPtr = std::unique_ptr<PauliStringXZWithCoefficient>;
	using PauliStringStorage = std::vector<PauliStringPtr>;

	struct PauliStringHash
	{
		std::size_t operator()(const PauliStringPtr& p) const
		{
			const auto h1 = std::hash<std::vector<bool>>{}(p->X);
			const auto h2 = std::hash<std::vector<bool>>{}(p->Z);
			return h1 ^ (h2 << 1);
		}
	};

	struct PauliStringEqual
	{
		bool operator()(const PauliStringPtr& a, const PauliStringPtr& b) const
		{
			return (a->X == b->X) && (a->Z == b->Z);
		}
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
			std::random_device rdev;
			rng.seed(rdev());
		}

		double Probability0(int qubit) const
		{
			PauliStringStorage pauliStrings;

			return Probability0(qubit, pauliStrings);
		}

		double Probability0(int qubit, PauliStringStorage& pauliStrings) const
		{
			PauliStringXZWithCoefficient pstr(nrQubits);
			pstr.Z[qubit] = true;

			return 0.5 * (1.0 + ExpectationValue(std::move(pstr), pauliStrings));
		}

		double Probability1(int qubit) const
		{
			PauliStringStorage pauliStrings;
			return Probability1(qubit, pauliStrings);
		}

		double Probability1(int qubit, PauliStringStorage& pauliStrings) const
		{
			PauliStringXZWithCoefficient pstr(nrQubits);
			pstr.Z[qubit] = true;

			return 0.5 * (1.0 - ExpectationValue(std::move(pstr), pauliStrings));
		}

		std::vector<bool> Measure(const std::vector<int>& qubits)
		{
			PauliStringStorage pauliStrings;
			return Measure(qubits, pauliStrings);
		}

		std::vector<bool> Measure(const std::vector<int>& qubits, PauliStringStorage& pauliStrings)
		{
			std::vector<bool> results;
			results.reserve(qubits.size());

			for (int qubit : qubits)
			{
				pauliStrings.clear();
				const double p1 = Probability1(qubit, pauliStrings);
				const double prob = uniformZeroOne(rng);
				const bool measuredOne = prob < p1;
				results.push_back(measuredOne);

				// update the circuit with the projector
				// the projector is P0 = (I + Z)/2 or P1 = (I - Z)/2 depending on what was measured

				// its action on pauli strings is:
				// Pout = 1 / (4 * probability) * (I +/- Z) Pin (I +/- Z)

				// there are two cases:
				// Pin commutes with Z (Z or I is present on that qubit)
				// Pin anticommutes with Z (X or Y is present on that qubit) - it turns out that in this case the result is 0

				// { sigma_i, sigma_j } = 2 delta_ij
				// also { Z, I } = 2 Z appears

				// the whole expression is:

				// Pout = 1 / (4 * probability) * (Pin +/- { Pin, Z } + Z Pin Z)
				// Z Pin Z doesn't do anything to the pauli string except changing the sign if X or Y is present on that qubit
				// in that case the last term cancels the first and the anticommutator is 0, so the result is 0

				// for the commuting case:
				// Pout = 1 / (4 * probability) * (2 * Pin +/- { Pin, Z })
				// { Pin, Z } = 2 Pmod, where Pmod is Pin with Z changed to I if present on that qubit, or I changed to Z if present

				// see the implementation for the execution of the projector for the details

				std::unique_ptr<Operator> proj = std::make_unique<Projector>(qubit, measuredOne, 0.5 / (measuredOne ? p1 : 1. - p1));
				operations.push_back(std::move(proj));
			}

			return results;
		}

		double ExpectationValue(const std::string& pauliString) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliString, pauliStrings);
		}

		double ExpectationValue(const std::string& pauliString, PauliStringStorage& pauliStrings) const
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

			return ExpectationValue(pauli, pauliStrings);
		}

		double ExpectationValue(const PauliStringXZWithCoefficient& pauliString) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliString, pauliStrings);
		}

		double ExpectationValue(const PauliStringXZWithCoefficient& pauliString, PauliStringStorage& pauliStrings) const
		{
			PauliStringPtr pstr = std::make_unique<PauliStringXZWithCoefficient>(pauliString);
			pstr->Resize(nrQubits);

			pauliStrings.push_back(std::move(pstr));

			Execute(pauliStrings);

			return ExpectationValueAtEnd(pauliStrings);
		}

		double ExpectationValue(PauliStringXZWithCoefficient&& pauliString) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliString, pauliStrings);
		}

		double ExpectationValue(PauliStringXZWithCoefficient&& pauliString, PauliStringStorage& pauliStrings) const
		{
			PauliStringPtr pstr = std::make_unique<PauliStringXZWithCoefficient>(std::move(pauliString));
			pauliStrings.push_back(std::move(pstr));

			Execute(pauliStrings);

			return ExpectationValueAtEnd(pauliStrings);
		}

		double ExpectationValue(const PauliStringStorage& pauliStringsInput) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliStringsInput, pauliStrings);
		}

		double ExpectationValue(const PauliStringStorage& pauliStringsInput, PauliStringStorage& pauliStrings) const
		{
			pauliStrings.reserve(pauliStringsInput.size());
			for (const auto& ps : pauliStringsInput)
			{
				PauliStringPtr pstr = std::make_unique<PauliStringXZWithCoefficient>(*ps);
				pauliStrings.push_back(std::move(pstr));
			}
			Execute(pauliStrings);
			return ExpectationValueAtEnd(pauliStrings);
		}

		std::vector<bool> Sample(const std::vector<int>& qubits)
		{
			PauliStringStorage pauliStrings;
			return Sample(qubits, pauliStrings);
		}

		std::vector<bool> Sample(const std::vector<int>& qubits, PauliStringStorage& pauliStrings)
		{
			std::vector<bool> results;
			results.reserve(qubits.size());

			PauliStringStorage pauliStringsStart;
			//pauliStringsStart.reserve(1ULL << qubits.size());
			
			// start with the identity on all qubits
			auto pauliStr = std::make_unique<PauliStringXZWithCoefficient>(nrQubits);
			pauliStringsStart.push_back(std::move(pauliStr));

			// NOTE: the scaling of the coefficients with 0.5 each time is not done, so the expectation values are not normalized
			// the division when computing expecZ below cancels this out

			double expectation = 1.0;
			for (int q = 0; q < qubits.size(); ++q)
			{
				const int qubit = qubits[q];
				pauliStrings.clear();

				PauliStringStorage newPauliStrings;
				newPauliStrings.reserve(pauliStringsStart.size());

				for (const auto& ps : pauliStringsStart)
				{
					pauliStr = std::make_unique<PauliStringXZWithCoefficient>(*ps);
					pauliStr->Z[qubit] = true;
					newPauliStrings.push_back(std::move(pauliStr));
				}

				const double newExpectation = ExpectationValue(newPauliStrings, pauliStrings);
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

				pauliStringsStart.reserve(pauliStringsStart.size() * 2);
				for (auto& ps : newPauliStrings)
				{
					if (measuredOne) // also adjust the sign for individual pauli strings if -Z is needed
						ps->Coefficient *= -1.0;

					pauliStringsStart.push_back(std::move(ps));
				}
			}

			return results;
		}

		void TrimWithDeduplication(PauliStringStorage& pauliStrings) const
		{
			// remove duplicate pauli strings
			std::unordered_set<PauliStringPtr, PauliStringHash, PauliStringEqual> uniqueSet;
			for (auto& ps : pauliStrings)
			{
				// also remove the ones with pauli weight over a threshold
				if (pauliWeightThreshold < nrQubits &&ps->PauliWeight() > pauliWeightThreshold)
					continue;

				auto it = uniqueSet.find(ps);
				if (it != uniqueSet.end())
					(*it)->Coefficient += ps->Coefficient;
				else
					uniqueSet.insert(std::move(ps));
			}
			pauliStrings.clear();
			pauliStrings.reserve(uniqueSet.size());
			
			for (auto& ps : uniqueSet)
			{
				if (std::abs(ps->Coefficient) > coefThreshold) // remove small coefficient pauli strings
					pauliStrings.push_back(std::move(const_cast<PauliStringPtr&>(ps)));
			}
		}

		void TrimWithoutDeduplication(PauliStringStorage& pauliStrings) const
		{
			// use the two pointer technique to remove small coefficient pauli strings in place
			// and also the ones with pauli weight over a threshold
			size_t writePos = 0;
			for (size_t readPos = 0; readPos < pauliStrings.size(); ++readPos)
			{
				auto& ps = pauliStrings[readPos];
				if (std::abs(ps->Coefficient) > coefThreshold && (pauliWeightThreshold >= nrQubits || ps->PauliWeight() <= pauliWeightThreshold))
				{
					if (writePos != readPos)
						pauliStrings[writePos] = std::move(pauliStrings[readPos]);
					++writePos;
				}
			}
			pauliStrings.resize(writePos);
		}

		int GetNrQubits() const
		{
			return nrQubits;
		}

		void SetNrQubits(int nQubits)
		{
			nrQubits = nQubits;
		}

		double GetCoefficientThreshold() const
		{
			return coefThreshold;
		}

		void SetCoefficientThreshold(double threshold)
		{
			coefThreshold = threshold;
		}

		size_t GetPauliWeightThreshold() const
		{
			return pauliWeightThreshold;
		}
		
		void SetPauliWeightThreshold(size_t threshold)
		{
			pauliWeightThreshold = threshold;
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
			std::unique_ptr<Operator> op = std::make_unique<OperatorX>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplyY(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorY>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplyZ(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorZ>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplyH(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorH>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplyK(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorK>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplyS(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorS>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplySDG(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorSDG>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplySX(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorSX>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplySXDG(int qubit)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorSXDG>(qubit);
			operations.push_back(std::move(op));
		}

		void ApplyCX(int target, int control)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorCX>(target, control);
			operations.push_back(std::move(op));
		}

		void ApplyCY(int target, int control)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorCY>(target, control);
			operations.push_back(std::move(op));
		}

		void ApplyCZ(int target, int control)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorCZ>(target, control);
			operations.push_back(std::move(op));
		}

		void ApplySWAP(int qubit1, int qubit2)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorSWAP>(qubit1, qubit2);
			operations.push_back(std::move(op));
		}

		void ApplyISWAP(int qubit1, int qubit2)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorISWAP>(qubit1, qubit2);
			operations.push_back(std::move(op));
		}

		void ApplyISWAPDG(int qubit1, int qubit2)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorISWAPDG>(qubit1, qubit2);
			operations.push_back(std::move(op));
		}

		void ApplyRX(int qubit, double angle)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorRX>(qubit, angle);
			operations.push_back(std::move(op));
		}

		void ApplyRY(int qubit, double angle)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorRY>(qubit, angle);
			operations.push_back(std::move(op));
		}

		void ApplyRZ(int qubit, double angle)
		{
			std::unique_ptr<Operator> op = std::make_unique<OperatorRZ>(qubit, angle);
			operations.push_back(std::move(op));
		}

	private:
		inline void Execute(PauliStringStorage& pauliStrings) const
		{
			for (int i = static_cast<int>(operations.size()) - 1; i >= 0; --i)
			{
				const auto& op = operations[i];
				// the operator needs to be applied on all current pauli strings
				// but not on the newly created ones during this operation - a projector or a non-clifford gate may create new pauli strings
				const size_t sizeBefore = pauliStrings.size();

				if (op->GetType() >= OperationType::PROJ)
					pauliStrings.reserve(sizeBefore * 2); // in case of projector, can double the number of pauli strings

				for (size_t j = 0; j < sizeBefore; ++j)
					op->Apply(pauliStrings[j], pauliStrings);
			}
		}

		static inline double ExpectationValueAtEnd(const PauliStringStorage& pauliStrings)
		{
			// this is for the future, when for example it will support non-clifford gates (maybe rotation gates)
			// such gates will expand the number of pauli strings
			// so it won't end up with a single pauli string at the end
			
			// also it's useful for measuring expectation values of sums of pauli strings directly
			// and for sampling too
			
			// TODO: can be parallelized if there are many pauli strings
			double expValue = 0.0;
			for (const auto& pstr : pauliStrings)
				expValue += pstr->ExpectationValue();

			return expValue;
		}

		int nrQubits = 0;

		std::vector<std::unique_ptr<Operator>> operations;
		size_t pos = 0;

		double coefThreshold = 0.0;
		size_t pauliWeightThreshold = std::numeric_limits<size_t>::max();

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne{ 0., 1. };
	};

}
