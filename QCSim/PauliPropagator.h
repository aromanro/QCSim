#pragma once

#include "PauliPropOper.h"
#include "ThreadPool.h"

#include <random>
#include <unordered_set>
#include <numeric>

namespace QC
{
	struct PauliStringHash
	{
		std::size_t operator()(const PauliStringXZWithCoefficient& p) const
		{
			const auto h1 = std::hash<std::vector<bool>>{}(p.X);
			const auto h2 = std::hash<std::vector<bool>>{}(p.Z);
			return h1 ^ (h2 << 1);
		}
	};

	struct PauliStringEqual
	{
		bool operator()(const PauliStringXZWithCoefficient& a, const PauliStringXZWithCoefficient& b) const
		{
			return (a.X == b.X) && (a.Z == b.Z);
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
			//auto nrThreads = QubitRegisterCalculator<>::GetNumberOfThreads();
			//threadPool = std::make_unique<ThreadPool<>>(nrThreads);
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

		double ExpectationValue(const std::string& pauliString) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliString, pauliStrings);
		}

		double ExpectationValue(const std::string& pauliString, PauliStringStorage& pauliStrings) const
		{
			PauliStringXZWithCoefficient pauli;
			pauli.Resize(nrQubits);

			for (int i = 0; i < static_cast<int>(pauliString.size()); ++i)
			{
				char c = pauliString[i];
				switch (c)
				{
				case 'X': case 'x':
					pauli.X[i] = true;
					break;
				case 'Y': case 'y':
					pauli.X[i] = true;
					pauli.Z[i] = true;
					break;
				case 'Z': case 'z':
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
			PauliStringXZWithCoefficient pstr = pauliString;
			pstr.Resize(nrQubits);
			pauliStrings.push_back(std::move(pstr));
			return Execute(pauliStrings);
		}

		double ExpectationValue(PauliStringXZWithCoefficient&& pauliString) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliString, pauliStrings);
		}

		double ExpectationValue(PauliStringXZWithCoefficient&& pauliString, PauliStringStorage& pauliStrings) const
		{
			pauliStrings.push_back(std::move(pauliString));
			return Execute(pauliStrings);
		}

		double ExpectationValue(const PauliStringStorage& pauliStringsInput) const
		{
			PauliStringStorage pauliStrings;
			return ExpectationValue(pauliStringsInput, pauliStrings);
		}

		
		std::vector<bool> Sample(const std::vector<int>& qubits)
		{
			PauliStringStorage pauliStrings;
			return Sample(qubits, pauliStrings);
		}

		int GetNrQubits() const { return nrQubits; }
		void SetNrQubits(int nQubits) { nrQubits = nQubits; }

		double GetCoefficientThreshold() const { return coefThreshold; }
		void SetCoefficientThreshold(double threshold) { coefThreshold = threshold; }

		size_t GetPauliWeightThreshold() const { return pauliWeightThreshold; }
		void SetPauliWeightThreshold(size_t threshold) { pauliWeightThreshold = threshold; }

		int StepsBetweenTrims() const { return stepsBetweenTrims; }
		void SetStepsBetweenTrims(int steps) { stepsBetweenTrims = steps; }

		int StepsBetweenDeduplication() const { return stepsBetweenDeduplication; }
		void SetStepsBetweenDeduplication(int steps) { stepsBetweenDeduplication = steps; }

		size_t GetParallelThreshold() const { return parallelThreshold; }
		void SetParallelThreshold(size_t threshold) { parallelThreshold = threshold; }

		size_t GetBatchSize() const { return batchSize; }
		void SetBatchSize(size_t size) { batchSize = size; }

		size_t GetParallelThresholdForSum() const { return parallelThresholdSum; }
		void SetParallelThresholdForSum(size_t threshold) { parallelThresholdSum = threshold; }

		size_t GetBatchSizeForSum() const { return batchSizeSum; }
		void SetBatchSizeForSum(size_t size) { batchSizeSum = size; }

		size_t GetSavePosition() const { return pos; }
		void SetSavePosition(size_t position) { pos = position; }

		void EnableParallel(size_t numThreads = 0)
		{
			threadPool = std::make_unique<ThreadPool<>>(numThreads);
		}

		void DisableParallel()
		{
			threadPool.reset();
		}

		bool IsParallelEnabled() const { return threadPool != nullptr; }

		void SaveState() { pos = operations.size(); }
		void RestoreState() { operations.resize(pos); }


		void ApplyX(int qubit) { operations.push_back(std::make_unique<OperatorX>(qubit)); }
		void ApplyY(int qubit) { operations.push_back(std::make_unique<OperatorY>(qubit)); }
		void ApplyZ(int qubit) { operations.push_back(std::make_unique<OperatorZ>(qubit)); }
		void ApplyH(int qubit) { operations.push_back(std::make_unique<OperatorH>(qubit)); }
		void ApplyK(int qubit) { operations.push_back(std::make_unique<OperatorK>(qubit)); }
		void ApplyS(int qubit) { operations.push_back(std::make_unique<OperatorS>(qubit)); }
		void ApplySDG(int qubit) { operations.push_back(std::make_unique<OperatorSDG>(qubit)); }
		void ApplySX(int qubit) { operations.push_back(std::make_unique<OperatorSX>(qubit)); }
		void ApplySXDG(int qubit) { operations.push_back(std::make_unique<OperatorSXDG>(qubit)); }

		void ApplyCX(int control, int target) { operations.push_back(std::make_unique<OperatorCX>(target, control)); }
		void ApplyCY(int control, int target) { operations.push_back(std::make_unique<OperatorCY>(target, control)); }
		void ApplyCZ(int control, int target) { operations.push_back(std::make_unique<OperatorCZ>(target, control)); }

		void ApplySWAP(int qubit1, int qubit2) { operations.push_back(std::make_unique<OperatorSWAP>(qubit1, qubit2)); }
		void ApplyISWAP(int qubit1, int qubit2) { operations.push_back(std::make_unique<OperatorISWAP>(qubit1, qubit2)); }
		void ApplyISWAPDG(int qubit1, int qubit2) { operations.push_back(std::make_unique<OperatorISWAPDG>(qubit1, qubit2)); }

		void ApplyRX(int qubit, double angle) { operations.push_back(std::make_unique<OperatorRX>(qubit, angle)); }
		void ApplyRY(int qubit, double angle) { operations.push_back(std::make_unique<OperatorRY>(qubit, angle)); }
		void ApplyRZ(int qubit, double angle) { operations.push_back(std::make_unique<OperatorRZ>(qubit, angle)); }

		void ClearOperations() { operations.clear(); }

		double Probability(size_t outcome)
		{
			if (nrQubits == 0 || outcome >= 1ULL << nrQubits)
				return 0.0;

			const size_t pos = operations.size();

			double res = 1.0;

			const size_t lastQubit = nrQubits - 1;
			for (size_t q = 0; q < lastQubit; ++q)
			{
				const bool measuredOne = (outcome & 1) == 1;
				const double prob = measuredOne ? Probability1(q) : Probability0(q);
				res *= prob;
				outcome >>= 1;
				std::unique_ptr<Operator> proj = std::make_unique<Projector>(q, measuredOne, 0.5 / prob);
				operations.push_back(std::move(proj));
			}

			const bool measuredOne = (outcome & 1) == 1;
			const double prob = measuredOne ? Probability1(lastQubit) : Probability0(lastQubit);
			res *= prob;
			
			operations.resize(pos);
			
			return res;
		}

		std::vector<std::unique_ptr<Operator>> GetOperations() const
		{
			std::vector<std::unique_ptr<Operator>> ops;
			ops.reserve(operations.size());
			for (const auto& op : operations)
				ops.push_back(std::move(op->Clone()));
			return ops;
		}

		void SetOperations(std::vector<std::unique_ptr<Operator>>&& newOps)
		{
			operations = std::move(newOps);
		}

	private:
		double ExpectationValue(const PauliStringStorage& pauliStringsInput, PauliStringStorage& pauliStrings) const
		{
			pauliStrings.reserve(pauliStringsInput.size());
			for (const auto& ps : pauliStringsInput)
				pauliStrings.push_back(ps);
			return Execute(pauliStrings);
		}

		std::vector<bool> Measure(const std::vector<int>& qubits, PauliStringStorage& pauliStrings)
		{
			std::vector<bool> results;
			results.reserve(qubits.size());

			for (int qubit : qubits)
			{
				const double p1 = Probability1(qubit, pauliStrings);
				pauliStrings.clear();
				const double prob = uniformZeroOne(rng);
				const bool measuredOne = prob < p1;
				results.push_back(measuredOne);

				std::unique_ptr<Operator> proj = std::make_unique<Projector>(qubit, measuredOne, 0.5 / (measuredOne ? p1 : 1. - p1));
				operations.push_back(std::move(proj));
			}

			return results;
		}

		std::vector<bool> Sample(const std::vector<int>& qubits, PauliStringStorage& pauliStrings)
		{
			std::vector<bool> results;
			results.reserve(qubits.size());

			PauliStringStorage pauliStringsStart;

			PauliStringXZWithCoefficient pauliStr(nrQubits);
			pauliStringsStart.push_back(std::move(pauliStr));

			double expectation = 1.0;
			for (int q = 0; q < static_cast<int>(qubits.size()); ++q)
			{
				const int qubit = qubits[q];
				
				PauliStringStorage newPauliStrings;
				newPauliStrings.reserve(pauliStringsStart.size());

				for (const auto& ps : pauliStringsStart)
				{
					pauliStr = PauliStringXZWithCoefficient(ps);
					pauliStr.Z[qubit] = true;
					newPauliStrings.push_back(std::move(pauliStr));
				}

				const double newExpectation = ExpectationValue(newPauliStrings, pauliStrings);
				pauliStrings.clear();
				const double expecZ = newExpectation / expectation;

				const double p1 = 0.5 * (1.0 - expecZ);
				const double prob = uniformZeroOne(rng);
				const bool measuredOne = prob < p1;
				results.push_back(measuredOne);

				if (q == static_cast<int>(qubits.size()) - 1)
					break;

				if (measuredOne)
					expectation -= newExpectation;
				else
					expectation += newExpectation;

				if (pauliWeightThreshold < static_cast<size_t>(nrQubits) && q % stepsBetweenTrims == 0)
					TrimWithoutDeduplication(newPauliStrings);

				pauliStringsStart.reserve(pauliStringsStart.size() + newPauliStrings.size());
				for (auto& ps : newPauliStrings)
				{
					if (measuredOne)
						ps.Coefficient *= -1.0;
					pauliStringsStart.push_back(std::move(ps));
				}
			}

			return results;
		}

		inline double ExecuteSequential(PauliStringStorage& pauliStrings) const
		{
			for (int i = static_cast<int>(operations.size()) - 1; i >= 0; --i)
			{
				const auto& op = operations[i];
				const size_t sizeBefore = pauliStrings.size();

				if (op->GetType() >= OperationType::PROJ)
					pauliStrings.reserve(sizeBefore * 2);

				for (size_t j = 0; j < sizeBefore; ++j)
					op->Apply(pauliStrings[j], pauliStrings);

				if (stepsBetweenDeduplication != std::numeric_limits<int>::max() && i % stepsBetweenDeduplication == 0)
					TrimWithDeduplication(pauliStrings);
				else if (stepsBetweenTrims != std::numeric_limits<int>::max() && i % stepsBetweenTrims == 0)
					TrimWithoutDeduplication(pauliStrings);
			}

			return ExpectationValueAtEnd(pauliStrings);
		}


		void TrimWithDeduplication(PauliStringStorage& pauliStrings) const
		{
			std::unordered_set<PauliStringXZWithCoefficient, PauliStringHash, PauliStringEqual> uniqueSet;
			for (auto& ps : pauliStrings)
			{
				if (pauliWeightThreshold < static_cast<size_t>(nrQubits) && ps.PauliWeight() > pauliWeightThreshold)
					continue;

				auto it = uniqueSet.find(ps);
				if (it != uniqueSet.end())
					(*it).Coefficient += ps.Coefficient;
				else
					uniqueSet.insert(std::move(ps));
			}
			pauliStrings.clear();
			pauliStrings.reserve(uniqueSet.size());

			for (auto& ps : uniqueSet)
			{
				if (std::abs(ps.Coefficient) > coefThreshold)
					pauliStrings.push_back(std::move(const_cast<PauliStringXZWithCoefficient&>(ps)));
			}
		}

		void TrimWithoutDeduplication(PauliStringStorage& pauliStrings) const
		{
			size_t writePos = 0;
			for (size_t readPos = 0; readPos < pauliStrings.size(); ++readPos)
			{
				auto& ps = pauliStrings[readPos];
				if (std::abs(ps.Coefficient) > coefThreshold && (pauliWeightThreshold >= static_cast<size_t>(nrQubits) || ps.PauliWeight() <= pauliWeightThreshold))
				{
					if (writePos != readPos)
						pauliStrings[writePos] = std::move(pauliStrings[readPos]);
					++writePos;
				}
			}
			pauliStrings.resize(writePos);
		}

		inline void ApplyOperatorsParallel(int opIndex, PauliStringStorage& pauliStrings, std::mutex& resmtx, std::vector<std::future<double>>& futures, int& inExec, std::condition_variable& cv) const
		{
			const size_t sizeBefore = pauliStrings.size();
			const size_t nJobs = (sizeBefore + batchSize - 1) / batchSize;

			const auto& op = operations[opIndex];

			size_t start = 0;
			for (size_t job = 0; job < nJobs; ++job)
			{
				const size_t end = std::min(start + batchSize, sizeBefore);

				PauliStringStorage localNewStrings;
				localNewStrings.reserve(end - start);

				for (size_t j = start; j < end; ++j)
					op->Apply(pauliStrings[j], localNewStrings);

				if (!localNewStrings.empty())
				{
					{
						std::lock_guard<std::mutex> lock(resmtx);
						++inExec;
					}
					auto futureJob = threadPool->Enqueue(
						[this, opIndex, localStrings = std::move(localNewStrings), &resmtx, &futures, &inExec, &cv]() -> double
						{
							auto ps = std::move(localStrings);
							const double res = ExecuteParallel(ps, opIndex - 1, resmtx, futures, inExec, cv);

							{
								std::lock_guard<std::mutex> lock(resmtx);
								--inExec;
							}
							cv.notify_one();
							return res;
						}
					);

					std::lock_guard<std::mutex> lock(resmtx);
					futures.push_back(std::move(futureJob));
				}
				start = end;
			}
		}

		inline double ExecuteParallel(PauliStringStorage& pauliStrings, int startIndex, std::mutex& resmtx, std::vector<std::future<double>>& futures, int& inExec, std::condition_variable& cv) const
		{
			for (int i = startIndex; i >= 0; --i)
			{
				const size_t sizeBefore = pauliStrings.size();

				// Only go parallel when there are enough strings to justify it
				const auto& op = operations[i];
				const bool isSplitting = op->GetType() >= OperationType::PROJ;
				if (sizeBefore >= parallelThreshold && isSplitting)
				{
					ApplyOperatorsParallel(i, pauliStrings, resmtx, futures, inExec, cv);
				}
				else
				{
					if (isSplitting)
						pauliStrings.reserve(sizeBefore * 2);

					for (size_t j = 0; j < sizeBefore; ++j)
						op->Apply(pauliStrings[j], pauliStrings);
				}

				if (stepsBetweenDeduplication != std::numeric_limits<int>::max() && i % stepsBetweenDeduplication == 0)
					TrimWithDeduplication(pauliStrings);
				else if (stepsBetweenTrims != std::numeric_limits<int>::max() && i % stepsBetweenTrims == 0)
					TrimWithoutDeduplication(pauliStrings);
			}

			const double res = ExpectationValueAtEnd(pauliStrings);

			return res;
		}

		inline double Execute(PauliStringStorage& pauliStrings) const
		{
			double res;
			if (threadPool)
			{
				std::vector<std::future<double>> futures;
				std::mutex resmtx;
				int inExec = 0;
				std::condition_variable cv;
				res = ExecuteParallel(pauliStrings, static_cast<int>(operations.size()) - 1, resmtx, futures, inExec, cv);

				for (;;)
				{
					std::unique_lock<std::mutex> lock(resmtx);

					cv.wait(lock, [&inExec] { return inExec == 0; });
					if (inExec == 0)
						break;
					std::this_thread::yield();
				}

				for (auto& f : futures)
					res += f.get();
			}
			else
				res = ExecuteSequential(pauliStrings);

			return res;
		}

		inline double ExpectationValueAtEnd(const PauliStringStorage& pauliStrings) const
		{
			const size_t size = pauliStrings.size();

			if (!threadPool || size < parallelThresholdSum)
			{
				double expValue = 0.0;
				for (const auto& pstr : pauliStrings)
					expValue += pstr.ExpectationValue();
				return expValue;
			}

			const size_t nJobs = (size + batchSizeSum - 1) / batchSizeSum;

			std::vector<std::future<double>> futures;
			futures.reserve(nJobs);

			size_t start = 0;
			for (size_t job = 0; job < nJobs; ++job)
			{
				const size_t end = std::min(start + batchSizeSum, size);

				futures.push_back(threadPool->Enqueue(
					[&pauliStrings, start, end]() -> double
					{
						double sum = 0.0;
						for (size_t i = start; i < end; ++i)
							sum += pauliStrings[i].ExpectationValue();
						return sum;
					}
				));
				start = end;
			}

			double result = 0.0;
			for (auto& f : futures)
				result += f.get();

			return result;
		}

		int nrQubits = 0;
		int stepsBetweenTrims = std::numeric_limits<int>::max();
		int stepsBetweenDeduplication = std::numeric_limits<int>::max();

		std::vector<std::unique_ptr<Operator>> operations;
		size_t pos = 0;

		double coefThreshold = 0.;
		size_t pauliWeightThreshold = std::numeric_limits<size_t>::max();

		std::mt19937_64 rng;
		std::uniform_real_distribution<double> uniformZeroOne{ 0., 1. };

		mutable std::unique_ptr<ThreadPool<>> threadPool;
		size_t parallelThresholdSum = 4096;
		size_t batchSizeSum = 2048; 
		size_t parallelThreshold = 256;
		size_t batchSize = 128;
	};
}
