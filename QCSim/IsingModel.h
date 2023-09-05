#pragma once

#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

// roughly based on Lesson 10 from
// "Fundamentals In Quantum Algorithms: A Tutorial Series Using Qiskit Continued" by Daniel Koch, Saahil Patel, Laura Wessing, Paul M. Alsing
// https://arxiv.org/abs/2008.10647

// but the code there seem to be wrong (they forgot to take into account the magnetic field values), so a good part of the chapter is based on wrong results

namespace Models {

	template <typename type1 = size_t, typename type2 = size_t> class PairHash {
	public:
		size_t operator()(const std::pair<type1, type2>& t) const
		{
			return 31 * t.first + t.second;
		}
	};

	class IsingModel
	{
	public:
		struct Interaction {
			size_t i = 0;
			size_t j = 0;
			double J = 1.;
		};

		void Clear()
		{
			spins.clear();
			h.clear();
			neighbours.clear();
			interactions.clear();
		}

		void AddInteraction(size_t i, size_t j, double J)
		{
			if (i == j) return;
			else if (i > j) std::swap(i, j);

			if (j >= spins.size()) return;

			neighbours[i].insert(j);
			interactions[std::make_pair(i, j)] = J;
		}

		void Set(const std::vector<bool>& s, const std::vector<Interaction>& ints = {}, const std::vector<double>& hvals = {})
		{
			Clear();
			spins = s;
			h = hvals;
			h.resize(s.size(), 0.);

			for (const auto& interaction : ints)
				AddInteraction(interaction.i, interaction.j, interaction.J);
		}

		void Set(const std::vector<Interaction>& ints = {}, const std::vector<double>& hvals = {})
		{
			Clear();

			h = hvals;
			spins.resize(h.size(), false);

			for (const auto& interaction : ints)
				AddInteraction(interaction.i, interaction.j, interaction.J);
		}

		void SetH(size_t i, double hval)
		{
			if (i >= spins.size()) return;

			h[i] = hval;
		}
		
		double GetH(size_t i) const
		{
			if (h.size() <= i) return 0.;

			return h[i];
		}

		void SetSpin(size_t i, bool s)
		{
			if (i >= spins.size()) return;

			spins[i] = s;
		}

		bool GetSpin(size_t i) const
		{
			if (i >= spins.size()) return false;

			return spins[i];
		}

		double Energy()
		{
			double energy = 0;

			for (size_t i = 0; i < spins.size(); ++i)
			{
				energy += GetH(i) * GetTrueSpin(spins[i]);

				for (const auto& j : neighbours[i])
					energy += interactions[std::make_pair(i, j)] * GetTrueSpin(spins[i]) * GetTrueSpin(spins[j]);
			}

			return -energy;
		}

		void SetState(unsigned int state)
		{
			for (size_t pos = 0; pos < spins.size(); ++pos)
			{
				spins[pos] = (state & 1) ? true : false;
				state >>= 1;
			}
		}

		double Energy(unsigned int state)
		{
			SetState(state);
			return Energy();
		}

		unsigned int GetMinEnergyState()
		{
			unsigned int maxStates = 1 << spins.size();

			double minEnergy = std::numeric_limits<double>::max();
			unsigned int minState = 0;

			for (unsigned int state = 0; state < maxStates; ++state)
			{
				const double stateEnergy = Energy(state);
				if (minEnergy > stateEnergy)
				{
					minEnergy = stateEnergy;
					minState = state;
				}
			}
			
			return minState;
		}

		std::set<unsigned int> GetMinEnergyStates()
		{
			std::set<unsigned int> res;
			unsigned int maxStates = 1 << spins.size();

			double minEnergy = std::numeric_limits<double>::max();
			unsigned int minState = 0;

			for (unsigned int state = 0; state < maxStates; ++state)
			{
				const double stateEnergy = Energy(state);

				if (minEnergy > stateEnergy && abs(minEnergy - stateEnergy) > 1E-9)
				{
					minEnergy = stateEnergy;
					minState = state;
					res.clear();
					res.insert(minState);
				}
				else if (abs(minEnergy - stateEnergy) <= 1E-9)
					res.insert(state);
			}

			return res;
		}

		const std::vector<bool>& GetSpins() const
		{
			return spins;
		}

		const std::vector<double>& GetH() const
		{
			return h;
		}

		const std::unordered_map<size_t, std::unordered_set<size_t>>& GetNeighbours() const
		{
			return neighbours;
		}

		const std::unordered_map<std::pair<size_t, size_t>, double, PairHash<>>& GetInteractions() const
		{
			return interactions;
		}

	protected:
		static int GetTrueSpin(bool s)
		{
			return s ? -1 : 1;
		}

		std::vector<bool> spins;
		std::vector<double> h;
		std::unordered_map<size_t, std::unordered_set<size_t>> neighbours;
		std::unordered_map<std::pair<size_t, size_t>, double, PairHash<>> interactions;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class QAOAIsingSubalgorithm : public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using RegisterClass = QC::QubitRegister<VectorClass, MatrixClass>;

		unsigned int Execute(RegisterClass& reg) override
		{
			double Eold = std::numeric_limits<double>::max();
			double Emin = Eold;
			double E = 0; // whatever
			
			double gamma1 = gamma1Start;
			double beta1 = beta1Start;
			double gamma2 = gamma2Start;
			double beta2 = beta2Start;

			double optGamma1 = gamma1;
			double optBeta1 = beta1;
			double optGamma2 = gamma2;
			double optBeta2 = beta2;

			for (int i = 0; abs(E - Eold) > deltaE && i < 100000; ++i)
			{
				Eold = E;

				// this is more like stochastic gradient descent, so it would benefit from
				// the methods used in the machine learning project (momentum/nesterov/adagrad/rmsprop/adam/whatever) 
				// I won't bother, for the curious the methods are implemented there (also in the python repo, the dft notebook)

				std::tie(gamma1, beta1, gamma2, beta2) = GradientDescentStep(reg, gamma1, beta1, gamma2, beta2, pVal, epsilon, stepSize, nrMeasurements);

				Exec(reg, gamma1, beta1, gamma2, beta2, pVal);
				E = EnergyExpectationValue(reg);

				if (E < Emin)
				{
					Emin = E;
					optGamma1 = gamma1;
					optBeta1 = beta1;
					optGamma2 = gamma2;
					optBeta2 = beta2;
				}

				//if (i % 10 == 0)
					std::cout << "E: " << E << " gamma1: " << gamma1 << " beta1: " << beta1 << " gamma2: " << gamma2 << " beta2: " << beta2 << std::endl;
			}

			Exec(reg, optGamma1, optBeta1, optGamma2, optBeta2, pVal);
			auto res = reg.RepeatedMeasure(nrMeasurements);

			Emin = std::numeric_limits<double>::max();
			unsigned int state = 0;
			for (const auto& r : res)
			{
				E = Energy(r.first) * static_cast<double>(r.second) / nrMeasurements;
				if (E < Emin)
				{
					Emin = E;
					state = r.first;
				}
			}

			return state;
		}

		void Clear()
		{
			model.Clear();
		}

		void AddInteraction(size_t i, size_t j, double J)
		{
			model.AddInteraction(i, j, J);
		}

		void Set(const std::vector<bool>& s, const std::vector<IsingModel::Interaction>& ints = {}, const std::vector<double>& hvals = {})
		{
			model.Set(s, ints, hvals);
		}

		void Set(const std::vector<IsingModel::Interaction>& ints = {}, const std::vector<double>& hvals = {})
		{
			model.Set(ints, hvals);
		}

		void SetH(size_t i, double hval)
		{
			model.SetH(i, hval);
		}

		double GetH(size_t i) const
		{
			return model.GetH(i);
		}

		void SetSpin(size_t i, bool s)
		{
			model.SetSpin(i, s);
		}

		bool GetSpin(size_t i) const
		{
			return model.GetSpin(i);
		}

		double Energy()
		{
			return model.Energy();
		}

		void SetState(unsigned int state)
		{
			model.SetState(state);
		}

		double Energy(unsigned int state)
		{
			return model.Energy(state);
		}

		unsigned int GetMinEnergyState()
		{
			return model.GetMinEnergyState();
		}

		std::set<unsigned int> GetMinEnergyStates()
		{
			return model.GetMinEnergyStates();
		}

		const std::vector<bool>& GetSpins() const
		{
			return model.GetSpins();
		}

		const std::vector<double>& GetH() const
		{
			return model.GetH();
		}

		const std::unordered_map<size_t, std::unordered_set<size_t>>& GetNeighbours() const
		{
			return model.GetNeighbours();
		}

		const std::unordered_map<std::pair<size_t, size_t>, double, PairHash<>>& GetInteractions() const
		{
			return model.GetInteractions();
		}

		double GetDeltaE() const
		{
			return deltaE;
		}

		void SetDeltaE(double dE)
		{
			deltaE = dE;
		}

		double GetGamma1Start() const
		{
			return gamma1Start;
		}

		void SetGamma1Start(double gamma)
		{
			gamma1Start = gamma;
		}

		double GetBeta1Start() const
		{
			return beta1Start;
		}

		void SetBeta1Start(double beta)
		{
			beta1Start = beta;
		}

		double GetGamma2Start() const
		{
			return gamma2Start;
		}

		void SetGamma2Start(double gamma)
		{
			gamma2Start = gamma;
		}

		double GetBeta2Start() const
		{
			return beta2Start;
		}

		void SetBeta2Start(double beta)
		{
			beta2Start = beta;
		}

		unsigned int GetP() const
		{
			return pVal;
		}

		void SetP(unsigned int p)
		{
			if (p != 1 && p != 2) return;

			pVal = p;
		}

		double GetEpsilon() const
		{
			return epsilon;
		}

		void SetEpsilon(double eps)
		{
			epsilon = eps;
		}

		double GetStepSize() const
		{
			return stepSize;
		}

		void SetStepSize(double step)
		{
			stepSize = step;
		}

		double GetNrMeasurements() const
		{
			return nrMeasurements;
		}

		void SetNrMeasurements(double nr)
		{
			nrMeasurements = nr;
		}

		void SetMixing(bool better)
		{
			betterMixing = better;
		}

		bool GetMixing() const
		{
			return betterMixing;
		}

	protected:
		void ApplyHadamardOnAll(RegisterClass& reg)
		{
			for (size_t q = 0; q < reg.getNrQubits(); ++q)
				reg.ApplyGate(h, q);
		}

		void ApplyIsingOperator(RegisterClass& reg, double gamma)
		{
			// interactions
			const auto& neighbors = model.GetNeighbours();
			const auto& interactions = model.GetInteractions();

			for (const auto& site : neighbors)
				for (const auto& n : site.second)
				{
					const unsigned int q1 = site.first;
					const unsigned int q2 = n;
					const auto& J = interactions.find({ q1, q2 });
					if (J != interactions.end())
					{
						const double Jval = J->second;

						rz.SetTheta(2 * gamma * Jval);
					}
					else
						rz.SetTheta(2 * gamma);

					reg.ApplyGate(cnot, q1, q2);
					reg.ApplyGate(rz, q1);
					reg.ApplyGate(cnot, q1, q2);
				}

			// on site
			for (size_t q = 0; q < reg.getNrQubits(); ++q)
			{
				const double hval = model.GetH(q);

				// no need to apply the gate if hval == 0
				// it's identity
				// as a note, hval can be 0 for all spins in some examples
				// for example for max/min cut problems
				if (hval == 0.) continue;

				rz.SetTheta(gamma * hval);
				reg.ApplyGate(rz, q);
			}
		}

		void ApplyMixingOperator(RegisterClass& reg, double beta, bool mixMore = false)
		{
			const double twobeta = 2 * beta;
			rx.SetTheta(twobeta);

			for (size_t q = 0; q < reg.getNrQubits(); ++q)
				reg.ApplyGate(rx, q);

			
			if (mixMore)
			{
				ry.SetTheta(twobeta);

				unsigned int lastQubit = reg.getNrQubits() - 1;
				for (size_t q = 0; q < lastQubit; ++q)
					reg.ApplyGate(cnot, q + 1, q);
				reg.ApplyGate(cnot, 0, lastQubit);

				for (size_t q = 0; q < reg.getNrQubits(); ++q)
					reg.ApplyGate(ry, q);
			}
		}


		void Exec(RegisterClass& reg, double gamma1, double beta1, double gamma2, double beta2, unsigned int p = 1)
		{
			reg.setToBasisState(0);
			ApplyHadamardOnAll(reg);

			for (unsigned int i = 0; i < p; ++i)
			{
				ApplyIsingOperator(reg, i ? gamma2 : gamma1);
				ApplyMixingOperator(reg, i ? beta2 : beta1, betterMixing);
			}
		}

		double EnergyExpectationValue(RegisterClass& reg, unsigned int nrShots = 100000)
		{
			double energy = 0;

			auto res = reg.RepeatedMeasure(nrShots);
			for (const auto& r : res)
				energy += Energy(r.first) * static_cast<double>(r.second) / nrShots;

			return energy;
		}

		std::tuple<double, double, double, double> GradientDescentStep(RegisterClass& reg, double gamma1, double beta1, double gamma2, double beta2, unsigned int p = 1, double eps = 0.0002, double step = 0.0001, unsigned int nrShots = 100000)
		{
			const double twoEps = 2. * eps;

			Exec(reg, gamma1 - eps, beta1, gamma2, beta2, p);
			const double E1 = EnergyExpectationValue(reg, nrShots);

			Exec(reg, gamma1 + eps, beta1, gamma2, beta2, p);
			const double E2 = EnergyExpectationValue(reg, nrShots);


			Exec(reg, gamma1, beta1 - eps, gamma2, beta2, p);
			const double E3 = EnergyExpectationValue(reg, nrShots);

			Exec(reg, gamma1, beta1 + eps, gamma2, beta2, p);
			const double E4 = EnergyExpectationValue(reg, nrShots);

			double newGamma1 = gamma1 - step * (E2 - E1) / twoEps;
			double newBeta1 = beta1 - step * (E4 - E3) / twoEps;

			if (pVal == 2)
			{
				Exec(reg, gamma1, beta1, gamma2 - eps, beta2, p);
				const double E5 = EnergyExpectationValue(reg, nrShots);

				Exec(reg, gamma1, beta1, gamma2 + eps, beta2, p);
				const double E6 = EnergyExpectationValue(reg, nrShots);

				Exec(reg, gamma1, beta1, gamma2, beta2 - eps, p);
				const double E7 = EnergyExpectationValue(reg, nrShots);

				Exec(reg, gamma1, beta1, gamma2, beta2 + eps, p);
				const double E8 = EnergyExpectationValue(reg, nrShots);

				return { newGamma1, newBeta1,
					 gamma2 - step * (E6 - E5) / twoEps, beta2 - step * (E8 - E7) / twoEps };
			}

			return { newGamma1, newBeta1, newGamma1, newBeta1 };
		}

		IsingModel model;

		QC::Gates::HadamardGate<MatrixClass> h;
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::RxGate<MatrixClass> rx;
		QC::Gates::RyGate<MatrixClass> ry;
		QC::Gates::RzGate<MatrixClass> rz;

		// TODO: make configurable
		double deltaE = 0.001;
		double gamma1Start = 1;
		double beta1Start = 1;
		double gamma2Start = 1;
		double beta2Start = 1;
		unsigned int pVal = 1;
		double epsilon = 0.0002;
		double stepSize = 0.0001;
		unsigned int nrMeasurements = 500000;
		bool betterMixing = false;
	};

} // namespace Models

