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

		void SetState(size_t state)
		{
			for (size_t pos = 0; pos < spins.size(); ++pos)
			{
				spins[pos] = (state & 1) ? true : false;
				state >>= 1;
			}
		}

		double Energy(size_t state)
		{
			SetState(state);
			return Energy();
		}

		size_t GetMinEnergyState()
		{
			size_t maxStates = 1 << spins.size();

			double minEnergy = std::numeric_limits<double>::max();
			size_t minState = 0;

			for (size_t state = 0; state < maxStates; ++state)
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

		std::set<size_t> GetMinEnergyStates()
		{
			std::set<size_t> res;
			size_t maxStates = 1 << spins.size();

			double minEnergy = std::numeric_limits<double>::max();
			size_t minState = 0;

			for (size_t state = 0; state < maxStates; ++state)
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

		size_t Execute(RegisterClass& reg) override
		{
			double Eold = std::numeric_limits<double>::max();
			double Emin = Eold;
			double E = 0; // whatever
			
			std::pair<std::vector<double>, std::vector<double>> gammaBeta = { gammaStart, betaStart };
			std::pair<std::vector<double>, std::vector<double>> optGammaBeta = gammaBeta;

			
			sGamma.resize(gammaStart.size(), 0.);
			mGamma.resize(gammaStart.size(), 0.);


			sBeta.resize(betaStart.size(), 0.);
			mBeta.resize(betaStart.size(), 0.);
			stepNr = 0;

			for (int i = 0; abs(E - Eold) > deltaE && i < 100000; ++i)
			{
				Eold = E;

				// TODO: Maybe I should change to Nelder-Mead - see the VQE implementation for an example
				gammaBeta = GradientDescentStep(reg, gammaBeta.first, gammaBeta.second, epsilon, stepSize, nrMeasurements);

				Exec(reg, gammaBeta.first, gammaBeta.second);
				E = EnergyExpectationValue(reg);

				if (E < Emin)
				{
					Emin = E;
					optGammaBeta = gammaBeta;
				}

				//if (i % 10 == 0)
					std::cout << "E: " << E << std::endl;
			}

			Exec(reg, optGammaBeta.first, optGammaBeta.second);
			auto res = reg.RepeatedMeasure(nrMeasurements);

			Emin = std::numeric_limits<double>::max();
			size_t state = 0;
			for (const auto& r : res)
			{
				E = Energy(r.first) * static_cast<double>(r.second) / nrMeasurements;
				if (E < Emin)
				{
					Emin = E;
					state = r.first;
				}
			}

			sGamma.clear();
			mGamma.clear();

			sBeta.clear();
			mBeta.clear();

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

		void SetState(size_t state)
		{
			model.SetState(state);
		}

		double Energy(size_t state)
		{
			return model.Energy(state);
		}

		size_t GetMinEnergyState()
		{
			return model.GetMinEnergyState();
		}

		std::set<size_t> GetMinEnergyStates()
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

		const std::vector<double>& GetGammaStart() const
		{
			return gammaStart;
		}

		void SetGammaStart(double gamma, size_t index = 0)
		{
			if (index >= gammaStart.size()) return;

			gammaStart[index] = gamma;
		}

		const std::vector<double>& GetBetaStart() const
		{
			return betaStart;
		}

		void SetBetaStart(double beta, size_t index = 0)
		{
			if (index >= betaStart.size()) return;

			betaStart[index] = beta;
		}

		size_t GetP() const
		{
			return gammaStart.size();
		}

		void SetP(size_t p)
		{
			if (p < 1) return;

			gammaStart.resize(p, 1.);
			betaStart.resize(p, 1.);
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

		void SetBeta1(double b1)
		{
			beta1 = b1;
		}

		double GetBeta1() const
		{
			return beta1;
		}

		void SetBeta2(double b2)
		{
			beta2 = b2;
		}

		double GetBeta2() const
		{
			return beta2;
		}

		void SetLambda(double l)
		{
			lambda = l;
		}

		double GetLambda() const
		{
			return lambda;
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
					const size_t q1 = site.first;
					const size_t q2 = n;
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

				size_t lastQubit = reg.getNrQubits() - 1;
				for (size_t q = 0; q < lastQubit; ++q)
					reg.ApplyGate(cnot, q + 1, q);
				reg.ApplyGate(cnot, 0, lastQubit);

				for (size_t q = 0; q < reg.getNrQubits(); ++q)
					reg.ApplyGate(ry, q);
			}
		}


		void Exec(RegisterClass& reg, const std::vector<double>& gamma, const std::vector<double>& beta)
		{
			reg.setToBasisState(0);
			ApplyHadamardOnAll(reg);

			for (size_t i = 0; i < gamma.size(); ++i)
			{
				ApplyIsingOperator(reg, gamma[i]);
				ApplyMixingOperator(reg, beta[i], betterMixing);
			}
		}

		double EnergyExpectationValue(RegisterClass& reg, size_t nrShots = 100000)
		{
			double energy = 0;

			auto res = reg.RepeatedMeasure(nrShots);
			for (const auto& r : res)
				energy += Energy(r.first) * static_cast<double>(r.second) / nrShots;

			return energy;
		}

		std::pair<std::vector<double>, std::vector<double>> GradientDescentStep(RegisterClass& reg, std::vector<double>& gamma, std::vector<double>& beta, double eps = 0.0002, double alpha = 0.0001, size_t nrShots = 100000)
		{
			const double twoEps = 2. * eps;

			std::pair<std::vector<double>, std::vector<double>> res;
			res.first.resize(gamma.size());
			res.second.resize(beta.size());

			++stepNr;
			const double div1 = 1. / (1. - pow(beta1, stepNr));
			const double div2 = 1. / (1. - pow(beta2, stepNr));

			for (size_t i = 0; i < gamma.size(); ++i)
			{
				const double gammaVal = gamma[i];
				const double betaVal = beta[i];

				gamma[i] = gammaVal - eps;
				Exec(reg, gamma, beta);
				const double E1 = EnergyExpectationValue(reg, nrShots);
				gamma[i] = gammaVal + eps;
				Exec(reg, gamma, beta);
				const double E2 = EnergyExpectationValue(reg, nrShots);
				gamma[i] = gammaVal;

				beta[i] = betaVal - eps;
				Exec(reg, gamma, beta);
				const double E3 = EnergyExpectationValue(reg, nrShots);
				beta[i] = betaVal + eps;
				Exec(reg, gamma, beta);
				const double E4 = EnergyExpectationValue(reg, nrShots);
				beta[i] = betaVal;

				// ADAMW
				double grad = (E2 - E1) / twoEps;

				mGamma[i] = beta1 * mGamma[i] - (1. - beta1) * grad;
				mGamma[i] *= div1;
				sGamma[i] = beta2 * sGamma[i] + (1. - beta2) * grad * grad;
				sGamma[i] *= div2;

				res.first[i] = gammaVal + alpha * (mGamma[i] / sqrt(sGamma[i] + 1E-8) - lambda * gammaVal);

				grad = (E4 - E3) / twoEps;

				mBeta[i] = beta1 * mBeta[i] - (1. - beta1) * grad;
				mBeta[i] *= div1;
				sBeta[i] = beta2 * sBeta[i] + (1. - beta2) * grad * grad;
				sBeta[i] *= div2;

				res.second[i] = betaVal + alpha * (mBeta[i] / sqrt(sBeta[i] + 1E-8) - lambda * betaVal);
			}

			return res;
		}

		IsingModel model;

		QC::Gates::HadamardGate<MatrixClass> h;
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::RxGate<MatrixClass> rx;
		QC::Gates::RyGate<MatrixClass> ry;
		QC::Gates::RzGate<MatrixClass> rz;

		double deltaE = 0.0005;

		std::vector<double> gammaStart = { 1. };
		std::vector<double> betaStart = { 1. };

		std::vector<double> sGamma;
		std::vector<double> mGamma;

		std::vector<double> sBeta;
		std::vector<double> mBeta;

		double beta1 = 0.8;
		double beta2 = 0.95;
		double lambda = 0.001; // with lambda 0 is the same as Adam
		int stepNr = 0;


		double epsilon = 0.001;
		double stepSize = 0.0005;
		size_t nrMeasurements = 100000;
		bool betterMixing = false;
	};

} // namespace Models

