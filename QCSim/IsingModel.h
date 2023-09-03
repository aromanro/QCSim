#pragma once

#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>

#include "QuantumGate.h"
#include "QuantumAlgorithm.h"

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
			for (size_t pos = 0; state && pos < spins.size(); ++pos)
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


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class IsingSubalgorithm : public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using RegisterClass = QC::QubitRegister<VectorClass, MatrixClass>;

		unsigned int Execute(RegisterClass& reg) override
		{
			// TODO: implement it

			return 0;
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

				rz.SetTheta(gamma * hval);
				reg.ApplyGate(rz, q);
			}
		}

		void ApplyMixingOperator(RegisterClass& reg, double beta, int nrTimes = 1)
		{
			const double twobeta = 2 * beta;
			rx.SetTheta(twobeta);
			ry.SetTheta(twobeta);

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				for (size_t q = 0; q < reg.getNrQubits(); ++q)
					reg.ApplyGate(rx, q);


				unsigned int lastQubit = reg.getNrQubits() - 1;
				for (size_t q = 0; q < lastQubit; ++q)
					reg.ApplyGate(cnot, q + 1, q);

				reg.ApplyGate(cnot, 0, lastQubit);

				for (size_t q = 0; q < reg.getNrQubits(); ++q)
					reg.ApplyGate(ry, q);
			}
		}


		void Exec(RegisterClass& reg, double gamma, double beta, unsigned int nrTimes = 1)
		{
			reg.setToBasisState(0);
			ApplyHadamardOnAll(reg);
			ApplyIsingOperator(reg, gamma);
			ApplyMixingOperator(reg, beta, nrTimes);
		}

		double EnergyExpectationValue(RegisterClass& reg, unsigned int nrShots = 1000)
		{
			double energy = 0;

			auto res = reg.RepeatedMeasure(nrShots);
			for (const auto& r : res)
				energy += Energy(r.first) * static_cast<double>(r.second) / nrShots;

			return energy;
		}


		std::pair<double, double> GradientDescentStep(RegisterClass& reg, double gamma, double beta, unsigned int nrTimes = 2, double eps = 0.001, double step = 0.01, unsigned int nrShots = 10000)
		{
			const double twoEps = 2. * eps;

			Exec(reg, gamma - eps, beta, nrTimes);
			const double E1 = EnergyExpectationValue(reg, nrShots);

			Exec(reg, gamma + eps, beta, nrTimes);
			const double E2 = EnergyExpectationValue(reg, nrShots);


			Exec(reg, gamma, beta - eps, nrTimes);
			const double E3 = EnergyExpectationValue(reg, nrShots);

			Exec(reg, gamma, beta + eps, nrTimes);
			const double E4 = EnergyExpectationValue(reg, nrShots);

			return { gamma - step * (E2 - E1) / twoEps, beta - step * (E4 - E3) / twoEps };
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

	protected:
		void ApplyHadamardOnAll(RegisterClass& reg)
		{
			for (size_t q = 0; q < reg.getNrQubits(); ++q)
				reg.ApplyGate(h, q);
		}

		IsingModel model;

		QC::Gates::HadamardGate<MatrixClass> h;
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::RxGate<MatrixClass> rx;
		QC::Gates::RyGate<MatrixClass> ry;
		QC::Gates::RzGate<MatrixClass> rz;
	};

} // namespace Models

