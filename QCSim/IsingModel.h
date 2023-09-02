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

			interactions[std::make_pair(i, j)] = J;
			neighbours[i].insert(j);
		}

		void Set(const std::vector<bool> s, const std::vector<Interaction>& ints = {}, const std::vector<double>& hvals = {})
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

		double SetState(unsigned int state)
		{
			size_t pos = 0;
			while (state && pos < spins.size())
			{
				spins[pos] = state & 1 ? true : false;
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
			return s ? 1 : -1;
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

					reg.Apply(cnot, q1, q2);
					reg.Apply(rz, q2);
					reg.Apply(cnot, q1, q2);
				}

			// on site
			for (size_t q = 0; q < reg.getNrQubits(); ++q)
			{
				const double hval = model.GetH(q);

				rz.SetTheta(2 * gamma * hval);
				reg.Apply(rz, q);
			}
		}

		void ApplyMixingOperator(RegisterClass& reg, double beta, int nrTimes = 1)
		{
			const double twobeta = 2 * beta;
			rx.SetTheta(twobeta);
			ry.SetTheta(twobeta);

			for (size_t q = 0; q < reg.getNrQubits(); ++q)
				reg.Apply(rx, q);

			for (unsigned int i = 0; i < nrTimes; ++i)
			{
				unsigned int lastQubit = reg.getNrQubits() - 1;
				for (size_t q = 0; q < lastQubit; ++q)
					reg.Apply(cnot, q + 1, q);

				reg.Apply(cnot, 0, lastQubit);
			}

			for (size_t q = 0; q < reg.getNrQubits(); ++q)
				reg.Apply(ry, q);
		}

	protected:
		void ApplyHadamardOnAll(RegisterClass& reg)
		{
			for (size_t q = 0; q < reg.getNrQubits(); ++q)
				reg.Apply(h, q);
		}

		IsingModel model;

		QC::Gates::HadamardGate<MatrixClass> h;
		QC::Gates::CNOTGate<MatrixClass> cnot;
		QC::Gates::RxGate<MatrixClass> rx;
		QC::Gates::RyGate<MatrixClass> ry;
		QC::Gates::RzGate<MatrixClass> rz;
	};

} // namespace Models

