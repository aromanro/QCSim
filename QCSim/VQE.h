#pragma once


#include "QuantumAlgorithm.h"
#include "QuantumGate.h"
#include "QubitRegister.h"
#include "Utils.h"
#include "PauliString.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include "PauliString.h"

// roughly based on Lesson 11 from
// "Fundamentals In Quantum Algorithms: A Tutorial Series Using Qiskit Continued" by Daniel Koch, Saahil Patel, Laura Wessing, Paul M. Alsing
// https://arxiv.org/abs/2008.10647

namespace VQE {

	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class PauliStringVQE :
		public QC::QuantumSubAlgorithm<VectorClass, MatrixClass>
	{
	public:
		PauliStringVQE()
		{
			const double theta = 0.5 * M_PI;
			rx.SetTheta(theta);
			ry.SetTheta(-theta);
		}

		unsigned int Execute(QC::QubitRegister<VectorClass, MatrixClass>& reg) override
		{
			const unsigned int nrQubits = reg.getNrQubits();
			const unsigned int lastQubit = nrQubits - 1;

			for (unsigned int qubit = 0; qubit < nrQubits; ++qubit)
				if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opX)
					reg.ApplyGate(ry, qubit);
				else if (getOperatorForQubit(qubit) == PauliString::PauliString::PauliOp::opY)
					reg.ApplyGate(rx, qubit);

			return 0; // not used here
		}

		void setOperatorForQubit(unsigned int qubit, PauliString::PauliString::PauliOp op)
		{
			pauliString.setOperatorForQubit(qubit, op);
		}

		PauliString::PauliString::PauliOp getOperatorForQubit(unsigned int qubit) const
		{
			return pauliString.getOperatorForQubit(qubit);
		}

		double getCoefficient() const
		{
			return pauliString.getCoefficient();
		}

		void setCoefficient(double c)
		{
			pauliString.setCoefficient(c);
		}

		void SingleQubitAnsatz(QC::QubitRegister<VectorClass, MatrixClass>& reg, unsigned int qubit, double theta, double phi)
		{
			rya.SetTheta(theta);
			reg.ApplyGate(rya, qubit);
			rza.SetTheta(phi);
			reg.ApplyGate(rza, qubit);
		}

		void Ansatz(QC::QubitRegister<VectorClass, MatrixClass>& reg, const std::vector<double>& params)
		{
			const unsigned int nrQubits = reg.getNrQubits();

			if (params.size() < nrQubits) return;

			unsigned int pos = 0;
			for (unsigned int q = 0; q < nrQubits; ++q)
			{
				SingleQubitAnsatz(reg, q, params[pos], params[pos + 1]);
				pos += 2;
			}

			bool flip = false;
			while (pos < params.size() - 1) {
				if (flip)
				{
					for (unsigned int q = nrQubits - 1; q > 0; --q)
						reg.ApplyGate(cnot, q - 1, q);

					if (nrQubits > 2)
						reg.ApplyGate(nrQubits - 1, 0);
				}
				else
				{
					for (unsigned int q = 0; q < nrQubits - 1; ++q)
						reg.ApplyGate(cnot, q + 1, q);

					if (nrQubits > 2)
						reg.ApplyGate(0, nrQubits - 1);
				}

				flip = !flip;

				for (unsigned int q = 0; q < nrQubits; ++q)
				{
					if (pos >= params.size() - 1) break;
					SingleQubitAnsatz(reg, q, params[pos], params[pos + 1]);
					pos += 2;
				}
			}
		}

		double EstimateEnergy(QC::QubitRegister<VectorClass, MatrixClass>& reg, const std::vector<double>& params, size_t nrMeasurements = 10000)
		{
			reg.setToBasisState(0);
			Ansatz(reg, params);
			Execute(reg);

			double energy = 0.0;

			const auto measurements = reg.RepeatedMeasure(nrMeasurements);

			for (const auto& m : measurements)
			{
				double e = 1.0;
				unsigned int val = m.first;
				for (unsigned int q = 0; val && q < reg.getNrQubits(); ++q)
				{
					if (val & 1)
						e *= -1;
					val >>= 1;
				}

				energy += e * m.second;
			}


			return getCoefficient() * energy / nrMeasurements;
		}

	protected:
		PauliString::PauliString pauliString;

		QC::Gates::RxGate<MatrixClass> rx;
		QC::Gates::RyGate<MatrixClass> ry;


		QC::Gates::RyGate<MatrixClass> rya;
		QC::Gates::RzGate<MatrixClass> rza;
		QC::Gates::CNOTGate<MatrixClass> cnot;
	};


	template<class VectorClass = Eigen::VectorXcd, class MatrixClass = Eigen::MatrixXcd> class VariationalQuantumEigensolver :
		public QC::QuantumAlgorithm<VectorClass, MatrixClass>
	{
	public:
		using BaseClass = QC::QuantumAlgorithm<VectorClass, MatrixClass>;

		VariationalQuantumEigensolver(unsigned int N = 2, bool addSeed = false)
			: BaseClass(N, addSeed), nrMeasurements(10000)
		{
		}

		void AddTerm(double coeff, const PauliStringVQE<VectorClass, MatrixClass>& term)
		{
			terms.push_back(term);
			terms.back().setCoefficient(coeff);
		}

		void Clear()
		{
			terms.clear();
		}

		unsigned int Execute() override
		{
			// TODO: implement it, using Nelder-Mead algorithm
			return 0;
		}

		static std::vector<double> Centroid(const std::vector<std::vector<double>>& points, int excludeIndex = -1)
		{
			if (points.empty()) return {};

			std::vector<double> centroid(points[0].size(), 0.0);

			for (int pi = 0; pi < points.size(); ++pi)
			{
				if (pi == excludeIndex) continue;

				const auto& p = points[pi];
				for (unsigned int i = 0; i < p.size(); ++i)
					centroid[i] += p[i];
			}

			size_t nrPoints = excludeIndex == -1 ? points.size() - 1 : points.size();
			for (auto& c : centroid)
				c /= nrPoints;

			return centroid;
		}

		static std::vector<double> ReflectionPoint(const std::vector<double>& p1, const std::vector<double>& p2, double alpha)
		{
			assert(p1.size() == p2.size());

			std::vector<double> rp(p1.size(), 0.0);

			for (unsigned int i = 0; i < p1.size(); ++i)
			{
				const double d = p2[i] - p1[i];
				rp[i] = p1[i] + alpha * d;
			}

			return rp;
		}

		double EstimateEnergy(const std::vector<double>& params, size_t nrMeasurements = 10000)
		{
			double energy = 0.0;

			for (const auto& term : terms)
				energy += term.EstimateEnergy(BaseClass::reg, params, nrMeasurements);

			return energy;
		}


		void NelderMeadStep(std::vector<std::vector<double>>& points)
		{
			if (points.empty()) return;

			double maxEnergy = EstimateEnergy(BaseClass::reg, points[0], nrMeasurements);
			int maxIndex = 0;

			double maxEnergy2 = maxEnergy;
			int maxIndex2 = maxIndex;

			double minEnergy = maxEnergy;
			int minIndex = maxIndex;

			for (int i = 1; i < points.size(); ++i)
			{
				const double e = EstimateEnergy(BaseClass::reg, points[i], nrMeasurements);
				if (e > maxEnergy)
				{
					maxEnergy2 = maxEnergy;
					maxIndex2 = maxIndex;
					maxEnergy = e;
					maxIndex = i;
				}
				else if (e < minEnergy)
				{
					minEnergy = e;
					minIndex = i;
				}
			}

			const auto centroid = Centroid(points, minIndex);
			auto newPoint = ReflectionPoint(centroid, points[minIndex], 2.0);

			const double reflectedEnergy = EstimateEnergy(BaseClass::reg, newPoint, nrMeasurements);

			// TODO: now using information from the 'reflected' point, decide what point should be chosen to replace the worst point
			// options are: reflect / expand / contract / shrink
			bool shrink = false;
			if (reflectedEnergy < minEnergy)
			{
				// we're on the right track, try to expand further
				const auto expandPoint = ReflectionPoint(centroid, newPoint, 2.0);
				const double expandEnergy = EstimateEnergy(BaseClass::reg, expandPoint, nrMeasurements);
				if (expandEnergy < reflectedEnergy)
					newPoint = expandPoint;
				points[maxIndex] = newPoint;
			} 
			else if (reflectedEnergy >= maxEnergy2)
			{
				if (reflectedEnergy < maxEnergy)
				{
					// with the reflected point we're worse than any other point (except the one that was the worst and was reflected)
					// try a point between the centroid and the reflected point, maybe it's better
					// if not, the reflected one is still better than the original one
					const auto contractPoint = ReflectionPoint(centroid, newPoint, 0.5);
					const double contractEnergy = EstimateEnergy(BaseClass::reg, contractPoint, nrMeasurements);
					if (contractEnergy < reflectedEnergy)
						points[maxIndex] = contractPoint;
					else
					{
						points[maxIndex] = newPoint;
						// also shrink
						shrink = true;
					}
				}
				else
				{
					// with the reflected point we're worse than all the points we were starting with
					// including the one that was reflected
					// try a point between the origianl point and the centroid, maybe it's better
					const auto contractPoint = ReflectionPoint(centroid, points[maxIndex], 0.5);
					const double contractEnergy = EstimateEnergy(BaseClass::reg, contractPoint, nrMeasurements);
					if (contractEnergy < minEnergy)
						points[maxIndex] = contractPoint;
					else
						// shrink
						shrink = true;
				}
			}
			else points[maxIndex] = newPoint;

			if (shrink)
			{
				// shrink all points (except the min one)
				// by reflecting them about the min point, with alpha = 0.5
				for (int i = 0; i < points.size(); ++i)
				{
					if (i == minIndex) continue;
					points[i] = ReflectionPoint(points[minIndex], points[i], 0.5);
				}
			}
		}

		void SetNrMeasurements(size_t n)
		{
			nrMeasurements = n;
		}

		size_t GetNrMeasurements() const
		{
			return nrMeasurements;
		}

	protected:
		std::vector<PauliStringVQE<VectorClass, MatrixClass>> terms;
		size_t nrMeasurements;
	};

}