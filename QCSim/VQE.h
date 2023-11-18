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
		PauliStringVQE(int nrQubits)
			: QC::QuantumSubAlgorithm<VectorClass, MatrixClass>(), pauliString(nrQubits)
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
			for (unsigned int q = 0; q < nrQubits && pos < params.size() - 1; ++q)
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
						reg.ApplyGate(cnot, nrQubits - 1, 0);
				}
				else
				{
					for (unsigned int q = 0; q < nrQubits - 1; ++q)
						reg.ApplyGate(cnot, q + 1, q);

					if (nrQubits > 2)
						reg.ApplyGate(cnot, 0, nrQubits - 1);
				}

				flip = !flip;

				for (unsigned int q = 0; q < nrQubits && pos < params.size() - 1; ++q)
				{
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
			: BaseClass(N, addSeed)
		{
		}

		void AddTerm(const PauliStringVQE<VectorClass, MatrixClass>& term)
		{
			terms.emplace_back(term);
		}

		void Clear()
		{
			terms.clear();
		}

		unsigned int Execute() override
		{
			if (vertices.empty()) return 0;

			EstimateEnergiesForVertices();
			double oldLowestEnergy = GetMinEnergy();

			int terminateCount = 0;
			for (int i = 0; i < 1000; ++i)
			{
				NelderMeadStep();

				double newLowestEnergy = GetMinEnergy();

				if (abs(newLowestEnergy - oldLowestEnergy) < delta)
					++terminateCount;
				else
					terminateCount = 0;

				oldLowestEnergy = newLowestEnergy;

				if (i % 5 == 0)
					std::cout << "Iteration " << i << ", Energy: " << newLowestEnergy << std::endl;

				if (terminateCount >= terminateLimit) break;
			}

			return 0;
		}

		double EstimateEnergy(const std::vector<double>& params, size_t nrMeasurements = 10000)
		{
			double energy = 0.0;

			for (auto& term : terms)
				energy += term.EstimateEnergy(BaseClass::reg, params, nrMeasurements);

			return energy;
		}

		void SetNrMeasurements(size_t n)
		{
			nrMeasurements = n;
		}

		size_t GetNrMeasurements() const
		{
			return nrMeasurements;
		}

		void SetVertices(const std::vector<std::vector<double>>& v)
		{
			vertices = v;
		}

		const std::vector<std::vector<double>>& GetVertices() const
		{
			return vertices;
		}

		void SetTerminateLimit(int limit)
		{
			terminateLimit = limit;
		}

		int GetTerminateLimit() const
		{
			return terminateLimit;
		}

		void SetDelta(double d)
		{
			delta = d;
		}

		double GetDelta() const
		{
			return delta;
		}

		double GetMinEnergy() const
		{
			if (vertexEnergies.empty()) return 0.0;

			double minEnergy = vertexEnergies[0];
			for (int i = 1; i < vertexEnergies.size(); ++i)
				if (vertexEnergies[i] < minEnergy)
					minEnergy = vertexEnergies[i];

			return minEnergy;
		}

		std::vector<double> GetMinVertex() const
		{
			if (vertexEnergies.empty()) return {};

			double minEnergy = vertexEnergies[0];
			int minIndex = 0;
			for (int i = 1; i < vertexEnergies.size(); ++i)
				if (vertexEnergies[i] < minEnergy)
				{
					minEnergy = vertexEnergies[i];
					minIndex = i;
				}

			return vertices[minIndex];
		}

	protected:
		void NelderMeadStep()
		{
			if (vertices.empty()) return;

			double maxEnergy = 0;
			int maxIndex = 0;

			double maxEnergy2 = maxEnergy;
			int maxIndex2 = maxIndex;

			double minEnergy = maxEnergy;
			int minIndex = maxIndex;
			
			// pick up the worst point

			// also identify the best one, need its value for comparisons and shrinking if needed
			// also keep the second worst, for comparisons

			GetMinMaxEnergy(minEnergy, maxEnergy, maxEnergy2, minIndex, maxIndex, maxIndex2);

			const auto centroid = Centroid(vertices, maxIndex);
			const auto reflectedPoint = ReflectionPoint(vertices[maxIndex], centroid, 2.0);
			const double reflectedEnergy = EstimateEnergy(reflectedPoint, nrMeasurements);

			// TODO: now using information from the 'reflected' point, decide what point should be chosen to replace the worst point
			// options are: reflect / expand / contract / shrink
			
			if (reflectedEnergy < minEnergy)
			{
				// we're on the right track, try to expand further
				const auto expandPoint = ReflectionPoint(centroid, reflectedPoint, 2.0);
				const double expandEnergy = EstimateEnergy(expandPoint, nrMeasurements);
				if (expandEnergy < reflectedEnergy)
				{
					vertices[maxIndex] = expandPoint;
					vertexEnergies[maxIndex] = expandEnergy;
				}
				else
				{
					vertices[maxIndex] = reflectedPoint;
					vertexEnergies[maxIndex] = reflectedEnergy;
				}
			}
			else if (reflectedEnergy > maxEnergy2)
			{
				bool shrink = false;

				if (reflectedEnergy < maxEnergy)
				{
					// with the reflected point we're worse than any other point (except the one that was the worst and was reflected)
					// try a point between the centroid and the reflected point, maybe it's better
					// if not, the reflected one is still better than the original one
					const auto contractPoint = ReflectionPoint(centroid, reflectedPoint, 0.5);
					const double contractEnergy = EstimateEnergy(contractPoint, nrMeasurements);
					if (contractEnergy < reflectedEnergy)
					{
						vertices[maxIndex] = contractPoint;
						vertexEnergies[maxIndex] = contractEnergy;
					}
					else
					{
						vertices[maxIndex] = reflectedPoint;
						vertexEnergies[maxIndex] = reflectedEnergy;
						shrink = true;
					}
				}
				else
				{
					// with the reflected point we're worse than all the points we were starting with
					// including the one that was reflected
					// try a point between the original point and the centroid, maybe it's better
					auto contractPoint = ReflectionPoint(centroid, vertices[maxIndex], 0.5);
					double contractEnergy = EstimateEnergy(contractPoint, nrMeasurements);
					if (contractEnergy < maxEnergy)
					{
						vertices[maxIndex] = contractPoint;
						vertexEnergies[maxIndex] = contractEnergy;
					}
					else
					{
						contractPoint = ReflectionPoint(centroid, vertices[maxIndex], 1.5);
						contractEnergy = EstimateEnergy(contractPoint, nrMeasurements);
						if (contractEnergy < maxEnergy)
						{
							vertices[maxIndex] = contractPoint;
							vertexEnergies[maxIndex] = contractEnergy;
						}

						// shrink
						shrink = true;
					}
				}


				if (shrink)
				{
					// shrink all points (except the min one)
					// by reflecting them about the min point, with alpha = 0.5
					for (int i = 0; i < vertices.size(); ++i)
					{
						if (i == minIndex) continue;
						vertices[i] = ReflectionPoint(vertices[minIndex], vertices[i], 0.5);
						vertexEnergies[i] = EstimateEnergy(vertices[i], nrMeasurements);
					}
				}
			}
			else 
			{
				vertices[maxIndex] = reflectedPoint;
				vertexEnergies[maxIndex] = reflectedEnergy;
			}
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

			size_t nrPoints = (excludeIndex == -1) ? points.size() : points.size() - 1;
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

		void EstimateEnergiesForVertices()
		{
			if (vertices.empty()) return;

			vertexEnergies.resize(vertices.size());
			for (int v = 0; v < vertices.size(); ++v)
				vertexEnergies[v] = EstimateEnergy(vertices[v], nrMeasurements);
		}


		void GetMinMaxEnergy(double& minEnergy, double& maxEnergy, double& maxEnergy2, int& minIndex, int& maxIndex, int& maxIndex2) const
		{
			if (vertexEnergies.empty()) return;
			maxEnergy = vertexEnergies[0];
			maxIndex = 0;

			maxEnergy2 = maxEnergy;
			maxIndex2 = 0;

			minEnergy = maxEnergy;
			minIndex = 0;

			// pick up the worst point

			// also identify the best one, need its value for comparisons and shrinking if needed
			// also keep the second worst, for comparisons

			for (int i = 1; i < vertexEnergies.size(); ++i)
			{
				const double e = vertexEnergies[i];
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
		}



		std::vector<PauliStringVQE<VectorClass, MatrixClass>> terms;
		size_t nrMeasurements = 10000;

		std::vector<std::vector<double>> vertices;
		int terminateLimit = 10;
		double delta = 0.0001;
		std::vector<double> vertexEnergies;
	};

}