#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cassert>

#include "QubitRegisterCalculator.h"

namespace QC {
	namespace Clifford {

		class Generator {
		public:
			Generator() {}

			Generator(size_t nQubits) : X(nQubits, false), Z(nQubits, false) { }

			Generator(const Generator& other) : X(other.X), Z(other.Z), PhaseSign(other.PhaseSign) { }

			Generator(Generator&& other) noexcept : X(std::move(other.X)), Z(std::move(other.Z)), PhaseSign(other.PhaseSign) { }

			Generator& operator=(const Generator& other)
			{
				if (this != &other)
				{
					X = other.X;
					Z = other.Z;
					PhaseSign = other.PhaseSign;
				}

				return *this;
			}

			Generator& operator=(Generator&& other) noexcept
			{
				if (this != &other)
				{
					X.swap(other.X);
					Z.swap(other.Z);
					PhaseSign = other.PhaseSign;
				}

				return *this;
			}

			void Resize(size_t nQubits) 
			{
				X.resize(nQubits, false);
				Z.resize(nQubits, false);
			}

			void Clear()
			{
				std::fill(X.begin(), X.end(), false);
				std::fill(Z.begin(), Z.end(), false);
				PhaseSign = false;
			}

			std::vector<bool> X;
			std::vector<bool> Z;
			bool PhaseSign = false;
		};

		class StabilizerSimulator {
		public:
			StabilizerSimulator() = delete;

			StabilizerSimulator(size_t nQubits)
				: destabilizerGenerators(nQubits), stabilizerGenerators(nQubits), gen(std::random_device{}()), rnd(0.5)
			{
				// this puts it in the |0> state
				for (size_t q = 0; q < nQubits; ++q) 
				{
					destabilizerGenerators[q].Resize(nQubits);
					stabilizerGenerators[q].Resize(nQubits);

					destabilizerGenerators[q].X[q] = true;
					stabilizerGenerators[q].Z[q] = true;
				}
			}

			void ApplyH(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyH(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyH(qubit, q);
				}
			}

			void ApplyS(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyS(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyS(qubit, q);
				}
			}

			void ApplySdg(size_t qubit)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplySx(size_t qubit)
			{
				if (stabilizerGenerators.size() < 512)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 128)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplySxDag(size_t qubit)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplyX(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyX(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyX(qubit, q);
				}
			}

			void ApplyY(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyY(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyY(qubit, q);
				}
			}

			void ApplyZ(size_t qubit)
			{
				if (stabilizerGenerators.size() < 2048)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyZ(qubit, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyZ(qubit, q);
				}
			}

			void ApplyCX(size_t target, size_t control)
			{
				if (stabilizerGenerators.size() < 1024)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
						ApplyCX(target, control, q);
				}
				else
				{
					const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(stabilizerGenerators.size()); ++q)
						ApplyCX(target, control, q);
				}
			}

			void ApplyCY(size_t target, size_t control)
			{
				ApplySdg(target);
				ApplyCX(target, control);
				ApplyS(target);
			}

			void ApplyCZ(size_t target, size_t control)
			{
				ApplyH(target);
				ApplyCX(target, control);
				ApplyH(target);
			}

			void ApplySwap(size_t qubit1, size_t qubit2)
			{
				ApplyCX(qubit1, qubit2);
				ApplyCX(qubit2, qubit1);
				ApplyCX(qubit1, qubit2);
			}

			bool MeasureQubit(size_t qubit)
			{
				// TODO: parellelize this
				bool case1 = false;
				size_t p = 0;
				for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
				{
					// Z anticommutes with X
					if (stabilizerGenerators[q].X[qubit])
					{
						p = q;
						case1 = true;
						break;
					}
				}

				if (case1)
				{
					for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
					{
						if (destabilizerGenerators[q].X[qubit]) 
							rowsum(destabilizerGenerators[q], p);

						if (stabilizerGenerators[q].X[qubit] && p != q) 
							rowsum(stabilizerGenerators[q], destabilizerGenerators.size() + p);
					}

					destabilizerGenerators[p] = stabilizerGenerators[p];
					
					stabilizerGenerators[p].Clear();
					stabilizerGenerators[p].Z[qubit] = true;
					stabilizerGenerators[p].PhaseSign = rnd(gen);

					return stabilizerGenerators[p].PhaseSign;
				}
				
				// case 2 - Z (on measured qubit) commutes with all generators
				// no change to generators, just need to compute the sign in order to get the measurement result
				Generator h(stabilizerGenerators.size());

				for (size_t q = 0; q < destabilizerGenerators.size(); ++q)
				{
					if (destabilizerGenerators[q].X[qubit])
						rowsum(h, destabilizerGenerators.size() + q);
				}

				return h.PhaseSign;
			}

		private:
			inline static bool XOR(bool a, bool b) 
			{ 
				//return ((a ? 1 : 0) ^ (b ? 1 : 0)) == 1;
				return a != b;
			}
			
			inline void ApplyH(size_t qubit, size_t q)
			{
				// how does this work?
				// looking at it might not reveal immediately what it does
				// the swaps below just switches the X with Z, because HZH^t = X and HXH^t = Z

				// if we have both X and Z, then it's an Y (with some global phase, given by the sign, Y = iXZ)
				// an Y is transformed to a -Y, so a singn change is needed
				
				if (destabilizerGenerators[q].X[qubit] && destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;
				
				bool t = destabilizerGenerators[q].X[qubit];
				destabilizerGenerators[q].X[qubit] = destabilizerGenerators[q].Z[qubit];
				destabilizerGenerators[q].Z[qubit] = t;

				if (stabilizerGenerators[q].X[qubit] && stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;

				t = stabilizerGenerators[q].X[qubit];
				stabilizerGenerators[q].X[qubit] = stabilizerGenerators[q].Z[qubit];
				stabilizerGenerators[q].Z[qubit] = t;
			}

			inline void ApplyS(size_t qubit, size_t q)
			{
				if (destabilizerGenerators[q].X[qubit] && destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				destabilizerGenerators[q].Z[qubit] = XOR(destabilizerGenerators[q].Z[qubit], destabilizerGenerators[q].X[qubit]);

				if (stabilizerGenerators[q].X[qubit] && stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
				
				stabilizerGenerators[q].Z[qubit] = XOR(stabilizerGenerators[q].Z[qubit], stabilizerGenerators[q].X[qubit]);
			}

			inline void ApplyX(size_t qubit, size_t q)
			{
				if (destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				if (stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
			}

			inline void ApplyY(size_t qubit, size_t q)
			{
				// can be done with ifs, can be done with XORs
				destabilizerGenerators[q].PhaseSign = XOR(destabilizerGenerators[q].PhaseSign, XOR(destabilizerGenerators[q].Z[qubit], destabilizerGenerators[q].X[qubit]));
				stabilizerGenerators[q].PhaseSign = XOR(stabilizerGenerators[q].PhaseSign, XOR(stabilizerGenerators[q].Z[qubit], stabilizerGenerators[q].X[qubit]));
			}

			inline void ApplyZ(size_t qubit, size_t q)
			{
				if (destabilizerGenerators[q].X[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				if (stabilizerGenerators[q].X[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
			}

			inline void ApplyCX(size_t target, size_t control, size_t q)
			{
				destabilizerGenerators[q].PhaseSign = XOR(destabilizerGenerators[q].PhaseSign, destabilizerGenerators[q].X[control] && destabilizerGenerators[q].Z[target] &&
					XOR(destabilizerGenerators[q].X[target], XOR(destabilizerGenerators[q].Z[control], true)));

				destabilizerGenerators[q].X[target] = XOR(destabilizerGenerators[q].X[target], destabilizerGenerators[q].X[control]);
				destabilizerGenerators[q].Z[control] = XOR(destabilizerGenerators[q].Z[control], destabilizerGenerators[q].Z[target]);

				stabilizerGenerators[q].PhaseSign = XOR(stabilizerGenerators[q].PhaseSign, stabilizerGenerators[q].X[control] && stabilizerGenerators[q].Z[target] &&
					XOR(stabilizerGenerators[q].X[target], XOR(stabilizerGenerators[q].Z[control], true)));

				stabilizerGenerators[q].X[target] = XOR(stabilizerGenerators[q].X[target], stabilizerGenerators[q].X[control]);
				stabilizerGenerators[q].Z[control] = XOR(stabilizerGenerators[q].Z[control], stabilizerGenerators[q].Z[target]);
			}

			static inline int g(int x1, int z1, int x2, int z2)
			{
				if (0 == x1 && 0 == z1) return 0;
				else if (1 == x1)
				{
					if (1 == z1) return z2 - x2;
					
					return z2 * (2 * x2 - 1);
				}

				return x2 * (1 - 2 * z2);
			}

			inline void rowsumDestabilizers(Generator& h, size_t j, long long int& m)
			{
				m += destabilizerGenerators[j].PhaseSign ? 2 : 0;

				for (size_t q = 0; q < destabilizerGenerators.size(); ++q)
				{
					const int x = h.X[q] ? 1 : 0;
					const int z = h.Z[q] ? 1 : 0;
					m += g(destabilizerGenerators[j].X[q] ? 1 : 0, destabilizerGenerators[j].Z[q] ? 1 : 0, x, z);

					h.X[q] = XOR(h.X[q], destabilizerGenerators[j].X[q]);
					h.Z[q] = XOR(h.Z[q], destabilizerGenerators[j].Z[q]);
				}
			}

			inline void rowsumStabilizers(Generator& h, size_t j, long long int& m)
			{
				m += stabilizerGenerators[j].PhaseSign ? 2 : 0;

				for (size_t q = 0; q < stabilizerGenerators.size(); ++q)
				{
					const int x = h.X[q] ? 1 : 0;
					const int z = h.Z[q] ? 1 : 0;
					m += g(stabilizerGenerators[j].X[q] ? 1 : 0, stabilizerGenerators[j].Z[q] ? 1 : 0, x, z);

					h.X[q] = XOR(h.X[q], stabilizerGenerators[j].X[q]);
					h.Z[q] = XOR(h.Z[q], stabilizerGenerators[j].Z[q]);
				}
			}

			inline void rowsum(Generator& h, size_t j)
			{
				long long int m = h.PhaseSign ? 2 : 0;

				if (j >= destabilizerGenerators.size())
					rowsumStabilizers(h, j - destabilizerGenerators.size(), m);
				else
					rowsumDestabilizers(h, j, m);

				m %= 4;
				assert(m == 0 || m == 2);

				h.PhaseSign = m != 0;
			}

			std::vector<Generator> destabilizerGenerators;
			std::vector<Generator> stabilizerGenerators;

			std::default_random_engine gen;
			std::bernoulli_distribution rnd;
		};
	}
}
