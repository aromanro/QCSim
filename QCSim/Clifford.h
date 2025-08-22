#pragma once

#include "StabilizerState.h"

namespace QC {
	namespace Clifford {

		class StabilizerSimulator : public StabilizerState {
		public:
			StabilizerSimulator() = delete;

			explicit StabilizerSimulator(size_t nQubits)
				: StabilizerState(nQubits)
			{
			}

			void ApplyH(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyH(qubit, q);
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyH(qubit, q);
				}
			}

			void ApplyK(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplyS(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyS(qubit, q);
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyS(qubit, q);
				}
			}

			void ApplySdg(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyZ(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplySx(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 512)
				{
					for (size_t q = 0; q < nrQubits; ++q)
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
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 128)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
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
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyS(qubit, q);
						ApplyH(qubit, q);
						ApplyS(qubit, q);
					}
				}
			}

			void ApplyX(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyX(qubit, q);
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyX(qubit, q);
				}
			}

			void ApplyY(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyY(qubit, q);
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyY(qubit, q);
				}
			}

			void ApplyZ(size_t qubit)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 2048)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyZ(qubit, q);
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 512)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyZ(qubit, q);
				}
			}

			void ApplyCX(size_t target, size_t control)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
						ApplyCX(target, control, q);
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
						ApplyCX(target, control, q);
				}
			}

			void ApplyCY(size_t target, size_t control)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyZ(target, q);
						ApplyS(target, q);
						ApplyCX(target, control, q);
						ApplyS(target, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyZ(target, q);
						ApplyS(target, q);
						ApplyCX(target, control, q);
						ApplyS(target, q);
					}
				}
			}

			void ApplyCZ(size_t target, size_t control)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyH(target, q);
						ApplyCX(target, control, q);
						ApplyH(target, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyH(target, q);
						ApplyCX(target, control, q);
						ApplyH(target, q);
					}
				}
			}

			void ApplySwap(size_t qubit1, size_t qubit2)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 1024)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 256)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
					}
				}
			}

			void ApplyISwap(size_t qubit1, size_t qubit2)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 512)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyS(qubit1, q);
						ApplyS(qubit2, q);
						ApplyH(qubit1, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
						ApplyH(qubit2, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 128)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyS(qubit1, q);
						ApplyS(qubit2, q);
						ApplyH(qubit1, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyCX(qubit1, qubit2, q);
						ApplyH(qubit2, q);
					}
				}
			}

			std::unique_ptr<StabilizerSimulator> Clone() const
			{
				auto sim = std::make_unique<StabilizerSimulator>(1);

				sim->destabilizerGenerators = destabilizerGenerators;
				sim->stabilizerGenerators = stabilizerGenerators;

				sim->savedDestabilizerGenerators = savedDestabilizerGenerators;
				sim->savedStabilizerGenerators = savedStabilizerGenerators;

				sim->enableMultithreading = enableMultithreading;

				return sim;
			}

		private:
			inline void ApplyH(size_t qubit, size_t q)
			{
				// how does this work?
				// looking at it might not reveal immediately what it does
				// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

				// if we have both X and Z, then it's an Y (with some global phase, given by the sign, Y = iXZ)
				// a Y is transformed to a -Y, so a sign change is needed

				if (destabilizerGenerators[q].X[qubit] && destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				// swap X and Z
				bool t = destabilizerGenerators[q].X[qubit];
				destabilizerGenerators[q].X[qubit] = destabilizerGenerators[q].Z[qubit];
				destabilizerGenerators[q].Z[qubit] = t;

				if (stabilizerGenerators[q].X[qubit] && stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;

				// swap X and Z
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
		};
	}
}
