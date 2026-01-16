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

			double ExpectationValue(const std::string& pauliString) const
			{
				// We compute this: <Psi|pauliString|Psi>, where |Psi> is defined by the stabilizers
				if (pauliString.empty()) return 1.0;

				Generator g(getNrQubits());
				std::vector<size_t> pos;
				pos.reserve(pauliString.size());
				size_t phase = 0;
				
				SetPauliString(pauliString, g, pos, phase);
				if (pos.empty()) return 1.0;

				// this is easy, if the pauli string transforming the stabilizers leads to an orthogonal state on the original one, the expectation value is zero
				if (CheckStabilizersAnticommutation(g, pos))
					return 0.0;

				// now we're left with the case when all stabilizers commute with the operator
				// we need to find if the result is 1 or -1

				// check destabilizers for that
				std::vector<bool> GZ(g.Z);

				for (size_t i = 0; i < getNrQubits(); ++i)
				{
					// check if the destabilizer anticommutes with the operator
					if (CommutesWithDestabilizer(g, destabilizerGenerators[i], pos))
						continue; // commutes, check next destabilizer

					// anticommutes with this destabilizer

					phase += stabilizerGenerators[i].PhaseSign ? 2 : 0;
					for (size_t q = 0; q < getNrQubits(); ++q)
					{
						if (stabilizerGenerators[i].X[q] && GZ[q])
							phase += 2;
						if (stabilizerGenerators[i].X[q] && stabilizerGenerators[i].Z[q])
							++phase;
						
						GZ[q] = XOR(GZ[q], stabilizerGenerators[i].Z[q]);
					}
				}

				return (phase % 4) ? -1.0 : 1.0;
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
			bool CommutesWithDestabilizer(const Generator& g, const Generator& destabilizer, const std::vector<size_t>& pos) const
			{
				bool anticommutes = false;
				for (size_t j = 0; j < pos.size(); ++j)
				{
					const size_t pauliOpQubit = pos[j];
					if (g.X[pauliOpQubit] && destabilizer.Z[pauliOpQubit])
						anticommutes = !anticommutes;
					if (g.Z[pauliOpQubit] && destabilizer.X[pauliOpQubit])
						anticommutes = !anticommutes;
				}

				return !anticommutes;
			}

			bool CheckStabilizersAnticommutation(const Generator& g, const std::vector<size_t>& pos) const
			{
				for (size_t i = 0; i < getNrQubits(); ++i)
				{
					// check if the stabilizer anticommutes with the operator
					bool anticommutes = false;
					for (size_t j = 0; j < pos.size(); ++j)
					{
						const size_t pauliOpQubit = pos[j];
						if (g.X[pauliOpQubit] && stabilizerGenerators[i].Z[pauliOpQubit])
							anticommutes = !anticommutes;
						if (g.Z[pauliOpQubit] && stabilizerGenerators[i].X[pauliOpQubit])
							anticommutes = !anticommutes;
					}
					if (anticommutes)
						return true;
				}

				return false;
			}

			static void SetPauliString(const std::string& pauliString, Generator& g, std::vector<size_t>& pos, size_t& phase)
			{
				for (size_t i = 0; i < pauliString.size(); ++i)
				{
					const char c = toupper(pauliString[i]);
					switch (c)
					{
					case 'I':
						break;
					case 'X':
						g.X[i] = true;
						pos.push_back(i);
						break;
					case 'Y':
						g.X[i] = true;
						g.Z[i] = true;
						// Y = iXZ, so we get a phase factor of i
						++phase;
						pos.push_back(i);
						break;
					case 'Z':
						g.Z[i] = true;
						pos.push_back(i);
						break;
					default:
						throw std::runtime_error("Invalid operator in the Pauli string");
					}
				}
			}

			inline void ApplyH(size_t qubit, size_t q)
			{
				// how does this work?
				// looking at it might not reveal immediately what it does
				// the swaps below just switch the X with Z, because HZH^t = X and HXH^t = Z

				// if we have both X and Z, then it's a Y (with some global phase, given by the sign, Y = iXZ)
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
				// X does nothing to X, but flips the sign for Z, as XZX^t = -Z
				if (destabilizerGenerators[q].Z[qubit])
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				if (stabilizerGenerators[q].Z[qubit])
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
			}

			inline void ApplyY(size_t qubit, size_t q)
			{
				// Y flips the sign for both X and Z, as YZY^t = -Z and YXY^t = -X
				// if both X and Z are present, the sign is flipped twice, so it remains unchanged

				// can be done with ifs, can be done with XORs
				if (XOR(destabilizerGenerators[q].Z[qubit], destabilizerGenerators[q].X[qubit]))
					destabilizerGenerators[q].PhaseSign = !destabilizerGenerators[q].PhaseSign;

				if (XOR(stabilizerGenerators[q].Z[qubit], stabilizerGenerators[q].X[qubit]))
					stabilizerGenerators[q].PhaseSign = !stabilizerGenerators[q].PhaseSign;
			}

			inline void ApplyZ(size_t qubit, size_t q)
			{
				// very similar with applying X, but now Z does nothing to Z, and flips the sign for X, as Z^tXZ = -X
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
