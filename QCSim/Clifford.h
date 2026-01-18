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


			void ApplyISwapDag(size_t qubit1, size_t qubit2)
			{
				const size_t nrQubits = getNrQubits();
				if (!enableMultithreading || nrQubits < 512)
				{
					for (size_t q = 0; q < nrQubits; ++q)
					{
						ApplyH(qubit2, q);
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyH(qubit1, q);
						ApplyZ(qubit2, q);
						ApplyS(qubit2, q);
						ApplyZ(qubit1, q);
						ApplyS(qubit1, q);
					}
				}
				else
				{
					//const auto processor_count = QC::QubitRegisterCalculator<>::GetNumberOfThreads();

#pragma omp parallel for 
					//num_threads(processor_count) schedule(static, 128)
					for (long long int q = 0; q < static_cast<long long int>(nrQubits); ++q)
					{
						ApplyH(qubit2, q);
						ApplyCX(qubit1, qubit2, q);
						ApplyCX(qubit2, qubit1, q);
						ApplyH(qubit1, q);
						ApplyZ(qubit2, q);
						ApplyS(qubit2, q);
						ApplyZ(qubit1, q);
						ApplyS(qubit1, q);
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
						
						GZ[q] = GZ[q] != stabilizerGenerators[i].Z[q];
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
				destabilizerGenerators[q].ApplyH(qubit);
				stabilizerGenerators[q].ApplyH(qubit);
			}

			inline void ApplyS(size_t qubit, size_t q)
			{
				destabilizerGenerators[q].ApplyS(qubit);
				stabilizerGenerators[q].ApplyS(qubit);
			}

			inline void ApplyX(size_t qubit, size_t q)
			{
				destabilizerGenerators[q].ApplyX(qubit);
				stabilizerGenerators[q].ApplyX(qubit);
			}

			inline void ApplyY(size_t qubit, size_t q)
			{
				destabilizerGenerators[q].ApplyY(qubit);
				stabilizerGenerators[q].ApplyY(qubit);
			}

			inline void ApplyZ(size_t qubit, size_t q)
			{
				destabilizerGenerators[q].ApplyZ(qubit);
				stabilizerGenerators[q].ApplyZ(qubit);
			}

			inline void ApplyCX(size_t target, size_t control, size_t q)
			{
				destabilizerGenerators[q].ApplyCX(target, control);
				stabilizerGenerators[q].ApplyCX(target, control);
			}
		};
	}
}
