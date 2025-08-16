#pragma once

#include <vector>
#include <iostream>

#include <unsupported/Eigen/CXX11/Tensor>

#include "MPSSimulatorInterface.h"
#include "Operators.h"

namespace QC {

	namespace TensorNetworks {

		class MPSSimulatorBaseState : public MPSSimulatorStateInterface
		{
		public:
			MPSSimulatorBaseState() = default;
			MPSSimulatorBaseState(const MPSSimulatorBaseState&) = default;
			MPSSimulatorBaseState(MPSSimulatorBaseState&&) = default;
			MPSSimulatorBaseState& operator=(const MPSSimulatorBaseState&) = default;
			MPSSimulatorBaseState& operator=(MPSSimulatorBaseState&&) = default;
			virtual ~MPSSimulatorBaseState() = default;
			
			std::vector<MPSSimulatorInterface::LambdaType> lambdas;
			std::vector<MPSSimulatorInterface::GammaType> gammas;
		};

		// this is separated from the actual simulator to reduce the class complexity
		// here there are the types definitions, the data structures used and some functions that are simpler and/or not so important for the implementation
		// for example, the code that converts the MPS to a state vector is here, but it wouldn't be needed for a simulation, it's needed just for comparing the results against the statevector simulator
		// also the initialization functions are here
		class MPSSimulatorBase : public MPSSimulatorInterface
		{
		public:
			MPSSimulatorBase() = delete;

			MPSSimulatorBase(size_t N, unsigned int addseed = 0)
				: lambdas(N - 1, LambdaType::Ones(1)), gammas(N, GammaType(1, 2, 1))
			{
				for (auto& gamma : gammas)
				{
					gamma(0, 0, 0) = 1.;
					gamma(0, 1, 0) = 0.;
				}

				const uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + addseed;
				std::seed_seq seed{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
				rng.seed(seed);
			}

			MPSSimulatorBase(const MPSSimulatorBase&) = default;
			MPSSimulatorBase(MPSSimulatorBase&&) = default;
			MPSSimulatorBase& operator=(const MPSSimulatorBase&) = default;
			MPSSimulatorBase& operator=(MPSSimulatorBase&&) = default;

			size_t getNrQubits() const override
			{
				return gammas.size();
			}

			void Clear() override
			{
				const size_t szm1 = lambdas.size();
				for (size_t i = 0; i < szm1; ++i)
				{
					gammas[i].resize(1, 2, 1);
					gammas[i](0, 0, 0) = 1.;
					gammas[i](0, 1, 0) = 0.;

					lambdas[i].resize(1);
					lambdas[i](0) = 1.;
				}

				gammas[szm1].resize(1, 2, 1);
				gammas[szm1](0, 0, 0) = 1.;
				gammas[szm1](0, 1, 0) = 0.;
			}

			void InitOnesState() override
			{
				const size_t szm1 = lambdas.size();
				for (size_t i = 0; i < szm1; ++i)
				{
					gammas[i].resize(1, 2, 1);
					gammas[i](0, 0, 0) = 0.;
					gammas[i](0, 1, 0) = 1.;

					lambdas[i].resize(1);
					lambdas[i](0) = 1.;
				}

				gammas[szm1].resize(1, 2, 1);
				gammas[szm1](0, 0, 0) = 0.;
				gammas[szm1](0, 1, 0) = 1.;
			}

			void setToQubitState(IndexType q) override
			{
				Clear();
				if (q >= static_cast<IndexType>(gammas.size()))
					return;

				gammas[q](0, 0, 0) = 0.;
				gammas[q](0, 1, 0) = 1.;
			}

			void setToBasisState(size_t State) override
			{
				const size_t NrBasisStates = gammas.size() > sizeof(size_t) * 8 ? 64 : (1ULL << gammas.size());
				if (State >= NrBasisStates) return;

				Clear();

				size_t pos = 0;
				while (State)
				{
					if (State & 1)
					{
						gammas[pos](0, 0, 0) = 0.;
						gammas[pos](0, 1, 0) = 1.;
					}
					State >>= 1;
					++pos;
				}
			}

			void setToBasisState(const std::vector<bool>& State) override
			{
				if (State.size() > gammas.size()) return;

				Clear();

				for (size_t i = 0; i < State.size(); ++i)
				{
					if (State[i] == 1)
					{
						gammas[i](0, 0, 0) = 0.;
						gammas[i](0, 1, 0) = 1.;
					}
				}
			}

			double GetProbability(IndexType qubit, bool zeroVal = true) const override
			{
				if (qubit < 0 || qubit >= static_cast<IndexType>(gammas.size()))
					throw std::invalid_argument("Qubit index out of bounds");

				const bool notFirst = qubit > 0;
				const bool notLast = qubit < static_cast<IndexType>(lambdas.size());
				if (notFirst && notLast)
					return GetProbabilityMiddleQubit(qubit, zeroVal);
				else if (notFirst)
					return GetProbabilityLastQubit(qubit, zeroVal);
				else if (notLast)
					return GetProbabilityFirstQubit(qubit, zeroVal);

				assert(qubit == 0);

				return GetProbabilitySingleQubit(zeroVal);
			}

			void setLimitBondDimension(IndexType chival) override
			{
				limitSize = true;
				chi = chival;
			}

			void setLimitEntanglement(double svdThreshold) override
			{
				limitEntanglement = true;
				singularValueThreshold = svdThreshold;
			}

			void dontLimitBondDimension() override
			{
				limitSize = false;
			}

			void dontLimitEntanglement() override
			{
				limitEntanglement = false;
			}

			// this is for 'compatibility' with the statevector simulator (QubitRegister)
			// it's not stored as this and it's costly to compute, it will throw an exception for more than 32 qubits
			// but don't call it for such a large number of qubits
			VectorClass getRegisterStorage() const override
			{
				const size_t sz = gammas.size();
				if (sz > sizeof(size_t) * 4) throw std::runtime_error("Too many qubits to compute the state vector");

				if (sz < 8)
					return getRegisterStorage8(sz);
				else if (sz < 16)
					return getRegisterStorage16(sz);
				else if (sz < 24)
					return getRegisterStorage24(sz);
				else if (sz < 32)
					return getRegisterStorage32(sz);

				return {};
			}

			std::complex<double> getBasisStateAmplitude(size_t State) const override
			{
				std::vector<bool> state(getNrQubits());

				for (size_t i = 0; i < state.size(); ++i)
				{
					state[i] = (State & 1) == 1;
					State >>= 1;
				}

				return getBasisStateAmplitude(state);
			}

			std::complex<double> getBasisStateAmplitude(std::vector<bool>& State) const override
			{
				const size_t nrQubits = getNrQubits();
				if (nrQubits == 0) return 0.;
				State.resize(nrQubits, false);

				static const Indexes product_dims{ IntIndexPair(1, 0) };
				MatrixTensorType res = gammas[0].chip(State[0] ? 1 : 0, 1);

				for (size_t q = 1; q < nrQubits; ++q)
				{
					const size_t q1 = q - 1;

					for (IndexType c = 0; c < res.dimension(1); ++c)
						for (IndexType r = 0; r < res.dimension(0); ++r)
							res(r, c) *= lambdas[q1][c];

					// why? Needs this intermediary variable here, not even calling eval() works if assigning directly to res
					MatrixTensorType tmp = res.contract(gammas[q].chip(State[q] ? 1 : 0, 1), product_dims);
					res = std::move(tmp);
				}

				return res(0, 0);
			}

			double getBasisStateProbability(size_t State) const override
			{
				return std::norm(getBasisStateAmplitude(State));
			}

			double getBasisStateProbability(std::vector<bool>& State) const override
			{
				return std::norm(getBasisStateAmplitude(State));
			}

			std::shared_ptr<MPSSimulatorStateInterface> getState() const override
			{
				auto state = std::make_shared<MPSSimulatorBaseState>();
				state->lambdas = lambdas;
				state->gammas = gammas;

				return state;
			}

			void setState(const std::shared_ptr<MPSSimulatorStateInterface>& state) override
			{
				auto stateRef = std::static_pointer_cast<MPSSimulatorBaseState>(state);
				lambdas = stateRef->lambdas;
				gammas = stateRef->gammas;
			}

			void print() const override
			{
				for (size_t i = 0; i < gammas.size() - 1; ++i)
				{
					std::cout << std::endl << "Gamma " << i << ":" << std::endl;
					PrintGamma(i);
					std::cout << "Lambda " << i << ":\n" << lambdas[i] << std::endl;
				}

				std::cout << std::endl << "Gamma " << gammas.size() - 1 << ":" << std::endl;
				PrintGamma(gammas.size() - 1);
			}

			// needed for some tests, ignore
			void SetSiteMatrices(size_t site, const Eigen::MatrixXcd& matrix0, const Eigen::MatrixXcd& matrix1)
			{
				assert(matrix0.rows() == matrix1.rows());
				assert(matrix0.cols() == matrix1.cols());

				gammas[site].resize(matrix0.rows(), 2, matrix0.cols());


				for (IndexType j = 0; j < matrix0.cols(); ++j)
					for (IndexType i = 0; i < matrix0.rows(); ++i) 
					{
						gammas[site](i, 0, j) = matrix0(i, j);
						gammas[site](i, 1, j) = matrix1(i, j);
					}
			}

			void SetLambdas(size_t pos, const Eigen::VectorXd& lambda)
			{
				assert(pos < lambdas.size());
				lambdas[pos] = lambda;
			}
			// end test functions

			void MoveAtBeginningOfChain(const std::set<IndexType>& qubits) override
			{
				// do nothing, it's here just to provide an implementation
			}

		protected:
			void PrintGamma(size_t i) const
			{
				assert(i < gammas.size());
				assert(gammas[i].dimension(1) == 2);

				std::cout << std::endl << "Phys index 0 matrix: " << std::endl;
				for (IndexType j = 0;  j < gammas[i].dimension(0); ++j)
				{
					for (IndexType k = 0; k < gammas[i].dimension(2); ++k)
						std::cout << gammas[i](j, 0, k) << " ";
					std::cout << std::endl;
				}
				std::cout << "Phys index 1 matrix: " << std::endl;
				for (IndexType j = 0; j < gammas[i].dimension(0); ++j)
				{
					for (IndexType k = 0; k < gammas[i].dimension(2); ++k)
						std::cout << gammas[i](j, 1, k) << " ";
					std::cout << std::endl;
				}
				std::cout << std::endl;
			}

			double GetProbabilitySingleQubit(bool zeroVal = true) const
			{
				double res = 0;

				const size_t physIndex = zeroVal ? 0 : 1;
				for (IndexType j = 0; j < gammas[0].dimension(2); ++j)
					for (IndexType i = 0; i < gammas[0].dimension(0); ++i)
						res += std::norm(gammas[0](i, physIndex, j));

				return res;
			}

			double GetProbabilityMiddleQubit(IndexType qubit, bool zeroVal = true) const
			{
				double res = 0;
				const size_t physIndex = zeroVal ? 0 : 1;

				const IndexType qbit1 = qubit - 1;
				for (IndexType j = 0; j < lambdas[qubit].size(); ++j)
					for (IndexType i = 0; i < lambdas[qbit1].size(); ++i)
						res += std::norm(lambdas[qbit1][i] * lambdas[qubit][j] * gammas[qubit](i, physIndex, j));
				
				return res;
			}

			double GetProbabilityLastQubit(IndexType qubit, bool zeroVal = true) const
			{
				double res = 0;
				const size_t physIndex = zeroVal ? 0 : 1;

				const IndexType qbit1 = qubit - 1;
				for (IndexType j = 0; j < gammas[qubit].dimension(2); ++j)
					for (IndexType i = 0; i < lambdas[qbit1].size(); ++i)
						res += std::norm(lambdas[qbit1][i] * gammas[qubit](i, physIndex, j));
				
				return res;
			}

			double GetProbabilityFirstQubit(IndexType qubit, bool zeroVal = true) const
			{
				double res = 0;
				const size_t physIndex = zeroVal ? 0 : 1;

				for (IndexType j = 0; j < lambdas[qubit].size(); ++j)
					for (IndexType i = 0; i < gammas[qubit].dimension(0); ++i)
						res += std::norm(lambdas[qubit][j] * gammas[qubit](i, physIndex, j));

				return res;
			}

			void ApplySingleQubitGate(const GateClass& gate, IndexType qubit)
			{
				ApplySingleQubitGate(gammas[qubit], gate);
			}

			void MultiplyMatrixWithLambda(IndexType qubit, MatrixClass& mat) const
			{
				const size_t nrQubits = gammas.size();

				if (qubit != static_cast<IndexType>(nrQubits) - 1)	
					for (IndexType col = 0; col < mat.cols(); ++col)
						for (IndexType row = 0; row < mat.rows(); ++row)
							mat(row, col) *= col < lambdas[qubit].size() ? lambdas[qubit][col] : 0.;
			}

			static void ApplySingleQubitGate(GammaType& gamma, const GateClass& gate)
			{
				ApplySingleQubitGate(gamma, gate.getRawOperatorMatrix());
			}

			static void ApplySingleQubitGate(GammaType& gamma, const MatrixClass& opMat)
			{
				// contract the gate tensor with the qubit tensor
				static const Indexes product_dims1{ IntIndexPair(1, 1) };
				static const std::array<int, 3> permute{ 0, 2, 1 };
				gamma = gamma.contract(Eigen::TensorMap<const OneQubitGateTensor>(opMat.data(), opMat.rows(), opMat.cols()), product_dims1).shuffle(permute);
			}

		private:
			template<int N> static Eigen::Tensor<std::complex<double>, N + 2> ContractNQubits(const Eigen::Tensor<std::complex<double>, N + 1>& left, const LambdaType& lambdaVal, const GammaType& nextQubit)
			{
				const IndexType dim1 = left.dimension(N);
				const IndexType dim2 = nextQubit.dimension(0);

				Eigen::Tensor<std::complex<double>, 2> lambdaTensor(dim1, dim2);
				lambdaTensor.setZero();

				for (IndexType i = 0; i < std::min(static_cast<IndexType>(lambdaVal.size()), std::min(dim1, dim2)); ++i)
					lambdaTensor(i, i) = lambdaVal(i);

				static const Indexes productDim{ IntIndexPair(N, 0) };

				return left.contract(lambdaTensor, productDim).contract(nextQubit, productDim);
			}

			template<int N> Eigen::Tensor<std::complex<double>, N + 2> GetContractedTensor() const
			{
				return ContractNQubits<N>(GetContractedTensor<N - 1>(), lambdas[N - 2], gammas[N - 1]);
			}

			template<int N> static VectorClass GenerateStatevector(const Eigen::Tensor<std::complex<double>, N + 2>& tensor)
			{
				const size_t NrBasisStates = 1ULL << N;
				VectorClass res(NrBasisStates);

				// index for tensor
				std::array<IndexType, N + 2> indices;
				indices[0] = 0;
				indices[N + 1] = 0;

				for (size_t state = 0; state < NrBasisStates; ++state)
				{
					size_t tmp = state;

					for (size_t q = 1; q <= N; ++q)
					{
						indices[q] = tmp & 1;
						tmp >>= 1;
					}

					res(state) = tensor(indices);
				}

				return res;
			}

			template<int N> VectorClass GenerateStatevector() const
			{
				return GenerateStatevector<N>(GetContractedTensor<N>());
			}

			inline VectorClass getRegisterStorage8(size_t sz) const
			{
				switch (sz)
				{
				case 0:
					return {};
				case 1:
					return GenerateStatevector<1>();
				case 2:
					return GenerateStatevector<2>();
				case 3:
					return GenerateStatevector<3>();
				case 4:
					return GenerateStatevector<4>();
				case 5:
					return GenerateStatevector<5>();
				case 6:
					return GenerateStatevector<6>();
				case 7:
					return GenerateStatevector<7>();
				default:
					break;
				}

				return {};
			}

			inline VectorClass getRegisterStorage16(size_t sz) const
			{
				switch (sz)
				{
				case 8:
					return GenerateStatevector<8>();
				case 9:
					return GenerateStatevector<9>();
				case 10:
					return GenerateStatevector<10>();
				case 11:
					return GenerateStatevector<11>();
				case 12:
					return GenerateStatevector<12>();
				case 13:
					return GenerateStatevector<13>();
				case 14:
					return GenerateStatevector<14>();
				case 15:
					return GenerateStatevector<15>();
				default:
					break;
				}

				return {};
			}

			inline VectorClass getRegisterStorage24(size_t sz) const
			{
				switch (sz)
				{
				case 16:
					return GenerateStatevector<16>();
				case 17:
					return GenerateStatevector<17>();
				case 18:
					return GenerateStatevector<18>();
				case 19:
					return GenerateStatevector<19>();
				case 20:
					return GenerateStatevector<20>();
				case 21:
					return GenerateStatevector<21>();
				case 22:
					return GenerateStatevector<22>();
				case 23:
					return GenerateStatevector<23>();
				default:
					break;
				}

				return {};
			}

			inline VectorClass getRegisterStorage32(size_t sz) const
			{
				switch (sz)
				{
				case 24:
					return GenerateStatevector<24>();
				case 25:
					return GenerateStatevector<25>();
				case 26:
					return GenerateStatevector<26>();
				case 27:
					return GenerateStatevector<27>();
				case 28:
					return GenerateStatevector<28>();
				case 29:
					return GenerateStatevector<29>();
				case 30:
					return GenerateStatevector<30>();
				case 31:
					return GenerateStatevector<31>();
				default:
					break;
				}

				return {};
			}

		protected:
			bool limitSize = false;
			bool limitEntanglement = false;
			IndexType chi = 10; // if limitSize is true
			double singularValueThreshold = 0.; // if limitEntanglement is true

			std::vector<LambdaType> lambdas;
			std::vector<GammaType> gammas;

			std::mt19937_64 rng;
			std::uniform_real_distribution<double> uniformZeroOne{ 0, 1 };

			const Operators::ZeroProjection<MatrixClass> zeroProjection;
			const Operators::OneProjection<MatrixClass> oneProjection;
		};

		template<> MPSSimulatorBase::GammaType MPSSimulatorBase::GetContractedTensor<1>() const
		{
			return gammas[0];
		}
	}
}