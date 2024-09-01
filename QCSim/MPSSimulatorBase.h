#pragma once

#include <vector>
#include <iostream>

#include <unsupported/Eigen/CXX11/Tensor>

#include "MPSSimulatorInterface.h"
#include "Operators.h"

namespace QC {

	namespace TensorNetworks {

		// this is separated from the actual simulator to reduce the class complexity
		// here there are the types definitions, the data structures used and some functions that are simpler and/or not so important for the implementation
		// for example, the code that converts the MPS to a state vector is here, but it wouldn't be needed for a simulation, it's needed just for comparing the results against the statevector simulator
		// also the initialization functions are here
		class MPSSimulatorBase : public MPSSimulatorInterface
		{
		public:
			MPSSimulatorBase() = delete;

			MPSSimulatorBase(size_t N, int addseed = 0)
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
			VectorClass getRegisterStorage() const
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

			void print() const override
			{
				for (size_t i = 0; i < gammas.size() - 1; ++i)
				{
					std::cout << "Lambda " << i << ":\n" << lambdas[i] << std::endl;
					std::cout << "Gamma " << i << ":\n" << gammas[i] << std::endl;
				}

				std::cout << "Gamma " << gammas.size() - 1 << ":\n" << gammas[gammas.size() - 1] << std::endl;
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