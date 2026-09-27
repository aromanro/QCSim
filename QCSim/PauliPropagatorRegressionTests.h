#pragma once

#include "PauliPropagator.h"
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
namespace QC
{
	namespace PauliRegression
	{
		inline void Require(bool ok, const char* message)
		{
			if (!ok)
				throw std::runtime_error(message);
		}
		inline void Close(double a, double b, const char* message, double eps = 2e-10)
		{
			Require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) < eps, message);
		}
		using Expansion = std::map<std::string, double>;
		inline Expansion Canonical(const PauliStringStorage& terms)
		{
			Expansion m;
			for (const auto& t : terms)
				m[t.PauliStringXZ::ToString()] += t.Coefficient;
			return m;
		}

		inline std::unique_ptr<Operator> Gate(int code, int a, int b, double angle = .37)
		{
			switch (code)
			{
#define G(n, name)                                                                                                                         \
	case n:                                                                                                                                \
		return std::make_unique<Operator##name>(a)
#define G2(n, name)                                                                                                                        \
	case n:                                                                                                                                \
		return std::make_unique<Operator##name>(a, b)
				G(0, X);
				G(1, Y);
				G(2, Z);
				G(3, H);
				G(4, K);
				G(5, S);
				G(6, SDG);
				G(7, SX);
				G(8, SXDG);
				G2(9, CX);
				G2(10, CY);
				G2(11, CZ);
				G2(12, SWAP);
				G2(13, ISWAP);
				G2(14, ISWAPDG);
#undef G
#undef G2
			case 15:
				return std::make_unique<Projector>(a, true, .73);
			case 16:
				return std::make_unique<OperatorRX>(a, angle);
			case 17:
				return std::make_unique<OperatorRY>(a, angle);
			default:
				return std::make_unique<OperatorRZ>(a, angle);
			}
		}
		inline void Compare(const Expansion& a, const Expansion& b)
		{
			for (const auto& kv : a)
			{
				auto i = b.find(kv.first);
				Close(kv.second, i == b.end() ? 0. : i->second, "expanded coefficient mismatch");
			}
			for (const auto& kv : b)
			{
				auto i = a.find(kv.first);
				Close(kv.second, i == a.end() ? 0. : i->second, "expanded coefficient mismatch");
			}
		}
		inline void LegacyTrim(PauliStringStorage& terms, int width, double cutoff, size_t weight, bool dedup)
		{
			if (dedup)
			{
				std::map<std::string, PauliStringXZWithCoefficient> unique;
				for (const auto& t : terms)
				{
					if (weight < size_t(width) && t.PauliWeight() > weight)
						continue;
					auto key = t.PauliStringXZ::ToString();
					auto i = unique.find(key);
					if (i == unique.end())
						unique.emplace(key, t);
					else
						i->second.Coefficient += t.Coefficient;
				}
				terms.clear();
				for (auto& kv : unique)
					if (std::abs(kv.second.Coefficient) > cutoff)
						terms.push_back(std::move(kv.second));
			}
			else
			{
				terms.erase(
					std::remove_if(terms.begin(), terms.end(), [&](const auto& t)
								   { return std::abs(t.Coefficient) <= cutoff || (weight < size_t(width) && t.PauliWeight() > weight); }),
					terms.end());
			}
		}
		inline void LegacyRun(PauliStringStorage& terms, const std::vector<std::unique_ptr<Operator>>& ops, int width, int dedup, int trims,
							  double cutoff, size_t weight)
		{
			for (int i = int(ops.size()) - 1; i >= 0; --i)
			{
				const size_t count = terms.size();
				terms.reserve(2 * count);
				for (size_t j = 0; j < count; ++j)
					ops[i]->Apply(terms[j], terms);
				if (dedup && i % dedup == 0)
					LegacyTrim(terms, width, cutoff, weight, true);
				else if (trims && i % trims == 0)
					LegacyTrim(terms, width, cutoff, weight, false);
			}
		}
		inline void Kernels()
		{
			for (int width : {2, 32, 63, 64, 65, 127, 128, 129, 255, 256, 257})
				for (int gate = 0; gate < 19; ++gate)
					for (double angle :
						 {0., .37, -1.5707963267948966, 1.5707963267948966, 3.141592653589793, std::nextafter(1.5707963267948966, 2.)})
						for (int bits = 0; bits < 16; ++bits)
						{
							int a = width - 1, b = width / 2 - 1;
							auto op = Gate(gate, a, b, angle);
							std::vector<std::unique_ptr<Operator>> ops;
							ops.push_back(op->Clone());
							PauliPropagator p;
							p.SetNrQubits(width);
							p.SetOperations(std::move(ops));
							PauliStringXZWithCoefficient t(width);
							t.X[a] = bits & 1;
							t.Z[a] = bits & 2;
							t.X[b] = bits & 4;
							t.Z[b] = bits & 8;
							PauliStringStorage expected{t}, actual;
							expected.reserve(2);
							op->Apply(expected[0], expected);
							p.ExpectationValue(t, actual);
							Compare(Canonical(expected), Canonical(actual));
						}
			std::cout << "PASS all local gate expansions across 32/64/128-bit boundaries and 257 qubits\n";
		}
		inline void RandomPropagation()
		{
			std::mt19937 gen(42);
			for (int trial = 0; trial < 80; ++trial)
			{
				int width = trial % 5 == 0 ? 129 : 6;
				PauliPropagator p;
				p.SetNrQubits(width);
				const double cutoff = trial % 3 == 0 ? .002 : 0.;
				const size_t weight = trial % 4 == 0 ? 3 : size_t(-1);
				p.SetStepsBetweenDeduplication(3);
				p.SetStepsBetweenTrims(2);
				p.SetCoefficientThreshold(cutoff);
				p.SetPauliWeightThreshold(weight);
				std::vector<std::unique_ptr<Operator>> ops;
				for (int j = 0; j < 28; ++j)
				{
					int a = int(gen() % width), b = int(gen() % width);
					if (a == b)
						b = (b + 1) % width;
					ops.push_back(Gate(int(gen() % 19), a, b, .1 + (gen() % 1000) * .001));
				}
				p.SetOperations(std::move(ops));
				auto legacyOps = p.GetOperations();
				PauliStringXZWithCoefficient initial(width);
				for (int j = 0; j < 6; ++j)
				{
					int q = gen() % width;
					initial.X[q] = gen() % 2;
					initial.Z[q] = gen() % 2;
				}
				PauliStringStorage expected{initial};
				LegacyRun(expected, legacyOps, width, 3, 2, cutoff, weight);
				for (int workers : {0, 1, 4})
				{
					if (workers)
						p.EnableParallel(workers);
					else
						p.DisableParallel();
					p.SetParallelThreshold(1);
					p.SetBatchSize(1);
					p.SetParallelThresholdForSum(1);
					p.SetBatchSizeForSum(1);
					PauliStringStorage actual;
					p.ExpectationValue(initial, actual);
					Compare(Canonical(expected), Canonical(actual));
				}
			}
			std::cout << "PASS random reference propagation, projectors and cutoffs with serial/1/4 workers\n";
		}
		inline void StatevectorChecks()
		{
			using C = std::complex<double>;
			const C I(0, 1);
			std::mt19937 gen(99);
			for (int trial = 0; trial < 60; ++trial)
			{
				constexpr int n = 4;
				std::vector<C> state(1 << n);
				state[0] = 1.;
				PauliPropagator p;
				p.SetNrQubits(n);
				p.SetStepsBetweenDeduplication(5);
				for (int j = 0; j < 24; ++j)
				{
					int kind = gen() % 5, q = gen() % n, b = (q + 1 + gen() % (n - 1)) % n;
					double a = (int(gen() % 200) - 100) * .013;
					if (kind == 4)
					{
						p.ApplyCX(b, q);
						for (int k = 0; k < (1 << n); ++k)
							if ((k & (1 << b)) && !(k & (1 << q)))
								std::swap(state[k], state[k | (1 << q)]);
						continue;
					}
					C u00, u01, u10, u11;
					if (kind == 0)
					{
						p.ApplyH(q);
						u00 = u01 = u10 = std::sqrt(.5);
						u11 = -std::sqrt(.5);
					}
					if (kind == 1)
					{
						p.ApplyRX(q, a);
						u00 = u11 = std::cos(a / 2);
						u01 = u10 = -I * std::sin(a / 2);
					}
					if (kind == 2)
					{
						p.ApplyRY(q, a);
						u00 = u11 = std::cos(a / 2);
						u01 = -std::sin(a / 2);
						u10 = std::sin(a / 2);
					}
					if (kind == 3)
					{
						p.ApplyRZ(q, a);
						u00 = std::exp(-I * a / 2.);
						u11 = std::exp(I * a / 2.);
						u01 = u10 = 0.;
					}
					for (int k = 0; k < (1 << n); ++k)
						if (!(k & (1 << q)))
						{
							auto x = state[k], y = state[k | (1 << q)];
							state[k] = u00 * x + u01 * y;
							state[k | (1 << q)] = u10 * x + u11 * y;
						}
				}
				for (int check = 0; check < 12; ++check)
				{
					std::string observable(n, 'I');
					for (char& c : observable)
						c = "IXYZ"[gen() % 4];
					C expected = 0.;
					for (int k = 0; k < (1 << n); ++k)
					{
						int out = k;
						C phase = 1.;
						for (int q = 0; q < n; ++q)
						{
							const bool bit = k & (1 << q);
							if (observable[q] == 'X' || observable[q] == 'Y')
								out ^= 1 << q;
							if (observable[q] == 'Y')
								phase *= bit ? -I : I;
							else if (observable[q] == 'Z' && bit)
								phase = -phase;
						}
						expected += std::conj(state[out]) * phase * state[k];
					}
					Close(p.ExpectationValue(observable), expected.real(), "statevector expectation");
				}
				for (size_t outcome = 0; outcome < state.size(); ++outcome)
					Close(p.Probability(outcome), std::norm(state[outcome]), "statevector probability");
			}
			std::cout << "PASS independent statevector expectations and all outcome probabilities\n";
		}
		inline void Sampling()
		{
			for (int width : {6, 65, 129, 257})
				for (int workers : {0, 4})
					for (size_t limit : {size_t(0), size_t(1), size_t(8), size_t(8192)})
					{
						PauliPropagator p;
						p.SetNrQubits(width);
						p.SetStepsBetweenDeduplication(3);
						if (workers)
							p.EnableParallel(workers);
						p.ApplyH(0);
						p.ApplyCX(0, width - 1);
						p.ApplyRY(2, .43);
						p.ApplyRZ(width - 1, .71);
						std::vector<int> qubits{width - 1, 2, 0};
						p.SetSamplingCacheMaxNodes(limit);
						for (bool approximate : {false, true})
						{
							p.SetCoefficientThreshold(approximate ? .015 : 0.);
							p.SetPauliWeightThreshold(approximate ? 2 : size_t(-1));
							p.SetStepsBetweenTrims(2);
							p.SetSeed(987);
							std::unordered_map<std::vector<bool>, size_t> expected;
							for (int shot = 0; shot < 300; ++shot)
								++expected[p.Sample(qubits)];
							p.SetSeed(987);
							Require(p.SampleCounts(qubits, 300) == expected, "seeded cached/repeated sampling mismatch");
						}
						Require(p.SampleCounts(qubits, 0).empty(), "zero shots");
						Require(p.SampleCounts({}, 3).at({}) == 3, "empty sample subset");
					}
			PauliPropagator ghz;
			ghz.SetNrQubits(65);
			ghz.ApplyH(0);
			ghz.ApplyCX(0, 64);
			ghz.SetSeed(777);
			const auto counts = ghz.SampleCounts({64, 0}, 5000);
			Require(counts.size() == 2, "GHZ correlation");
			for (const auto& kv : counts)
			{
				Require(kv.first[0] == kv.first[1], "GHZ bits");
				Require(kv.second > 2300 && kv.second < 2700, "GHZ frequencies");
			}
			ghz.SaveState();
			ghz.Measure({0});
			auto collapse = ghz.SampleCounts({64, 0}, 100);
			Require(collapse.size() == 1, "projector sampling");
			ghz.RestoreState();
			Require(ghz.SampleCounts({64, 0}, 100).size() == 2, "sampling invalidation after restore");
			std::cout << "PASS seeded batch sampling, bounded/disabled caches, cutoffs, large widths and measurement correlations\n";
		}
		inline void RegressionAndLifetime()
		{
			for (int workers : {0, 1, 4})
			{
				PauliPropagator p;
				p.SetNrQubits(9);
				if (workers)
					p.EnableParallel(workers);
				p.SetParallelThreshold(256);
				p.SetBatchSize(128);
				p.SetCoefficientThreshold(.002);
				p.SetStepsBetweenDeduplication(10);
				for (int q = 0; q < 9; ++q)
					p.ApplyRY(q, .5);
				Close(p.ExpectationValue("XXXXXXXXX"), 0., "spawn checkpoint regression");
			}
			PauliPropagator p;
			p.SetNrQubits(1);
			p.EnableParallel(1);
			p.SetParallelThreshold(1);
			p.SetBatchSize(1);
			p.SetParallelThresholdForSum(1);
			p.SetBatchSizeForSum(1);
			p.ApplyRY(0, .5);
			Close(p.ExpectationValue("X"), std::sin(.5), "one-worker nested-sum regression");
			p.EnableParallel(1);
			Require(p.GetThreadCount() == 1, "worker budget");
			p.SaveState();
			PauliPropagator clone;
			clone.SetNrQubits(1);
			clone.ShareOperationsFrom(p);
			clone.SetSavePosition(p.GetSavePosition());
			clone.ApplyX(0);
			Close(p.ExpectationValue("Z"), std::cos(.5), "clone modified source");
			clone.RestoreState();
			Close(clone.ExpectationValue("Z"), p.ExpectationValue("Z"), "clone restore");
			ThreadPool<> pool(2);
			std::atomic<int> finished{0};
			bool caught = false;
			try
			{
				PauliDetail::Chunks(9, 3, &pool,
									[&](size_t job, size_t, size_t)
									{
										++finished;
										if (job == 1)
											throw std::runtime_error("expected worker exception");
									});
			}
			catch (const std::runtime_error&)
			{
				caught = true;
			}
			Require(caught && finished == 3, "failed to drain exceptional jobs");
			auto nested = pool.Enqueue(
				[&]
				{
					Require(PauliDetail::Parts(100, 1, 1, &pool) == 1, "nested scheduling guard");
					return 1.;
				});
			Close(nested.get(), 1., "pool reuse");
			std::cout << "PASS checkpoint/deadlock regressions, copy-on-write, exception draining and pool reuse\n";
		}
		class CustomScale : public OperatorX
		{
		public:
			explicit CustomScale(std::shared_ptr<int> clones) : OperatorX(0), clones(std::move(clones)) {}
			void Apply(PauliStringXZWithCoefficient& p, PauliStringStorage&) const override
			{
				p.Coefficient *= .25;
			}
			std::unique_ptr<Operator> Clone() const override
			{
				++*clones;
				return std::make_unique<CustomScale>(*this);
			}

		private:
			std::shared_ptr<int> clones;
		};
		inline void CustomAndConcurrent()
		{
			PauliPropagator p;
			p.SetNrQubits(1);
			auto clones = std::make_shared<int>(0);
			std::vector<std::unique_ptr<Operator>> ops;
			ops.push_back(std::make_unique<CustomScale>(clones));
			p.SetOperations(std::move(ops));
			Close(p.ExpectationValue("Z"), .25, "custom operator subtype lost");
			PauliPropagator copy;
			copy.SetNrQubits(1);
			copy.ShareOperationsFrom(p);
			Require(*clones == 1, "custom clone contract");
			Close(copy.ExpectationValue("Z"), .25, "custom cloned behavior");
			PauliPropagator shared;
			shared.SetNrQubits(10);
			shared.EnableParallel(2);
			shared.SetStepsBetweenDeduplication(3);
			shared.SetParallelThreshold(1);
			shared.SetBatchSize(32);
			for (int q = 0; q < 10; ++q)
				shared.ApplyRY(q, .37);
			std::vector<std::future<double>> queries;
			for (int i = 0; i < 8; ++i)
				queries.push_back(std::async(std::launch::async, [&] { return shared.ExpectationValue("XXXXXXXXXX"); }));
			for (auto& query : queries)
				Close(query.get(), std::pow(std::sin(.37), 10), "concurrent readonly query");
			std::cout << "PASS custom operator import/clone and concurrent readonly queries on one pool\n";
		}

		template <class Exception, class F> void Throws(F&& fn, const char* message)
		{
			bool caught = false;
			try
			{
				fn();
			}
			catch (const Exception&)
			{
				caught = true;
			}
			Require(caught, message);
		}

		class ReviveZero : public OperatorX
		{
		public:
			ReviveZero() : OperatorX(0) {}
			void Apply(PauliStringXZWithCoefficient& p, PauliStringStorage&) const override
			{
				if (p.Coefficient == 0.)
				{
					p.Coefficient = 1.;
					p.X[0] = false;
				}
			}
			std::unique_ptr<Operator> Clone() const override
			{
				return std::make_unique<ReviveZero>(*this);
			}
		};
		class ProbabilityEstimate : public OperatorX
		{
		public:
			explicit ProbabilityEstimate(double scale, int failQubit = -1) : OperatorX(0), scale(scale), failQubit(failQubit) {}
			void Apply(PauliStringXZWithCoefficient& p, PauliStringStorage&) const override
			{
				if (failQubit < 0 || p.Z[failQubit])
					p.Coefficient *= scale;
			}
			std::unique_ptr<Operator> Clone() const override
			{
				return std::make_unique<ProbabilityEstimate>(*this);
			}

		private:
			double scale;
			int failQubit;
		};

		inline void ReviewRegressions()
		{
			PauliPropagator source, copy;
			source.SetNrQubits(2);
			copy.SetNrQubits(2);
			source.ApplyX(0);
			source.SaveState();
			source.ApplyH(1);
			copy.ApplyX(1);
			copy.ApplyX(1);
			copy.ApplyX(1);
			copy.SaveState();
			copy.ShareOperationsFrom(source);
			Require(copy.GetSavePosition() == 1, "sharing must copy the source checkpoint");
			copy.RestoreState();
			copy.RestoreState();
			Require(copy.GetOperations().size() == 1 && source.GetOperations().size() == 2, "restore grew or modified a shared circuit");
			Close(copy.ExpectationValue("ZI"), -1., "phantom gate after restore");
			copy.ShareOperationsFrom(copy);
			copy.RestoreState();
			Require(copy.GetOperations().size() == 1, "self sharing checkpoint");
			copy.ClearOperations();
			copy.RestoreState();
			Require(copy.GetSavePosition() == 0 && copy.GetOperations().empty(), "clear checkpoint");

			auto cloneCount = std::make_shared<int>(0);
			std::vector<std::unique_ptr<Operator>> custom;
			custom.push_back(std::make_unique<CustomScale>(cloneCount));
			source.SetOperations(std::move(custom));
			source.SaveState();
			source.ApplyX(1);
			copy.ShareOperationsFrom(source);
			copy.RestoreState();
			Require(*cloneCount == 1 && copy.GetOperations().size() == 1, "custom sharing checkpoint");
			Close(copy.ExpectationValue("ZI"), .25, "custom checkpoint behavior");

			PauliPropagator imports;
			imports.SetNrQubits(2);
			imports.ApplyX(0);
			imports.SaveState();
			for (int gate = 9; gate <= 14; ++gate)
			{
				std::vector<std::unique_ptr<Operator>> invalid;
				invalid.push_back(Gate(3, 0, 1));
				invalid.push_back(Gate(gate, 1, 1));
				Throws<std::invalid_argument>([&] { imports.SetOperations(std::move(invalid)); }, "duplicate gate import accepted");
				Require(imports.GetSavePosition() == 1 && imports.GetOperations().size() == 1, "failed import changed circuit");
				Close(imports.ExpectationValue("ZI"), -1., "failed import changed results");
			}
			PauliPropagator empty;
			Throws<std::out_of_range>([&] { empty.ApplyH(0); }, "gate before width accepted");
			Throws<std::out_of_range>(
				[&]
				{
					auto ops = imports.GetOperations();
					empty.SetOperations(std::move(ops));
				},
				"import before width accepted");
			Throws<std::invalid_argument>([&] { imports.SetNrQubits(0); }, "invalid circuit shrink accepted");
			Throws<std::invalid_argument>([&] { imports.ExpectationValue("III"); }, "oversize observable accepted");
			Throws<std::invalid_argument>([&] { imports.Sample({0, 0}); }, "duplicate sampled qubit accepted");
			for (int value : {0, -1})
			{
				Throws<std::invalid_argument>([&] { imports.SetStepsBetweenTrims(value); }, "trim interval accepted");
				Throws<std::invalid_argument>([&] { imports.SetStepsBetweenDeduplication(value); }, "dedup interval accepted");
				Throws<std::invalid_argument>([&] { imports.SetParallelThreshold(value); }, "parallel threshold accepted");
				Throws<std::invalid_argument>([&] { imports.SetBatchSize(value); }, "batch size accepted");
				Throws<std::invalid_argument>([&] { imports.SetParallelThresholdForSum(static_cast<long long>(value)); },
											  "sum threshold accepted");
				Throws<std::invalid_argument>([&] { imports.SetBatchSizeForSum(static_cast<long>(value)); }, "sum batch size accepted");
			}
			Throws<std::invalid_argument>([&] { imports.SetBatchSize(size_t(0)); }, "unsigned zero batch accepted");
			imports.SetBatchSize(size_t(7));
			imports.SetParallelThreshold(9LL);
			Require(imports.GetBatchSize() == 7 && imports.GetParallelThreshold() == 9, "valid integral settings rejected");
			for (double value : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
				Throws<std::invalid_argument>([&] { imports.SetCoefficientThreshold(value); }, "invalid coefficient accepted");
			std::cout << "PASS save/share/restore, transactional import and argument validation regressions\n";

			for (int workers : {0, 1, 4})
				for (int width : {2, 65, 129, 257})
				{
					PauliPropagator p;
					p.SetNrQubits(width);
					if (workers)
						p.EnableParallel(workers);
					p.SetParallelThreshold(1);
					p.SetBatchSize(1);
					std::vector<std::unique_ptr<Operator>> ops;
					ops.push_back(std::make_unique<OperatorH>(width - 1));
					ops.push_back(std::make_unique<Projector>(0, false, .5));
					p.SetOperations(std::move(ops));
					PauliStringXZWithCoefficient term(width);
					term.X[0] = true;
					PauliStringStorage output(257, term);
					Close(p.ExpectationValue(term, output), 0., "projector zero expectation");
					Require(output.empty(), "projector zero terms were retained");
					ops.clear();
					ops.push_back(std::make_unique<ReviveZero>());
					ops.push_back(std::make_unique<Projector>(0, false, .5));
					p.SetOperations(std::move(ops));
					output.clear();
					Close(p.ExpectationValue(term, output), 1., "zero pruning changed custom operator behavior");
				}
			std::cout << "PASS immediate projector compaction, wide terms and custom zero-term behavior\n";

			Require(PauliDetail::DefaultWorkerCount(0) == 0 && PauliDetail::DefaultWorkerCount(1) == 0 &&
						PauliDetail::DefaultWorkerCount(2) == 1 && PauliDetail::DefaultWorkerCount(32) == 31,
					"automatic worker budget");
			imports.EnableParallel();
			Require(imports.GetThreadCount() == PauliDetail::DefaultWorkerCount(std::thread::hardware_concurrency()),
					"caller omitted from default thread budget");
			imports.EnableParallel(2);
			imports.EnableParallel(2);
			Require(imports.GetThreadCount() == 2, "explicit worker budget changed");
			ThreadPool<> pool(1);
			std::promise<void> tail;
			auto completed = tail.get_future();
			std::atomic<size_t> visited{0};
			PauliDetail::Chunks(
				64, 2, &pool,
				[&](size_t job, size_t begin, size_t end)
				{
					Require(begin == job && end == job + 1, "dynamic grain ignored");
					++visited;
					if (job == 0)
						Require(completed.wait_for(std::chrono::seconds(5)) == std::future_status::ready, "slow chunk blocked other work");
					if (job == 63)
						tail.set_value();
				},
				1);
			Require(visited == 64, "dynamic chunks lost work");
			visited = 0;
			Throws<std::runtime_error>(
				[&]
				{
					PauliDetail::Chunks(
						97, 2, &pool,
						[&](size_t job, size_t, size_t)
						{
							++visited;
							if (job == 14 || job == 25)
								throw std::runtime_error("expected dynamic exception");
						},
						1);
				},
				"dynamic exception lost");
			Require(visited == 97, "dynamic work not drained after exception");
			std::cout << "PASS caller-aware default budget, dynamic pickup and exception draining\n";

			for (double estimate : {2., -2., std::nextafter(1., 2.)})
			{
				PauliPropagator p;
				p.SetNrQubits(2);
				p.SetSeed(7);
				std::vector<std::unique_ptr<Operator>> ops;
				ops.push_back(std::make_unique<ProbabilityEstimate>(estimate));
				p.SetOperations(std::move(ops));
				Close(p.Probability1(0), estimate > 0 ? 0. : 1., "unbounded probability one");
				Close(p.Probability0(0), estimate > 0 ? 1. : 0., "unbounded probability zero");
				const auto counts = p.SampleCounts({0, 1}, 40);
				for (const auto& entry : counts)
					Require(entry.first[0] == (estimate < 0), "invalid bounded sample");
				for (size_t outcome = 0; outcome < 4; ++outcome)
				{
					const double probability = p.Probability(outcome);
					Require(std::isfinite(probability) && probability >= 0. && probability <= 1., "invalid joint probability");
				}
				Require(p.Measure({0})[0] == (estimate < 0), "invalid bounded measurement");
			}
			for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
			{
				PauliPropagator p;
				p.SetNrQubits(2);
				std::vector<std::unique_ptr<Operator>> ops;
				ops.push_back(std::make_unique<ProbabilityEstimate>(invalid, 1));
				p.SetOperations(std::move(ops));
				Throws<std::domain_error>([&] { p.Probability1(1); }, "nonfinite probability accepted");
				Throws<std::domain_error>([&] { p.SampleCounts({0, 1}, 2); }, "nonfinite sample accepted");
				Throws<std::domain_error>([&] { p.Measure({1}); }, "nonfinite measurement accepted");
				Throws<std::domain_error>([&] { p.Probability(0); }, "nonfinite joint probability accepted");
				Require(p.GetOperations().size() == 1, "probability failure left temporary projectors");
			}
			std::cout << "PASS bounded probabilities, nonfinite rejection and probability exception cleanup\n";
		}
		inline void Run()
		{
			Kernels();
			RandomPropagation();
			StatevectorChecks();
			Sampling();
			RegressionAndLifetime();
			CustomAndConcurrent();
			ReviewRegressions();
			std::cout << "ALL CHECKS PASSED\n";
		}
	} // namespace PauliRegression
} // namespace QC
