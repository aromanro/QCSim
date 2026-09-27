#pragma once

#include "PauliPropOper.h"
#include "ThreadPool.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <typeinfo>

namespace QC
{
	namespace PauliDetail
	{

		// The public vector<bool> Pauli API is unchanged. Propagation uses word masks.
		template <size_t N> struct Term
		{
			std::array<uint64_t, N> X{}, Z{};
			double Coefficient = 1.;
			explicit Term(size_t = 0) {}
		};
		template <> struct Term<0>
		{
			std::vector<uint64_t> X, Z;
			double Coefficient = 1.;
			explicit Term(size_t qubits = 0) : X((qubits + 63) / 64), Z(X.size()) {}
		};
		inline size_t Popcount(uint64_t x)
		{
#if defined(_MSC_VER) && defined(_M_X64)
			return static_cast<size_t>(__popcnt64(x));
#elif defined(__GNUC__) || defined(__clang__)
			return static_cast<size_t>(__builtin_popcountll(x));
#else
			size_t n = 0;
			for (; x; x &= x - 1)
				++n;
			return n;
#endif
		}
		template <class T> size_t Weight(const T& t)
		{
			size_t result = 0;
			for (size_t w = 0; w < t.X.size(); ++w)
				result += Popcount(t.X[w] | t.Z[w]);
			return result;
		}
		template <class T> double Expectation(const T& t)
		{
			for (auto x : t.X)
				if (x)
					return 0.;
			return t.Coefficient;
		}
		template <class T> T Pack(const PauliStringXZWithCoefficient& p, size_t qubits)
		{
			T t(qubits);
			t.Coefficient = p.Coefficient;
			for (size_t q = 0; q < std::min(qubits, p.X.size()); ++q)
			{
				if (p.X[q])
					t.X[q / 64] |= uint64_t(1) << (q % 64);
				if (p.Z[q])
					t.Z[q / 64] |= uint64_t(1) << (q % 64);
			}
			return t;
		}
		template <class T> PauliStringXZWithCoefficient Unpack(const T& t, size_t qubits)
		{
			PauliStringXZWithCoefficient p(qubits);
			p.Coefficient = t.Coefficient;
			for (size_t q = 0; q < qubits; ++q)
			{
				p.X[q] = (t.X[q / 64] >> (q % 64)) & 1;
				p.Z[q] = (t.Z[q / 64] >> (q % 64)) & 1;
			}
			return p;
		}

		// Recorded operations have no per-gate heap allocation. Unknown user-defined
		// Operator subclasses retain their virtual Apply/Clone behavior.
		struct Operation
		{
			OperationType type;
			int q0, q1;
			double angle = 0., sine = 0., cosine = 1., coefficient = 1.;
			bool projectOne = false;
			std::shared_ptr<const Operator> custom;
			Operation(OperationType t = OperationType::X, int a = 0, int b = 0) : type(t), q0(a), q1(b) {}
			bool Clifford() const
			{
				return !custom && type < OperationType::PROJ;
			}
			static Operation Rotation(OperationType t, int q, double angle)
			{
				Operation op(t, q);
				op.angle = angle;
				constexpr double halfPi = 1.57079632679489661923132169163975144;
				if (std::isfinite(angle) && std::remainder(angle, halfPi) == 0.)
				{
					const int quadrant = static_cast<int>(std::remainder(angle, 4. * halfPi) / halfPi);
					op.sine = quadrant == 1 ? 1. : quadrant == -1 ? -1. : 0.;
					op.cosine = quadrant == 0 ? 1. : (quadrant == 2 || quadrant == -2) ? -1. : 0.;
				}
				else
				{
					op.sine = std::sin(angle);
					op.cosine = std::cos(angle);
				}
				return op;
			}
			std::unique_ptr<Operator> Legacy() const
			{
				if (custom)
					return custom->Clone();
#define QC_PP_OP(name)                                                                                                                     \
	case OperationType::name:                                                                                                              \
		return std::make_unique<Operator##name>(q0)
#define QC_PP_OP2(name)                                                                                                                    \
	case OperationType::name:                                                                                                              \
		return std::make_unique<Operator##name>(q0, q1)
				switch (type)
				{
					QC_PP_OP(X);
					QC_PP_OP(Y);
					QC_PP_OP(Z);
					QC_PP_OP(H);
					QC_PP_OP(K);
					QC_PP_OP(S);
					QC_PP_OP(SDG);
					QC_PP_OP(SX);
					QC_PP_OP(SXDG);
					QC_PP_OP2(CX);
					QC_PP_OP2(CY);
					QC_PP_OP2(CZ);
					QC_PP_OP2(SWAP);
					QC_PP_OP2(ISWAP);
					QC_PP_OP2(ISWAPDG);
				case OperationType::PROJ:
					return std::make_unique<Projector>(q0, projectOne, coefficient);
				case OperationType::RX:
					return std::make_unique<OperatorRX>(q0, angle);
				case OperationType::RY:
					return std::make_unique<OperatorRY>(q0, angle);
				case OperationType::RZ:
					return std::make_unique<OperatorRZ>(q0, angle);
				}
#undef QC_PP_OP
#undef QC_PP_OP2
				throw std::invalid_argument("Unknown Pauli operation");
			}
			static Operation Import(std::unique_ptr<Operator> p)
			{
				if (!p)
					throw std::invalid_argument("Null Pauli operation");
				Operation op(p->GetType(), p->GetQubit(0), p->GetNrQubits() > 1 ? p->GetQubit(1) : 0);
				bool builtin = false;
#define QC_PP_MATCH(name)                                                                                                                  \
	case OperationType::name:                                                                                                              \
		builtin = typeid(*p) == typeid(Operator##name);                                                                                    \
		break
				switch (op.type)
				{
					QC_PP_MATCH(X);
					QC_PP_MATCH(Y);
					QC_PP_MATCH(Z);
					QC_PP_MATCH(H);
					QC_PP_MATCH(K);
					QC_PP_MATCH(S);
					QC_PP_MATCH(SDG);
					QC_PP_MATCH(SX);
					QC_PP_MATCH(SXDG);
					QC_PP_MATCH(CX);
					QC_PP_MATCH(CY);
					QC_PP_MATCH(CZ);
					QC_PP_MATCH(SWAP);
					QC_PP_MATCH(ISWAP);
					QC_PP_MATCH(ISWAPDG);
					QC_PP_MATCH(RX);
					QC_PP_MATCH(RY);
					QC_PP_MATCH(RZ);
				case OperationType::PROJ:
					builtin = typeid(*p) == typeid(Projector);
					break;
				}
#undef QC_PP_MATCH
				if (!builtin)
					op.custom = std::move(p);
				else if (op.type >= OperationType::RX)
					op = Rotation(op.type, op.q0, static_cast<const OperatorRotation&>(*p).GetAngle());
				else if (op.type == OperationType::PROJ)
				{
					const auto& proj = static_cast<const Projector&>(*p);
					op.projectOne = proj.IsProjectOne();
					op.coefficient = proj.GetCoefficient();
				}
				return op;
			}
		};

		// Derive signed local Clifford permutations from the established public kernels
		// once. This also pins the adjoint and target/control conventions in one place.
		inline const std::array<std::array<unsigned char, 16>, 15>& CliffordMaps()
		{
			static const auto maps = []
			{
				std::array<std::array<unsigned char, 16>, 15> result{};
				for (size_t gate = 0; gate < result.size(); ++gate)
				{
					const auto op = Operation(static_cast<OperationType>(gate), 0, 1).Legacy();
					for (unsigned p = 0; p < 16; ++p)
					{
						PauliStringXZWithCoefficient t(2);
						PauliStringStorage unused;
						t.X[0] = p & 1;
						t.Z[0] = p & 2;
						t.X[1] = p & 4;
						t.Z[1] = p & 8;
						op->Apply(t, unused);
						result[gate][p] = static_cast<unsigned char>((t.X[0] ? 1 : 0) | (t.Z[0] ? 2 : 0) | (t.X[1] ? 4 : 0) |
																	 (t.Z[1] ? 8 : 0) | (t.Coefficient < 0 ? 16 : 0));
					}
				}
				return result;
			}();
			return maps;
		}

		template <OperationType Axis, class T>
		void Rotate(const Operation& op, std::vector<T>& terms, size_t begin, size_t end, std::vector<T>& extra)
		{
			if (op.sine == 0. && op.cosine == 1.)
				return;
			const size_t w = static_cast<size_t>(op.q0) / 64;
			const uint64_t bit = uint64_t(1) << (op.q0 % 64);
			for (size_t i = begin; i < end; ++i)
			{
				auto& t = terms[i];
				if (t.Coefficient == 0.)
					continue;
				const bool x = (t.X[w] & bit) != 0, z = (t.Z[w] & bit) != 0;
				const bool anti = Axis == OperationType::RX ? z : Axis == OperationType::RZ ? x : x != z;
				if (!anti)
					continue;
				if (op.sine == 0.)
				{
					t.Coefficient *= op.cosine;
					continue;
				}
				const double sign = Axis == OperationType::RX	? (x ? -1. : 1.)
									: Axis == OperationType::RZ ? (z ? 1. : -1.)
																: (x ? 1. : -1.);
				if (op.cosine == 0.)
				{
					t.Coefficient *= sign * op.sine;
					if constexpr (Axis != OperationType::RZ)
						t.X[w] ^= bit;
					if constexpr (Axis != OperationType::RX)
						t.Z[w] ^= bit;
				}
				else
				{
					const double coefficient = t.Coefficient * sign * op.sine;
					if (coefficient != 0.)
					{
						extra.push_back(t);
						auto& branch = extra.back();
						branch.Coefficient = coefficient;
						if constexpr (Axis != OperationType::RZ)
							branch.X[w] ^= bit;
						if constexpr (Axis != OperationType::RX)
							branch.Z[w] ^= bit;
					}
					t.Coefficient *= op.cosine;
				}
			}
		}

		template <class T>
		void Apply(const Operation& op, std::vector<T>& terms, size_t begin, size_t end, std::vector<T>& extra, size_t qubits)
		{
			if (op.custom)
			{
				for (size_t i = begin; i < end; ++i)
				{
					auto p = Unpack(terms[i], qubits);
					PauliStringStorage output;
					op.custom->Apply(p, output);
					terms[i] = Pack<T>(p, qubits);
					for (const auto& child : output)
						extra.push_back(Pack<T>(child, qubits));
				}
				return;
			}
			const size_t w = static_cast<size_t>(op.q0) / 64;
			const uint64_t bit = uint64_t(1) << (op.q0 % 64);
			switch (op.type)
			{
			case OperationType::X:
				for (size_t i = begin; i < end; ++i)
					if (terms[i].Z[w] & bit)
						terms[i].Coefficient = -terms[i].Coefficient;
				return;
			case OperationType::Z:
				for (size_t i = begin; i < end; ++i)
					if (terms[i].X[w] & bit)
						terms[i].Coefficient = -terms[i].Coefficient;
				return;
			case OperationType::Y:
				for (size_t i = begin; i < end; ++i)
					if ((terms[i].X[w] ^ terms[i].Z[w]) & bit)
						terms[i].Coefficient = -terms[i].Coefficient;
				return;
			case OperationType::H:
				for (size_t i = begin; i < end; ++i)
				{
					auto& t = terms[i];
					if (t.X[w] & t.Z[w] & bit)
						t.Coefficient = -t.Coefficient;
					const uint64_t delta = (t.X[w] ^ t.Z[w]) & bit;
					t.X[w] ^= delta;
					t.Z[w] ^= delta;
				}
				return;
			case OperationType::S:
			case OperationType::SDG:
				for (size_t i = begin; i < end; ++i)
				{
					auto& t = terms[i];
					const uint64_t phaseZ = op.type == OperationType::S ? ~t.Z[w] : t.Z[w];
					if (t.X[w] & phaseZ & bit)
						t.Coefficient = -t.Coefficient;
					t.Z[w] ^= t.X[w] & bit;
				}
				return;
			case OperationType::CX:
			{
				const size_t controlWord = static_cast<size_t>(op.q1) / 64;
				const uint64_t controlBit = uint64_t(1) << (op.q1 % 64);
				for (size_t i = begin; i < end; ++i)
				{
					auto& t = terms[i];
					const bool xc = (t.X[controlWord] & controlBit) != 0, zt = (t.Z[w] & bit) != 0;
					const bool xt = (t.X[w] & bit) != 0, zc = (t.Z[controlWord] & controlBit) != 0;
					if (xc && zt && xt == zc)
						t.Coefficient = -t.Coefficient;
					t.X[w] ^= xc ? bit : 0;
					t.Z[controlWord] ^= zt ? controlBit : 0;
				}
				return;
			}
			case OperationType::RX:
				Rotate<OperationType::RX>(op, terms, begin, end, extra);
				return;
			case OperationType::RY:
				Rotate<OperationType::RY>(op, terms, begin, end, extra);
				return;
			case OperationType::RZ:
				Rotate<OperationType::RZ>(op, terms, begin, end, extra);
				return;
			case OperationType::PROJ:
				for (size_t i = begin; i < end; ++i)
				{
					auto& t = terms[i];
					if (t.X[w] & bit)
					{
						t.Coefficient = 0.;
						continue;
					}
					t.Coefficient *= op.coefficient;
					if (t.Coefficient == 0.)
						continue;
					extra.push_back(t);
					auto& branch = extra.back();
					branch.Z[w] ^= bit;
					if (op.projectOne)
						branch.Coefficient = -branch.Coefficient;
				}
				return;
			default:
				break;
			}
			const bool two = op.type >= OperationType::CX;
			const size_t w1 = two ? static_cast<size_t>(op.q1) / 64 : w;
			const uint64_t bit1 = two ? uint64_t(1) << (op.q1 % 64) : 0;
			const auto& map = CliffordMaps()[static_cast<size_t>(op.type)];
			for (size_t i = begin; i < end; ++i)
			{
				auto& t = terms[i];
				const unsigned input =
					((t.X[w] & bit) ? 1 : 0) | ((t.Z[w] & bit) ? 2 : 0) | ((t.X[w1] & bit1) ? 4 : 0) | ((t.Z[w1] & bit1) ? 8 : 0);
				const unsigned out = map[input];
				if (out & 16)
					t.Coefficient = -t.Coefficient;
				t.X[w] = (t.X[w] & ~bit) | ((out & 1) ? bit : 0);
				t.Z[w] = (t.Z[w] & ~bit) | ((out & 2) ? bit : 0);
				if (two)
				{
					t.X[w1] = (t.X[w1] & ~bit1) | ((out & 4) ? bit1 : 0);
					t.Z[w1] = (t.Z[w1] & ~bit1) | ((out & 8) ? bit1 : 0);
				}
			}
		}

		inline uint64_t Mix(uint64_t x)
		{
			x ^= x >> 30;
			x *= UINT64_C(0xbf58476d1ce4e5b9);
			x ^= x >> 27;
			x *= UINT64_C(0x94d049bb133111eb);
			return x ^ (x >> 31);
		}
		template <class T> size_t Hash(const T& t)
		{
			uint64_t h = UINT64_C(0x9e3779b97f4a7c15);
			for (size_t w = 0; w < t.X.size(); ++w)
				h = Mix(h ^ Mix(t.X[w]) ^ (Mix(t.Z[w]) + UINT64_C(0x9e3779b97f4a7c15)));
			return static_cast<size_t>(h);
		}
		template <class T> struct Workspace
		{
			std::vector<T> terms;
			std::vector<std::vector<T>> extra;
			std::vector<size_t> hashSlots;
		};
		struct Settings
		{
			size_t qubits = 0, weight = std::numeric_limits<size_t>::max();
			double coefficient = 0.;
			int trims = std::numeric_limits<int>::max(), dedup = std::numeric_limits<int>::max();
			size_t parallelThreshold = 16384, batch = 4096, sumThreshold = 65536, sumBatch = 16384;
			ThreadPool<>* pool = nullptr;
			bool Dedup(size_t index) const
			{
				return dedup != std::numeric_limits<int>::max() && index % dedup == 0;
			}
			bool Trim(size_t index) const
			{
				return trims != std::numeric_limits<int>::max() && index % trims == 0;
			}
		};
		template <class T> void Trim(std::vector<T>& terms, const Settings& s)
		{
			size_t write = 0;
			for (size_t read = 0; read < terms.size(); ++read)
			{
				auto& t = terms[read];
				if (std::abs(t.Coefficient) > s.coefficient && (s.weight >= s.qubits || Weight(t) <= s.weight))
				{
					if (read != write)
						terms[write] = std::move(t);
					++write;
				}
			}
			terms.resize(write);
		}
		template <class T> void Deduplicate(Workspace<T>& ws, const Settings& s)
		{
			auto& terms = ws.terms;
			auto& slots = ws.hashSlots;
			const size_t empty = std::numeric_limits<size_t>::max();
			if (slots.empty())
				slots.resize(64);
			std::fill(slots.begin(), slots.end(), empty);
			size_t write = 0;
			for (size_t read = 0; read < terms.size(); ++read)
			{
				auto& t = terms[read];
				if (s.weight < s.qubits && Weight(t) > s.weight)
					continue;
				if ((write + 1) * 2 > slots.size())
				{
					slots.assign(slots.size() * 2, empty);
					for (size_t j = 0; j < write; ++j)
					{
						size_t h = Hash(terms[j]) & (slots.size() - 1);
						while (slots[h] != empty)
							h = (h + 1) & (slots.size() - 1);
						slots[h] = j;
					}
				}
				size_t h = Hash(t) & (slots.size() - 1);
				while (slots[h] != empty && !(terms[slots[h]].X == t.X && terms[slots[h]].Z == t.Z))
					h = (h + 1) & (slots.size() - 1);
				if (slots[h] != empty)
					terms[slots[h]].Coefficient += t.Coefficient;
				else
				{
					slots[h] = write;
					if (read != write)
						terms[write] = std::move(t);
					++write;
				}
			}
			terms.resize(write);
			Trim(terms, s);
		}

		inline size_t Parts(size_t size, size_t threshold, size_t grain, ThreadPool<>* pool)
		{
			if (!pool || pool->IsWorkerThread() || size < threshold)
				return 1;
			return std::max<size_t>(1, std::min(pool->GetThreadCount() + 1, size / grain));
		}
		inline size_t DefaultWorkerCount(size_t hardwareThreads)
		{
			return hardwareThreads > 1 ? hardwareThreads - 1 : 0;
		}
		// Up to four work items per participant balance variable core speeds without
		// allocating a branch buffer/future for every term when batch is very small.
		inline size_t Grain(size_t size, size_t parts, size_t batch)
		{
			return parts == 1 ? std::max<size_t>(1, size) : std::max(batch, 1 + (size - 1) / (parts * 4));
		}
		inline size_t ChunkCount(size_t size, size_t grain)
		{
			return size ? 1 + (size - 1) / grain : 1;
		}
		// Every submitted job is drained before references to the query's scratch die,
		// including allocation/enqueue errors and exceptions thrown by a custom gate.
		// Only participants are enqueued; they pick up indexed work items dynamically.
		// Output buffers belong to work-item indices, never to completion order.
		template <class F> void Chunks(size_t size, size_t parts, ThreadPool<>* pool, F&& fn, size_t grain = 0)
		{
			if (parts == 1)
			{
				fn(0, 0, size);
				return;
			}
			std::vector<std::future<double>> jobs;
			jobs.reserve(parts - 1);
			std::exception_ptr error;
			if (!grain)
				grain = 1 + (size - 1) / parts;
			const size_t count = ChunkCount(size, grain);
			std::atomic<size_t> next{0};
			const auto run = [&]
			{
				std::exception_ptr failure;
				for (;;)
				{
					const size_t job = next.fetch_add(1, std::memory_order_relaxed);
					if (job >= count)
						break;
					const size_t begin = job * grain;
					const size_t end = begin + std::min(grain, size - begin);
					try
					{
						fn(job, begin, end);
					}
					catch (...)
					{
						if (!failure)
							failure = std::current_exception();
					}
				}
				if (failure)
					std::rethrow_exception(failure);
			};
			try
			{
				for (size_t j = 1; j < parts; ++j)
				{
					jobs.push_back(pool->Enqueue(
						[&]
						{
							run();
							return 0.;
						}));
				}
				run();
			}
			catch (...)
			{
				error = std::current_exception();
			}
			for (auto& job : jobs)
				try
				{
					job.get();
				}
				catch (...)
				{
					if (!error)
						error = std::current_exception();
				}
			if (error)
				std::rethrow_exception(error);
		}

		template <class T> double Execute(Workspace<T>& ws, const std::vector<Operation>& operations, const Settings& s)
		{
			auto& terms = ws.terms;
			size_t firstCustom = operations.size();
			bool checkedCustom = false;
			for (size_t next = operations.size(); next && !terms.empty();)
			{
				const size_t high = next - 1;
				size_t low = high;
				if (operations[high].Clifford())
					while (low && !s.Dedup(low) && !s.Trim(low) && operations[low - 1].Clifford())
						--low;
				const auto& op = operations[high];
				const size_t size = terms.size();
				const size_t threshold = op.Clifford() ? std::max<size_t>(1, s.parallelThreshold / (high - low + 1)) : s.parallelThreshold;
				const size_t parts = op.custom ? 1 : Parts(size, threshold, s.batch, s.pool);
				// Avoid chunk-count division on the common sequential path.
				const size_t grain = parts == 1 ? size : Grain(size, parts, s.batch);
				const size_t chunks = parts == 1 ? 1 : ChunkCount(size, grain);
				if (ws.extra.size() < chunks)
					ws.extra.resize(chunks);
				for (size_t j = 0; j < chunks; ++j)
					ws.extra[j].clear();
				const auto applyChunk = [&](size_t job, size_t begin, size_t end)
				{
					auto& extra = ws.extra[job];
					if (!op.Clifford())
						extra.reserve(end - begin);
					for (size_t pos = high + 1; pos > low;)
						Apply(operations[--pos], terms, begin, end, extra, s.qubits);
				};
				// Keep the serial gate loop independent of parallel dispatch and
				// exception-draining machinery, including in compiler inlining decisions.
				if (parts == 1)
					applyChunk(0, 0, size);
				else
					Chunks(size, parts, s.pool, applyChunk, grain);
				if (!op.custom && op.type == OperationType::PROJ)
				{
					// Compact only after every worker has released its vector range.
					// Custom Apply implementations may observe or revive zero terms;
					// retain those terms while such an operation remains ahead.
					if (!checkedCustom)
					{
						firstCustom = static_cast<size_t>(
							std::find_if(operations.begin(), operations.end(), [](const Operation& gate) { return bool(gate.custom); }) -
							operations.begin());
						checkedCustom = true;
					}
					if (high < firstCustom)
						terms.erase(std::remove_if(terms.begin(), terms.end(), [](const T& term) { return term.Coefficient == 0.; }),
									terms.end());
				}
				size_t total = terms.size();
				for (size_t j = 0; j < chunks; ++j)
					total += ws.extra[j].size();
				terms.reserve(total);
				// Original terms followed by branches in input order, regardless of
				// worker count. Global coefficient aggregation precedes truncation.
				for (size_t j = 0; j < chunks; ++j)
					for (auto& child : ws.extra[j])
						terms.push_back(std::move(child));
				if (s.Dedup(low))
					Deduplicate(ws, s);
				else if (s.Trim(low))
					Trim(terms, s);
				next = low;
			}
			const size_t parts = Parts(terms.size(), s.sumThreshold, s.sumBatch, s.pool);
			if (parts == 1)
			{
				double result = 0.;
				for (const auto& t : terms)
					result += Expectation(t);
				return result;
			}
			std::vector<double> sums(parts);
			Chunks(terms.size(), parts, s.pool,
				   [&](size_t job, size_t begin, size_t end)
				   {
					   double sum = 0.;
					   for (size_t i = begin; i < end; ++i)
						   sum += Expectation(terms[i]);
					   sums[job] = sum;
				   });
			double result = 0.;
			for (double sum : sums)
				result += sum;
			return result;
		}

	} // namespace PauliDetail
} // namespace QC
