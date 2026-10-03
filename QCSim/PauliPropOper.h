#pragma once

#include <memory>
#include <cmath>
#include <utility>
#include "PauliStringXZCoeff.h"
#include "PauliTransfer.h"

namespace QC
{
	using PauliStringStorage = std::vector<PauliStringXZWithCoefficient>;

	class Operator {
	public:
		virtual ~Operator() = default;

		Operator() : type(OperationType::X), qubits(1, 0) {}

		Operator(OperationType type, int q1 = 0, int q2 = 0, int q3 = 0)
			: type(type), qubits(PauliOperationArity(type))
		{
			qubits[0] = q1;
			if (GetNrQubits() > 1)
				qubits[1] = q2;
			if (GetNrQubits() > 2)
				qubits[2] = q3;
		}

		int GetNrQubits() const
		{
			return static_cast<int>(qubits.size());
		}

		int GetQubit(size_t index) const
		{
			return qubits[index];
		}

		OperationType GetType() const
		{
			return type;
		}

		// the second parameter is for the case when the operator expands the number of pauli strings - should add them at the end
		// the first is changed in place
		virtual void Apply(PauliStringXZWithCoefficient& /*pauliString*/, PauliStringStorage& /*pauliStrings*/) const
		{
		}

		virtual std::unique_ptr<Operator> Clone() const
		{
			return nullptr;
		}

	private:
		OperationType type;
		std::vector<int> qubits;
	};

	// Public operation snapshots preserve the compiled action and share its immutable
	// table. Exact-type imports use the packed kernel; subclasses retain virtual Apply.
	class OperatorLocal : public Operator {
	public:
		OperatorLocal(OperationType type, int q0, int q1, int q2,
			std::shared_ptr<const PauliDetail::LocalTransfer> table)
			: Operator(type, q0, q1, q2), table(std::move(table))
		{
			if (!IsPauliLocalGate(type) || !this->table || this->table->qubits != PauliOperationArity(type))
				throw std::invalid_argument("Invalid local Pauli operation");
		}
		const std::shared_ptr<const PauliDetail::LocalTransfer>& GetTransfer() const { return table; }
		std::unique_ptr<Operator> Clone() const override { return std::make_unique<OperatorLocal>(*this); }
		void Apply(PauliStringXZWithCoefficient& term, PauliStringStorage& extra) const override
		{
			if (term.Coefficient == 0.) return;
			unsigned input = 0;
			for (unsigned q = 0; q < table->qubits; ++q)
				input |= unsigned(term.X[GetQubit(q)]) << (2*q) | unsigned(term.Z[GetQubit(q)]) << (2*q+1);
			const auto first = table->offsets[input], last = table->offsets[input+1];
			const auto set = [&](PauliStringXZWithCoefficient& out, const PauliDetail::LocalTransfer::Entry& e) {
				out.Coefficient *= e.coefficient;
				for (unsigned q = 0; q < table->qubits; ++q) {
					out.X[GetQubit(q)] = (e.pauli >> (2*q)) & 1;
					out.Z[GetQubit(q)] = (e.pauli >> (2*q+1)) & 1;
				}
			};
			// The caller may pass a member of extra as term; retain the source before
			// appending, and finish the in-place update before a possible reallocation.
			const auto original = term;
			if (first == last) { term.Coefficient = 0.; return; }
			set(term, table->entries[first]);
			for (auto i = first+1; i < last; ++i) {
				extra.push_back(original);
				set(extra.back(), table->entries[i]);
			}
		}
	private:
		std::shared_ptr<const PauliDetail::LocalTransfer> table;
	};

	class Projector : public Operator {
	public:
		Projector(int qubit, bool projectOne, double coefficient)
			: Operator(OperationType::PROJ, qubit), projectOne(projectOne), coefficient(coefficient)
		{
		}

		bool IsProjectOne() const
		{
			return projectOne;
		}

		double GetCoefficient() const
		{
			return coefficient;
		}

		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& pauliStrings) const override
		{
			if (pauliString.Coefficient == 0.0)
				return;
			
			const int qubit = GetQubit(0);
			if (pauliString.X[qubit]) // X or Y present - P anticommutes with Z
			{
				pauliString.Coefficient = 0.0;
				return;
			}

			pauliString.Coefficient *= coefficient; // <P>

			auto pstrNew = pauliString;
			// +/- {P, Z}
			// I or Z present - P commutes with Z
			pstrNew.Z[qubit] = !pstrNew.Z[qubit]; // Z becomes I, I becomes Z

			if (projectOne) // P1 = (I - Z)/2
				pstrNew.Coefficient *= -1.0;

			pauliStrings.push_back(std::move(pstrNew));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<Projector>(GetQubit(0), projectOne, coefficient);
		}
		
	private:
		bool projectOne;
		double coefficient;
	};

	class OperatorX : public Operator {
	public:
		OperatorX(int q1 = 0)
			: Operator(OperationType::X, q1)
		{
		}

		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplyX(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorX>(GetQubit(0));
		}
	};

	class OperatorY : public Operator {
	public:
		OperatorY(int q1 = 0)
			: Operator(OperationType::Y, q1)
		{
		}

		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplyY(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorY>(GetQubit(0));
		}
	};

	class OperatorZ : public Operator {
	public:
		OperatorZ(int q1 = 0)
			: Operator(OperationType::Z, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplyZ(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorZ>(GetQubit(0));
		}
	};

	class OperatorH : public Operator {
	public:
		OperatorH(int q1 = 0)
			: Operator(OperationType::H, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplyH(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorH>(GetQubit(0));
		}
	};

	class OperatorK : public Operator {
	public:
		OperatorK(int q1 = 0)
			: Operator(OperationType::K, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplyK(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorK>(GetQubit(0));
		}
	};

	class OperatorS : public Operator {
	public:
		OperatorS(int q1 = 0)
			: Operator(OperationType::S, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplySdag(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorS>(GetQubit(0));
		}
	};

	class OperatorSDG : public Operator {
	public:
		OperatorSDG(int q1 = 0)
			: Operator(OperationType::SDG, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplyS(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorSDG>(GetQubit(0));
		}
	};

	class OperatorSX : public Operator {
	public:
		OperatorSX(int q1 = 0)
			: Operator(OperationType::SX, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplySxDag(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorSX>(GetQubit(0));
		}
	};

	class OperatorSXDG : public Operator {
	public:
		OperatorSXDG(int q1 = 0)
			: Operator(OperationType::SXDG, q1)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString.ApplySx(static_cast<size_t>(qubit));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorSXDG>(GetQubit(0));
		}
	};

	class OperatorCX : public Operator {
	public:
		OperatorCX(int target = 0, int control = 0)
			: Operator(OperationType::CX, target, control)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int target = GetQubit(0);
			const int control = GetQubit(1);
			pauliString.ApplyCX(static_cast<size_t>(target), static_cast<size_t>(control));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorCX>(GetQubit(0), GetQubit(1));
		}
	};

	class OperatorCY : public Operator {
	public:
		OperatorCY(int target = 0, int control = 0)
			: Operator(OperationType::CY, target, control)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int target = GetQubit(0);
			const int control = GetQubit(1);
			pauliString.ApplyCY(static_cast<size_t>(target), static_cast<size_t>(control));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorCY>(GetQubit(0), GetQubit(1));
		}
	};

	class OperatorCZ : public Operator {
	public:
		OperatorCZ(int target = 0, int control = 0)
			: Operator(OperationType::CZ, target, control)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int target = GetQubit(0);
			const int control = GetQubit(1);
			pauliString.ApplyCZ(static_cast<size_t>(target), static_cast<size_t>(control));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorCZ>(GetQubit(0), GetQubit(1));
		}
	};

	class OperatorSWAP : public Operator {
	public:
		OperatorSWAP(int q1 = 0, int q2 = 0)
			: Operator(OperationType::SWAP, q1, q2)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit1 = GetQubit(0);
			const int qubit2 = GetQubit(1);
			pauliString.ApplySwap(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorSWAP>(GetQubit(0), GetQubit(1));
		}
	};

	class OperatorISWAP : public Operator {
	public:
		OperatorISWAP(int q1 = 0, int q2 = 0)
			: Operator(OperationType::ISWAP, q1, q2)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit1 = GetQubit(0);
			const int qubit2 = GetQubit(1);
			pauliString.ApplyISwapDag(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorISWAP>(GetQubit(0), GetQubit(1));
		}
	};

	class OperatorISWAPDG : public Operator {
	public:
		OperatorISWAPDG(int q1 = 0, int q2 = 0)
			: Operator(OperationType::ISWAPDG, q1, q2)
		{
		}
		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& /*pauliStrings*/) const override
		{
			const int qubit1 = GetQubit(0);
			const int qubit2 = GetQubit(1);
			pauliString.ApplyISwap(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorISWAPDG>(GetQubit(0), GetQubit(1));
		}
	};

	class OperatorRotation : public Operator {
	public:
		OperatorRotation(OperationType type, int q1 = 0, double angle = 0.0)
			: Operator(type, q1), angle(angle), sine(std::sin(angle)), cosine(std::cos(angle))
		{
		}

		double GetAngle() const
		{
			return angle;
		}

		double GetSin() const { return sine; }
		double GetCos() const { return cosine; }
	private:
		double angle;
		double sine, cosine;
	};

	class OperatorRZ : public OperatorRotation {
	public:
		OperatorRZ(int q1 = 0, double angle = 0.0)
			: OperatorRotation(OperationType::RZ, q1, angle)
		{
		}

		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& pauliStrings) const override
		{
			if (pauliString.Coefficient == 0.0)
				return;

			const int qubit = GetQubit(0);
			// if I or Z, nothing changes
			if (!pauliString.X[qubit])
				return;

			// the Pauli string is split in two, make a copy for the second term
			PauliStringXZWithCoefficient pstrNew = pauliString;

			// the first term is multiplied by cos(angle) and preserves X or Y on the qubit position, so we're done with it
			pauliString.Coefficient *= GetCos();

			// now deal with the second term
			// X is set, check Y
			if (pauliString.Z[qubit]) // Y present
			{
				pstrNew.Coefficient *= GetSin();
				pstrNew.Z[qubit] = false; // Y becomes X
			}
			else // only X present
			{
				pstrNew.Coefficient *= -GetSin();
				pstrNew.Z[qubit] = true; // X becomes Y	
			}
			pauliStrings.push_back(std::move(pstrNew));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorRZ>(GetQubit(0), GetAngle());
		}
	};


	class OperatorRX : public OperatorRotation {
	public:
		OperatorRX(int q1 = 0, double angle = 0.0)
			: OperatorRotation(OperationType::RX, q1, angle)
		{
		}

		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& pauliStrings) const override
		{
			if (pauliString.Coefficient == 0.0)
				return;

			const int qubit = GetQubit(0);

			// if I or X, nothing changes
			if (!pauliString.Z[qubit])
				return;

			// the Pauli string is split in two, make a copy for the second term
			PauliStringXZWithCoefficient pstrNew = pauliString;

			// the first term is multiplied by cos(angle) and preserves Z or Y on the qubit position, so we're done with it
			pauliString.Coefficient *= GetCos();

			// now deal with the second term
			// Z is set, check X
			if (pauliString.X[qubit]) // Y present
			{
				pstrNew.Coefficient *= -GetSin();
				pstrNew.X[qubit] = false; // Y becomes Z
			}
			else // only Z present
			{
				pstrNew.Coefficient *= GetSin();
				pstrNew.X[qubit] = true; // Z becomes Y
			}
			pauliStrings.push_back(std::move(pstrNew));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorRX>(GetQubit(0), GetAngle());
		}
	};

	class OperatorRY : public OperatorRotation {
	public:
		OperatorRY(int q1 = 0, double angle = 0.0)
			: OperatorRotation(OperationType::RY, q1, angle)
		{
		}

		void Apply(PauliStringXZWithCoefficient& pauliString, PauliStringStorage& pauliStrings) const override
		{
			if (pauliString.Coefficient == 0.0)
				return;

			const int qubit = GetQubit(0);

			// if I or Y, nothing changes
			if (pauliString.X[qubit] == pauliString.Z[qubit])
				return;

			// the Pauli string is split in two, make a copy for the second term
			PauliStringXZWithCoefficient pstrNew = pauliString;

			// the first term is multiplied by cos(angle) and preserves X or Z on the qubit position, so we're done with it
			pauliString.Coefficient *= GetCos();

			// now deal with the second term
			// any can be checked, as only one is set
			if (pauliString.X[qubit]) // X present
			{
				pstrNew.Coefficient *= GetSin();
				// X becomes Z
				pstrNew.X[qubit] = false;
				pstrNew.Z[qubit] = true;
			}
			else // Z case
			{
				pstrNew.Coefficient *= -GetSin();
				// Z becomes X
				pstrNew.X[qubit] = true;
				pstrNew.Z[qubit] = false;
			}
			pauliStrings.push_back(std::move(pstrNew));
		}

		std::unique_ptr<Operator> Clone() const override
		{
			return std::make_unique<OperatorRY>(GetQubit(0), GetAngle());
		}
	};

}
