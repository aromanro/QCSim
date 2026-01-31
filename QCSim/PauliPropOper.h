#pragma once

#include "PauliStringXZCoeff.h"

namespace QC
{
	enum class OperationType : unsigned char
	{
		X,
		Y,
		Z,
		H,
		K,
		S,
		SDG,
		SX,
		SXDG,
		CX,
		CY,
		CZ,
		SWAP,
		ISWAP,
		ISWAPDG,
		PROJ,
		RX,
		RY,
		RZ
	};

	class Operator {
	public:
		Operator() : type(OperationType::X), qubits(1, 0) {}

		Operator(OperationType type, int q1 = 0, int q2 = 0, int q3 = 0)
			: type(type), qubits(GetNrQubitsForType(type))
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
		virtual void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& /*pauliString*/, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const
		{
		}

	private:
		static int GetNrQubitsForType(OperationType type)
		{
			if (static_cast<int>(type) >= static_cast<int>(OperationType::CX) && type != OperationType::PROJ)
				return 2;

			return 1;
		}

		OperationType type;
		std::vector<int> qubits;
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

		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& pauliStrings) const override
		{
			if (pauliString->Coefficient == 0.0)
				return;
			
			const int qubit = GetQubit(0);
			if (pauliString->X[qubit]) // X or Y present - P anticommutes with Z
			{
				pauliString->Coefficient = 0.0;
				return;
			}

			pauliString->Coefficient *= coefficient; // <P>

			auto pstrNew = std::make_unique<PauliStringXZWithCoefficient>(*pauliString);
			// +/- {P, Z}
			// I or Z present - P commutes with Z
			pstrNew->Z[qubit] = !pstrNew->Z[qubit]; // Z becomes I, I becomes Z

			if (projectOne) // P1 = (I - Z)/2
				pstrNew->Coefficient *= -1.0;

			pauliStrings.push_back(std::move(pstrNew));
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

		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplyX(static_cast<size_t>(qubit));
		}
	};

	class OperatorY : public Operator {
	public:
		OperatorY(int q1 = 0)
			: Operator(OperationType::Y, q1)
		{
		}

		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplyY(static_cast<size_t>(qubit));
		}
	};

	class OperatorZ : public Operator {
	public:
		OperatorZ(int q1 = 0)
			: Operator(OperationType::Z, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplyZ(static_cast<size_t>(qubit));
		}
	};

	class OperatorH : public Operator {
	public:
		OperatorH(int q1 = 0)
			: Operator(OperationType::H, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplyH(static_cast<size_t>(qubit));
		}
	};

	class OperatorK : public Operator {
	public:
		OperatorK(int q1 = 0)
			: Operator(OperationType::K, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplyK(static_cast<size_t>(qubit));
		}
	};

	class OperatorS : public Operator {
	public:
		OperatorS(int q1 = 0)
			: Operator(OperationType::S, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplySdag(static_cast<size_t>(qubit));
		}
	};

	class OperatorSDG : public Operator {
	public:
		OperatorSDG(int q1 = 0)
			: Operator(OperationType::SDG, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplyS(static_cast<size_t>(qubit));
		}
	};

	class OperatorSX : public Operator {
	public:
		OperatorSX(int q1 = 0)
			: Operator(OperationType::SX, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplySxDag(static_cast<size_t>(qubit));
		}
	};

	class OperatorSXDG : public Operator {
	public:
		OperatorSXDG(int q1 = 0)
			: Operator(OperationType::SXDG, q1)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit = GetQubit(0);
			pauliString->ApplySx(static_cast<size_t>(qubit));
		}
	};

	class OperatorCX : public Operator {
	public:
		OperatorCX(int control = 0, int target = 0)
			: Operator(OperationType::CX, control, target)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int control = GetQubit(0);
			const int target = GetQubit(1);
			pauliString->ApplyCX(static_cast<size_t>(control), static_cast<size_t>(target));
		}
	};

	class OperatorCY : public Operator {
	public:
		OperatorCY(int control = 0, int target = 0)
			: Operator(OperationType::CY, control, target)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int control = GetQubit(0);
			const int target = GetQubit(1);
			pauliString->ApplyCY(static_cast<size_t>(control), static_cast<size_t>(target));
		}
	};

	class OperatorCZ : public Operator {
	public:
		OperatorCZ(int control = 0, int target = 0)
			: Operator(OperationType::CZ, control, target)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int control = GetQubit(0);
			const int target = GetQubit(1);
			pauliString->ApplyCZ(static_cast<size_t>(control), static_cast<size_t>(target));
		}
	};

	class OperatorSWAP : public Operator {
	public:
		OperatorSWAP(int q1 = 0, int q2 = 0)
			: Operator(OperationType::SWAP, q1, q2)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit1 = GetQubit(0);
			const int qubit2 = GetQubit(1);
			pauliString->ApplySwap(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
		}
	};

	class OperatorISWAP : public Operator {
	public:
		OperatorISWAP(int q1 = 0, int q2 = 0)
			: Operator(OperationType::ISWAP, q1, q2)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit1 = GetQubit(0);
			const int qubit2 = GetQubit(1);
			pauliString->ApplyISwapDag(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
		}
	};

	class OperatorISWAPDG : public Operator {
	public:
		OperatorISWAPDG(int q1 = 0, int q2 = 0)
			: Operator(OperationType::ISWAPDG, q1, q2)
		{
		}
		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& /*pauliStrings*/) const override
		{
			const int qubit1 = GetQubit(0);
			const int qubit2 = GetQubit(1);
			pauliString->ApplyISwap(static_cast<size_t>(qubit1), static_cast<size_t>(qubit2));
		}
	};

	class OperatorRotation : public Operator {
	public:
		OperatorRotation(OperationType type, int q1 = 0, double angle = 0.0)
			: Operator(type, q1), angle(angle)
		{
		}

		double GetAngle() const
		{
			return angle;
		}

	private:
		double angle;
	};

	class OperatorRZ : public OperatorRotation {
	public:
		OperatorRZ(int q1 = 0, double angle = 0.0)
			: OperatorRotation(OperationType::RZ, q1, angle)
		{
		}

		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& pauliStrings) const override
		{
			if (pauliString->Coefficient == 0.0)
				return;

			const int qubit = GetQubit(0);
			// if I or Z, nothing changes
			if (!pauliString->X[qubit])
				return;

			// the Pauli string is split in two, make a copy for the second term
			std::unique_ptr<PauliStringXZWithCoefficient> pstrNew = std::make_unique<PauliStringXZWithCoefficient>(*pauliString);

			const double ang = GetAngle();
			// the first term is multiplied by cos(angle) and preserves X or Y on the qubit position, so we're done with it
			pauliString->Coefficient *= std::cos(ang);

			// now deal with the second term
			// X is set, check Y
			if (pauliString->Z[qubit]) // Y present
			{
				pstrNew->Coefficient *= std::sin(ang);
				pstrNew->Z[qubit] = false; // Y becomes X
			}
			else // only X present
			{
				pstrNew->Coefficient *= -std::sin(ang);
				pstrNew->Z[qubit] = true; // X becomes Y	
			}
			pauliStrings.push_back(std::move(pstrNew));
		}
	};


	class OperatorRX : public OperatorRotation {
	public:
		OperatorRX(int q1 = 0, double angle = 0.0)
			: OperatorRotation(OperationType::RX, q1, angle)
		{
		}

		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& pauliStrings) const override
		{
			if (pauliString->Coefficient == 0.0)
				return;

			const int qubit = GetQubit(0);

			// if I or X, nothing changes
			if (!pauliString->Z[qubit])
				return;

			// the Pauli string is split in two, make a copy for the second term
			std::unique_ptr<PauliStringXZWithCoefficient> pstrNew = std::make_unique<PauliStringXZWithCoefficient>(*pauliString);

			const double ang = GetAngle();
			// the first term is multiplied by cos(angle) and preserves Z or Y on the qubit position, so we're done with it
			pauliString->Coefficient *= std::cos(ang);

			// now deal with the second term
			// Z is set, check X
			if (pauliString->X[qubit]) // Y present
			{
				pstrNew->Coefficient *= -std::sin(ang);
				pstrNew->X[qubit] = false; // Y becomes Z
			}
			else // only Z present
			{
				pstrNew->Coefficient *= std::sin(ang);
				pstrNew->X[qubit] = true; // Z becomes Y
			}
			pauliStrings.push_back(std::move(pstrNew));
		}
	};

	class OperatorRY : public OperatorRotation {
	public:
		OperatorRY(int q1 = 0, double angle = 0.0)
			: OperatorRotation(OperationType::RY, q1, angle)
		{
		}

		void Apply(std::unique_ptr<PauliStringXZWithCoefficient>& pauliString, std::vector<std::unique_ptr<PauliStringXZWithCoefficient>>& pauliStrings) const override
		{
			if (pauliString->Coefficient == 0.0)
				return;

			const int qubit = GetQubit(0);

			// if I or Y, nothing changes
			if (pauliString->X[qubit] == pauliString->Z[qubit])
				return;

			// the Pauli string is split in two, make a copy for the second term
			std::unique_ptr<PauliStringXZWithCoefficient> pstrNew = std::make_unique<PauliStringXZWithCoefficient>(*pauliString);

			const double ang = GetAngle();
			// the first term is multiplied by cos(angle) and preserves X or Z on the qubit position, so we're done with it
			pauliString->Coefficient *= std::cos(ang);

			// now deal with the second term
			// any can be checked, as only one is set
			if (pauliString->X[qubit]) // X present
			{
				pstrNew->Coefficient *= std::sin(ang);
				// X becomes Z
				pstrNew->X[qubit] = false;
				pstrNew->Z[qubit] = true;
			}
			else // Z case
			{
				pstrNew->Coefficient *= -std::sin(ang);
				// Z becomes X
				pstrNew->X[qubit] = true;
				pstrNew->Z[qubit] = false;
			}
			pauliStrings.push_back(std::move(pstrNew));
		}
	};

}




