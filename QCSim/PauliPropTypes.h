#pragma once

namespace QC {

// Keep the original values stable for callers storing operation identifiers.
enum class OperationType : unsigned char {
    X, Y, Z, H, K, S, SDG, SX, SXDG, CX, CY, CZ, SWAP, ISWAP, ISWAPDG,
    PROJ, RX, RY, RZ,
    U, CU, CRX, CRY, CRZ, CP, CS, CSDAG, CSX, CSXDAG, CH, CCX, CSWAP
};

inline int PauliOperationArity(OperationType type) {
    switch (type) {
    case OperationType::CCX: case OperationType::CSWAP: return 3;
    case OperationType::CX: case OperationType::CY: case OperationType::CZ:
    case OperationType::SWAP: case OperationType::ISWAP: case OperationType::ISWAPDG:
    case OperationType::CU: case OperationType::CRX: case OperationType::CRY:
    case OperationType::CRZ: case OperationType::CP: case OperationType::CS:
    case OperationType::CSDAG: case OperationType::CSX: case OperationType::CSXDAG:
    case OperationType::CH: return 2;
    default: return 1;
    }
}

inline bool IsPauliClifford(OperationType type) {
    switch (type) {
    case OperationType::X: case OperationType::Y: case OperationType::Z:
    case OperationType::H: case OperationType::K: case OperationType::S:
    case OperationType::SDG: case OperationType::SX: case OperationType::SXDG:
    case OperationType::CX: case OperationType::CY: case OperationType::CZ:
    case OperationType::SWAP: case OperationType::ISWAP: case OperationType::ISWAPDG:
        return true;
    default: return false;
    }
}

inline bool IsPauliRotation(OperationType type) {
    return type == OperationType::RX || type == OperationType::RY || type == OperationType::RZ;
}

inline bool IsPauliLocalGate(OperationType type) {
    switch (type) {
    case OperationType::U: case OperationType::CU: case OperationType::CRX:
    case OperationType::CRY: case OperationType::CRZ: case OperationType::CP:
    case OperationType::CS: case OperationType::CSDAG: case OperationType::CSX:
    case OperationType::CSXDAG: case OperationType::CH: case OperationType::CCX:
    case OperationType::CSWAP: return true;
    default: return false;
    }
}

}
