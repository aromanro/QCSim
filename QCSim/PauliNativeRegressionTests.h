#pragma once
#include "PauliPropagator.h"
#include <map>
#include <future>
#include <Eigen/Dense>
#include <iostream>
#include <limits>
#include <random>

namespace QC { namespace PauliNativeRegression {
using Matrix = Eigen::MatrixXcd;
constexpr double pi = 3.141592653589793238462643383279502884;
inline void Require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
inline Matrix Pauli(unsigned code, unsigned qubits) {
    Matrix p = Matrix::Zero(Eigen::Index(1) << qubits, Eigen::Index(1) << qubits);
    for (unsigned col = 0; col < (1U << qubits); ++col) {
        unsigned row = col;
        std::complex<double> coefficient = 1.;
        for (unsigned q = 0; q < qubits; ++q) {
            const unsigned factor = (code >> (2*q)) & 3, bit = (col >> q) & 1;
            if (factor & 1) row ^= 1U << q;
            if (factor & 2) coefficient *= bit ? -1. : 1.;
            if (factor == 3) coefficient *= std::complex<double>(0.,1.);
        }
        p(row,col) = coefficient;
    }
    return p;
}
inline Matrix U(double theta, double phi, double lambda, double gamma = 0.) {
    Matrix u(2,2);
    u << std::cos(theta/2), -std::sin(theta/2)*std::polar(1., lambda),
        std::sin(theta/2)*std::polar(1., phi), std::cos(theta/2)*std::polar(1., phi+lambda);
    return std::polar(1., gamma) * u;
}
inline Matrix Controlled(const Matrix& u) {
    Matrix m = Matrix::Identity(4,4); m.bottomRightCorner(2,2) = u; return m;
}
inline Matrix Rotation(unsigned axis, double angle) {
    return std::cos(angle/2)*Matrix::Identity(2,2) - std::complex<double>(0.,std::sin(angle/2))*Pauli(axis,1);
}
inline Matrix Fixed(QC::OperationType type) {
    using T = QC::OperationType;
    if (type == T::CS || type == T::CSDAG) {
        Matrix u = Matrix::Identity(2,2); u(1,1) = type == T::CS ? std::complex<double>(0,1) : std::complex<double>(0,-1);
        return Controlled(u);
    }
    if (type == T::CSX || type == T::CSXDAG) {
        Matrix u = std::complex<double>(.5,.5)*Matrix::Identity(2,2) + std::complex<double>(.5,-.5)*Pauli(1,1);
        if (type == T::CSXDAG) u = u.adjoint().eval();
        return Controlled(u);
    }
    if (type == T::CH) return Controlled((Pauli(1,1)+Pauli(2,1))/std::sqrt(2.));
    Matrix u = Matrix::Zero(8,8);
    for (unsigned col=0;col<8;++col) {
        unsigned row=col;
        if (type == T::CCX) { if ((col & 6) == 6) row ^= 1; }
        else if ((col & 4) && ((col & 1) != ((col >> 1) & 1))) row ^= 3;
        u(row,col)=1.;
    }
    return u;
}
inline size_t comparisons = 0;
inline void Check(const std::shared_ptr<const QC::PauliDetail::LocalTransfer>& map, const Matrix& u, unsigned bound) {
    Require(map->maxOutputs <= bound, "nonminimal gate expansion");
    for (unsigned p=0; p < (1U << (2*map->qubits)); ++p) {
        Matrix got = Matrix::Zero(u.rows(),u.cols());
        for (auto i=map->offsets[p];i<map->offsets[p+1];++i)
            got += map->entries[i].coefficient * Pauli(map->entries[i].pauli,map->qubits);
        const Matrix want = u.adjoint()*Pauli(p,map->qubits)*u;
        Require((got-want).cwiseAbs().maxCoeff() < 3e-14, "dense adjoint mismatch");
        const auto first = map->offsets[p];
        const bool unchanged = map->offsets[p+1] - first == 1 && map->entries[first].pauli == p
            && map->entries[first].coefficient == 1.;
        Require(((map->unchanged >> p) & 1U) == unchanged, "wrong unchanged-label mask");
        ++comparisons;
    }
}
inline void Algebra() {
        using namespace QC::PauliDetail;
        using T = QC::OperationType;
        for (auto gate : {T::CS,T::CSDAG,T::CSX,T::CSXDAG,T::CH,T::CCX,T::CSWAP}) {
            const auto map=FixedTransfer(gate);
            Require(map == FixedTransfer(gate), "fixed map is not shared");
            Check(map, Fixed(gate), 4);
        }
        std::mt19937 rng(1349);
        std::uniform_real_distribution<double> angle(-4.,4.);
        for (int trial=0;trial<120;++trial) {
            double t=angle(rng),p=angle(rng),l=angle(rng),g=angle(rng);
            if (trial<64) { t=(trial%8-4)*pi/2; p=(trial/8-4)*pi/2; l=(trial%3-1)*pi/2; g=(trial%5-2)*pi/2; }
            if (trial==64) t=p=l=g=1e-12;
            if (trial==65) t=p=l=g=std::nextafter(pi/2,2.);
            const auto one = CompileU(t,p,l);
            Check(one,U(t,p,l),3);
            if (trial<64) Require(one->clifford, "Clifford U was not recognized");
            Check(CompileCU(t,p,l,g), Controlled(U(t,p,l,g)), 8);
            for (int axis=0;axis<3;++axis) Check(CompileCR(axis,t),Controlled(Rotation(AxisCode(axis),t)),4);
            Matrix phase = Matrix::Identity(2,2); phase(1,1)=std::polar(1.,t);
            Check(CompileCP(t), Controlled(phase), 4);
        }
        for (const auto& map : {CompileCR(2,1e-12),CompileCP(1e-12),CompileCU(0.,0.,1e-12,0.)}) {
            double tiny=0.;
            for (auto i=map->offsets[1];i<map->offsets[2];++i)
                if (map->entries[i].pauli==9) tiny=map->entries[i].coefficient;
            Require(std::abs(tiny/2.5e-25-1.)<1e-12,"small quadratic branch was discarded");
        }
        Require(CompileCU(pi,0.,pi,0.)->clifford, "controlled X was not recognized");
        Require(CompileCP(pi)->clifford, "CZ was not recognized");
        std::cout << "PASS " << comparisons << " independent dense local adjoints, sparsity, Clifford angles and tiny coefficients\n";
}

struct Case {
    OperationType gate;
    double theta = .713, phi = .417, lambda = -.923, gamma = .319;
};
template<class Sim> inline void Record(Sim& sim, const Case& c, const std::array<int,3>& q) {
    using T = OperationType;
    switch(c.gate) {
    case T::U: sim.ApplyU(q[0],c.theta,c.phi,c.lambda,c.gamma); break;
    case T::CU: sim.ApplyCU(q[1],q[0],c.theta,c.phi,c.lambda,c.gamma); break;
    case T::CRX: sim.ApplyCRX(q[1],q[0],c.theta); break;
    case T::CRY: sim.ApplyCRY(q[1],q[0],c.theta); break;
    case T::CRZ: sim.ApplyCRZ(q[1],q[0],c.theta); break;
    case T::CP: sim.ApplyCP(q[1],q[0],c.theta); break;
    case T::CS: sim.ApplyCS(q[1],q[0]); break;
    case T::CSDAG: sim.ApplyCSDAG(q[1],q[0]); break;
    case T::CSX: sim.ApplyCSX(q[1],q[0]); break;
    case T::CSXDAG: sim.ApplyCSXDAG(q[1],q[0]); break;
    case T::CH: sim.ApplyCH(q[1],q[0]); break;
    case T::CCX: sim.ApplyCCX(q[1],q[2],q[0]); break;
    case T::CSWAP: sim.ApplyCSwap(q[2],q[0],q[1]); break;
    default: throw std::logic_error("Not a native gate");
    }
}
inline Matrix CaseMatrix(const Case& c) {
    using T = OperationType;
    switch(c.gate) {
    case T::U: return U(c.theta,c.phi,c.lambda,c.gamma);
    case T::CU: return Controlled(U(c.theta,c.phi,c.lambda,c.gamma));
    case T::CRX: return Controlled(Rotation(1,c.theta));
    case T::CRY: return Controlled(Rotation(3,c.theta));
    case T::CRZ: return Controlled(Rotation(2,c.theta));
    case T::CP: { Matrix m=Matrix::Identity(2,2); m(1,1)=std::polar(1.,c.theta); return Controlled(m); }
    default: return Fixed(c.gate);
    }
}
inline std::vector<Case> Cases() {
    using T = OperationType;
    std::vector<Case> result;
    for(auto gate:{T::U,T::CU,T::CRX,T::CRY,T::CRZ,T::CP,T::CS,T::CSDAG,T::CSX,T::CSXDAG,T::CH,T::CCX,T::CSWAP})
        result.push_back({gate});
    result.push_back({T::CU,pi,0.,pi,0.});
    result.push_back({T::U,pi/2,0.,pi,0.});
    return result;
}
inline void Packed() {
    size_t count=0;
    for (int width:{3,64,65,128,129,256,257,513}) {
        const std::array<int,3> q{{width-1,std::max(1,width/2-1),0}};
        for(const auto& c:Cases()) {
            PauliPropagator sim; sim.SetNrQubits(width); Record(sim,c,q);
            const auto ops=sim.GetOperations();
            Require(ops.size()==1 && ops[0]->GetType()==c.gate,"gate was decomposed");
            PauliPropagator imported; imported.SetNrQubits(width); imported.SetOperations(sim.GetOperations());
            const Matrix u=CaseMatrix(c);
            const unsigned arity=PauliOperationArity(c.gate);
            for(unsigned input=0;input<(1U<<(2*arity));++input) {
                PauliStringXZWithCoefficient term(width);
                for(unsigned j=0;j<arity;++j) { term.X[q[j]]=(input>>(2*j))&1; term.Z[q[j]]=(input>>(2*j+1))&1; }
                if(width>3) term.X[2]=term.Z[2]=true;
                const Matrix want=u.adjoint()*Pauli(input,arity)*u;
                const auto verify=[&](const PauliStringStorage& output) {
                    Matrix actual=Matrix::Zero(u.rows(),u.cols());
                    for(const auto& t:output) {
                        unsigned label=0;
                        for(int bit=0;bit<width;++bit) {
                            bool local=false;
                            for(unsigned j=0;j<arity;++j) if(bit==q[j]) local=true;
                            if(!local) Require(t.X[bit]==term.X[bit] && t.Z[bit]==term.Z[bit],"spectator Pauli changed");
                        }
                        for(unsigned j=0;j<arity;++j) label|=unsigned(t.X[q[j]])<<(2*j) | unsigned(t.Z[q[j]])<<(2*j+1);
                        actual+=t.Coefficient*Pauli(label,arity);
                    }
                    if ((actual-want).cwiseAbs().maxCoeff()>=3e-14)
                        throw std::runtime_error("packed/native import mismatch: width="+std::to_string(width)
                            +" gate="+std::to_string(int(c.gate))+" input="+std::to_string(input));
                };
                PauliStringStorage output;
                sim.ExpectationValue(term,output); verify(output);
                output.clear(); imported.ExpectationValue(term,output); verify(output);
                PauliStringStorage legacy{term}; ops[0]->Apply(legacy[0],legacy); verify(legacy);
                ++count;
            }
        }
    }
    std::cout<<"PASS "<<count<<" packed expansions, legacy round trips and spectator bits through 513 qubits\n";
}
inline Matrix Embed(const Matrix& u, const std::array<int,3>& q, int arity, int width) {
    Matrix m=Matrix::Zero(Eigen::Index(1)<<width,Eigen::Index(1)<<width);
    for(unsigned col=0;col<(1U<<width);++col) {
        unsigned local=0;
        for(int j=0;j<arity;++j) local|=((col>>q[j])&1)<<j;
        for(unsigned r=0;r<(1U<<arity);++r) {
            unsigned row=col;
            for(int j=0;j<arity;++j) row=(row&~(1U<<q[j]))|(((r>>j)&1)<<q[j]);
            m(row,col)=u(r,local);
        }
    }
    return m;
}
inline void CircuitsAndSampling() {
    auto cases=Cases();
    for(int workers:{0,3}) {
        PauliPropagator sim; sim.SetNrQubits(3); sim.SetStepsBetweenDeduplication(1);
        if(workers) { sim.EnableParallel(workers); sim.SetParallelThreshold(1); sim.SetBatchSize(1); }
        Eigen::VectorXcd state=Eigen::VectorXcd::Zero(8); state[0]=1.;
        for(int step=0;step<30;++step) {
            const auto& c=cases[step%cases.size()];
            const std::array<int,3> q{{step%3,(step+1)%3,(step+2)%3}};
            Record(sim,c,q); state=Embed(CaseMatrix(c),q,PauliOperationArity(c.gate),3)*state;
            for(unsigned p=0;p<64;++p) {
                std::string observable(3,'I');
                const char labels[]={'I','X','Z','Y'};
                for(unsigned j=0;j<3;++j) observable[j]=labels[(p>>(2*j))&3];
                const double want=(state.adjoint()*Pauli(p,3)*state)(0,0).real();
                Require(std::abs(sim.ExpectationValue(observable)-want)<2e-12,"mixed circuit mismatch");
            }
        }
        for(size_t outcome=0;outcome<8;++outcome)
            Require(std::abs(sim.Probability(outcome)-std::norm(state[outcome]))<2e-12,"native probability mismatch");
        sim.SaveState();
        PauliPropagator clone; clone.SetNrQubits(3); clone.ShareOperationsFrom(sim); clone.SetStepsBetweenDeduplication(1);
        const double saved=sim.ExpectationValue("YZX");
        sim.ApplyCCX(0,1,2); sim.RestoreState();
        Require(std::abs(sim.ExpectationValue("YZX")-saved)<1e-14,"native restore failed");
        clone.ApplyCSwap(0,1,2);
        Require(sim.GetOperations().size()+1==clone.GetOperations().size(),"native clone mutation leaked");
        sim.SetSeed(93); const auto counts=sim.SampleCounts({2,0,1},128);
        sim.SetSeed(93); std::unordered_map<std::vector<bool>,size_t> repeated;
        for(int shot=0;shot<128;++shot) ++repeated[sim.Sample({2,0,1})];
        Require(counts==repeated,"native batch sampling changed seeded draws");
        Require(sim.GetSavePosition()==30 && sim.GetOperations().size()==30,"sampling changed native circuit");
    }
    // Explicit interval semantics: a native gate is one step and truncation is
    // evaluated after its complete expansion, not inside its old decomposition.
    PauliPropagator p; p.SetNrQubits(3); p.ApplyCCX(0,1,2);
    p.SetStepsBetweenDeduplication(1); p.SetCoefficientThreshold(.51);
    PauliStringStorage output; p.ExpectationValue("XII",output);
    Require(output.empty(),"native cutoff was not applied at gate boundary");
    std::cout<<"PASS mixed dense circuits, threaded execution, probabilities, sampling, snapshots and gate-boundary cutoffs\n";
}
inline void ValidationAndCancellation() {
    PauliPropagator p; p.SetNrQubits(3); p.ApplyH(0); p.SaveState();
    const auto rejects=[&](auto f) {
        bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; }
        Require(rejected,"invalid native operation accepted");
        Require(p.GetOperations().size()==1 && p.GetSavePosition()==1,"invalid gate changed circuit");
    };
    rejects([&]{p.ApplyCCX(0,0,2);}); rejects([&]{p.ApplyCCX(0,1,3);});
    rejects([&]{p.ApplyCSwap(0,1,0);}); rejects([&]{p.ApplyCU(-1,0,0.,0.,0.);});
    rejects([&]{p.ApplyU(0,std::numeric_limits<double>::infinity(),0.,0.);});
    rejects([&]{p.ApplyCU(0,1,0.,0.,0.,std::numeric_limits<double>::quiet_NaN());});
    std::vector<std::unique_ptr<Operator>> invalid;
    invalid.push_back(std::make_unique<OperatorLocal>(OperationType::CCX,0,1,3,PauliDetail::FixedTransfer(OperationType::CCX)));
    rejects([&]{p.SetOperations(std::move(invalid));});
    p.ApplyCCX(0,1,2);
    bool rejected=false; try { p.SetNrQubits(2); } catch(const std::exception&) { rejected=true; }
    Require(rejected && p.GetNrQubits()==3,"native third qubit escaped width validation");
    for(auto gate:{OperationType::CCX,OperationType::CSWAP}) {
        PauliPropagator twice; twice.SetNrQubits(3); twice.SetStepsBetweenDeduplication(1);
        Record(twice,{gate},{0,1,2}); Record(twice,{gate},{0,1,2});
        for(unsigned k=0;k<64;++k) {
            PauliStringXZWithCoefficient in(3);
            for(unsigned q=0;q<3;++q) {in.X[q]=(k>>(2*q))&1; in.Z[q]=(k>>(2*q+1))&1;}
            PauliStringStorage out; twice.ExpectationValue(in,out);
            Require(out.size()==1 && out[0].X==in.X && out[0].Z==in.Z && out[0].Coefficient==1.,"involution left residual terms");
        }
    }
    std::cout<<"PASS transactional native validation, third-qubit bounds and exact involution cancellation\n";
}
template<class T> inline void ParallelBatches(int width, const std::array<int,3>& q) {
    using namespace PauliDetail;
    std::mt19937_64 random(913);
    std::vector<T> input;
    for(int i=0;i<1025;++i) {
        PauliStringXZWithCoefficient p(width);
        p.Coefficient=i%17 ? double(i%13-6)/7. : 0.;
        for(int bit=0;bit<width;++bit) { p.X[bit]=random()&1; p.Z[bit]=random()&1; }
        input.push_back(Pack<T>(p,width));
    }
    ThreadPool<> pool(3);
    for(const auto& c:Cases()) for(int cadence:{0,1}) {
        PauliPropagator sim; sim.SetNrQubits(width); Record(sim,c,q);
        std::vector<Operation> operations;
        for(auto& op:sim.GetOperations()) operations.push_back(Operation::Import(std::move(op)));
        Settings settings; settings.qubits=width;
        if(cadence) settings.dedup=cadence;
        Workspace<T> serial, parallel; serial.terms=parallel.terms=input;
        Execute(serial,operations,settings);
        settings.pool=&pool; settings.parallelThreshold=1; settings.batch=31;
        Execute(parallel,operations,settings);
        Require(serial.terms.size()==parallel.terms.size(),"parallel native branch count changed");
        for(size_t i=0;i<serial.terms.size();++i) {
            const auto& a=serial.terms[i]; const auto& b=parallel.terms[i];
            Require(a.X==b.X && a.Z==b.Z && a.Coefficient==b.Coefficient,"parallel native output order or coefficient changed");
        }
    }
}
class CustomNative : public OperatorLocal {
public:
    CustomNative() : OperatorLocal(OperationType::CCX,0,1,2,PauliDetail::FixedTransfer(OperationType::CCX)) {}
    std::unique_ptr<Operator> Clone() const override { return std::make_unique<CustomNative>(*this); }
    void Apply(PauliStringXZWithCoefficient& p, PauliStringStorage&) const override { p.Coefficient*=2.; }
};
inline void ParallelAndCustom() {
    ParallelBatches<PauliDetail::Term<1>>(64,{63,0,32});
    ParallelBatches<PauliDetail::Term<2>>(65,{64,0,32});
    ParallelBatches<PauliDetail::Term<4>>(129,{128,63,64});
    ParallelBatches<PauliDetail::Term<0>>(513,{512,63,256});
    PauliPropagator p; p.SetNrQubits(3);
    std::vector<std::unique_ptr<Operator>> operations;
    operations.push_back(std::make_unique<CustomNative>()); p.SetOperations(std::move(operations));
    p.ApplyCH(0,2); p.EnableParallel(3); p.SetParallelThreshold(1); p.SetBatchSize(1);
    Require(p.ExpectationValue("ZII")==2.,"native subclass bypassed virtual Apply");
    PauliPropagator clone; clone.SetNrQubits(3); clone.ShareOperationsFrom(p);
    Require(dynamic_cast<CustomNative*>(clone.GetOperations()[0].get())!=nullptr,"native subclass lost clone type");
    Require(clone.ExpectationValue("ZII")==2.,"custom native clone changed behavior");
    p.ApplyCU(1,2,.7,.3,-.2,.1); p.ApplyU(0,.4,.5,.1); p.SetStepsBetweenDeduplication(1);
    const double expected=p.ExpectationValue("XYZ");
    std::vector<std::future<double>> readers;
    for(int i=0;i<4;++i) readers.push_back(std::async(std::launch::async,[&]{ return p.ExpectationValue("XYZ"); }));
    for(auto& reader:readers) Require(reader.get()==expected,"concurrent native queries changed result");
    std::cout<<"PASS deterministic native batches across packed widths, custom subclasses and concurrent queries\n";
}
// The thread-local recording cache must return the table a fresh compile
// gives, share it on a repeat, and never confuse gates, angles or evicted slots.
inline void Recording() {
    using namespace QC::PauliDetail;
    using T = QC::OperationType;
    const auto same = [](const std::shared_ptr<const LocalTransfer>& a, const std::shared_ptr<const LocalTransfer>& b) {
        if (a->qubits != b->qubits || a->offsets != b->offsets || a->entries.size() != b->entries.size()
            || a->unchanged != b->unchanged || a->clifford != b->clifford || a->maxOutputs != b->maxOutputs) return false;
        for (size_t i = 0; i < a->entries.size(); ++i)
            if (a->entries[i].pauli != b->entries[i].pauli || a->entries[i].coefficient != b->entries[i].coefficient) return false;
        return true;
    };
    const auto fresh = [](T type, double a, double b, double c, double d) {
        switch (type) {
        case T::U: return CompileU(a, b, c);
        case T::CU: return CompileCU(a, b, c, d);
        case T::CRX: return CompileCR(0, a);
        case T::CRY: return CompileCR(1, a);
        case T::CRZ: return CompileCR(2, a);
        default: return CompileCP(a);
        }
    };
    std::mt19937 rng(2718);
    std::uniform_real_distribution<double> angle(-4., 4.);
    std::vector<double> angles{0., -0., pi, -pi/2, 1e-12};
    for (int i = 0; i < 59; ++i) angles.push_back(angle(rng));
    size_t checks = 0;
    for (int pass = 0; pass < 2; ++pass)
        for (auto type : {T::U, T::CU, T::CRX, T::CRY, T::CRZ, T::CP})
            for (size_t i = 0; i < angles.size(); ++i) {
                const bool many = type == T::U || type == T::CU;
                const double a = angles[i], b = many ? angles[(i+1) % angles.size()] : 0.;
                const double c = many ? angles[(i+2) % angles.size()] : 0., d = type == T::CU ? angles[(i+3) % angles.size()] : 0.;
                const auto table = ParameterizedTransfer(type, a, b, c, d);
                Require(same(table, fresh(type, a, b, c, d)), "cached table differs from a fresh compile");
                Require(table == ParameterizedTransfer(type, a, b, c, d), "repeated angles did not share the table");
                ++checks;
            }
    Require(ParameterizedTransfer(T::CRX, .5) != ParameterizedTransfer(T::CRY, .5), "gates with equal angles shared a table");
    bool threw = false;
    try { ParameterizedTransfer(T::CP, std::numeric_limits<double>::quiet_NaN()); } catch (const std::invalid_argument&) { threw = true; }
    Require(threw, "non-finite cached angle was accepted");

    // One quarter-turn rule: rounded products k * pi/2 are exact for rotations,
    // controlled phases and the extended stabilizer alike.
    for (int k = -64; k <= 64; ++k) {
        long long turns = 0;
        Require(detail::TryGetQuarterTurns(k * (pi/2), turns) && turns == k, "quarter turn was not recognized");
        Require(Operation::Rotation(T::RZ, 0, k * (pi/2)).Clifford(), "quarter-turn rotation was not Clifford");
        if (k % 2) Require(ParameterizedTransfer(T::CP, k * pi)->clifford, "rounded CZ angle was not Clifford");
    }
    long long turns = 0;
    Require(!detail::TryGetQuarterTurns(std::nextafter(pi/2, 2.), turns), "near quarter turn was snapped");
    std::cout << "PASS " << checks << " cached gate tables, eviction, sharing and the quarter-turn rule\n";
}
inline void Run() { Algebra(); Recording(); Packed(); CircuitsAndSampling(); ValidationAndCancellation(); ParallelAndCustom(); }
}}
