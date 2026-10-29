#include "closure.hpp"
#include <soplex.h>
#include <iostream>

namespace vor {
inline soplex::Rational to_solver(const Q& q) {soplex::Rational r;mpq_set(r.backend().data(),q.get_mpq_t());return r;}
inline Q from_solver(const soplex::Rational& r) {Q q;mpq_set(q.get_mpq_t(),r.backend().data());return q;}

std::optional<Vec> ClosureOracle::linear_program(const std::vector<Vec>& equalities,const Vec* objective,bool margin) {
        budget.check();if(++budget.lp>budget.max_lp)throw Limit("LP solve limit reached");
        if(budget.trace)std::cerr<<"LP "<<budget.lp<<": "<<equalities.size()<<" equalities, "<<cuts.size()<<" cuts, margin="<<margin<<'\n';
        using namespace soplex;
        int k=frame.n*(frame.n-1)/2;
        if(k==0){for(auto& v:equalities)if(constraint_rhs(v)!=0)return {};for(auto& v:cuts)if(constraint_rhs(v)>0)return {};return Vec{};}
        SoPlex lp;
        lp.setIntParam(SoPlex::VERBOSITY,0);
        lp.setIntParam(SoPlex::READMODE,SoPlex::READMODE_RATIONAL);
        lp.setIntParam(SoPlex::SOLVEMODE,SoPlex::SOLVEMODE_RATIONAL);
        lp.setIntParam(SoPlex::CHECKMODE,SoPlex::CHECKMODE_RATIONAL);
        lp.setIntParam(SoPlex::SYNCMODE,SoPlex::SYNCMODE_AUTO);
        lp.setIntParam(SoPlex::OBJSENSE,SoPlex::OBJSENSE_MAXIMIZE);
        lp.setRealParam(SoPlex::FEASTOL,0);lp.setRealParam(SoPlex::OPTTOL,0);
        if(budget.seconds){double used=std::chrono::duration<double>(std::chrono::steady_clock::now()-budget.start).count();lp.setRealParam(SoPlex::TIMELIMIT,std::max(0.001,budget.seconds-used));}
        Vec obj=objective?offdiag_constraint(*objective):Vec(k);
        DSVectorRational empty(0);
        for(int j=0;j<k;++j)lp.addColRational(LPColRational(to_solver(obj[j]),empty,to_solver(frac(1,2)),to_solver(frac(-1,2))));
        if(margin)lp.addColRational(LPColRational(1,empty,1,0));
        Span affine;
        for(const auto& v:equalities){auto a=offdiag_constraint(v);auto rhs=constraint_rhs(v);DSVectorRational row(k);
            for(int j=0;j<k;++j)if(a[j]!=0)row.add(j,to_solver(a[j]));
            lp.addRowRational(LPRowRational(to_solver(rhs),row,to_solver(rhs)));affine.add(augmented_constraint(v));
        }
        for(const auto& v:cuts){auto a=offdiag_constraint(v);DSVectorRational row(k+1);
            for(int j=0;j<k;++j)if(a[j]!=0)row.add(j,to_solver(a[j]));
            if(margin&&!affine.contains(augmented_constraint(v)))row.add(k,-1);
            lp.addRowRational(LPRowRational(to_solver(constraint_rhs(v)),row,infinity));
        }
        auto status=lp.optimize();
        if(budget.trace)std::cerr<<"LP "<<budget.lp<<" returned status "<<int(status)<<'\n';
        if(status==SPxSolver::INFEASIBLE)return {};
        if(status!=SPxSolver::OPTIMAL)throw Limit("exact LP did not reach an optimum, status "+std::to_string(int(status)));
        DVectorRational primal(k+int(margin));
        if(!lp.getPrimalRational(primal))throw std::runtime_error("missing rational LP primal");
        Vec y;for(int j=0;j<k;++j)y.push_back(from_solver(primal[j]));
        // Independent exact check of the returned primal, before trusting it.
        for(auto& q:y)if(absq(q)>frac(1,2))throw std::runtime_error("LP bound violation");
        auto g=gram(frame.n,y);
        for(auto& v:equalities)if(value(g,v)!=1)throw std::runtime_error("LP equality violation");
        for(auto& v:cuts)if(value(g,v)<1){std::ostringstream message;message<<"LP inequality violation: value="<<value(g,v)<<", vector=";for(auto& x:v)message<<x<<',';message<<", y=";for(auto& x:y)message<<x<<',';throw std::runtime_error(message.str());}
        return y;
    }

} // namespace vor
