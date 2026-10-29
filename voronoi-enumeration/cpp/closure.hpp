#pragma once
#include "frames.hpp"
#include <optional>
#include <iostream>

namespace vor {
struct LDL {
    bool positive=true;
    Vec diagonal, bad;
    Mat upper;
};
inline LDL decompose(const Mat& g) {
    int n=g.size();LDL l;l.diagonal.resize(n);l.upper=eye(n);
    for(int i=0;i<n;++i){
        Q d=g[i][i];for(int k=0;k<i;++k)d-=l.diagonal[k]*l.upper[k][i]*l.upper[k][i];
        if(d<=0){
            l.positive=false;l.bad=Vec(n);l.bad[i]=1;
            for(int k=i-1;k>=0;--k){l.bad[k]=-l.upper[k][i];for(int j=k+1;j<i;++j)l.bad[k]-=l.upper[k][j]*l.bad[j];}
            Z denominator=1;for(auto& x:l.bad)mpz_lcm(denominator.get_mpz_t(),denominator.get_mpz_t(),x.get_den_mpz_t());
            for(auto& x:l.bad)x*=denominator;
            if(value(g,l.bad)>0)throw std::runtime_error("invalid nonpositive-direction witness");
            return l;
        }
        l.diagonal[i]=d;
        for(int j=i+1;j<n;++j){Q x=g[i][j];for(int k=0;k<i;++k)x-=l.diagonal[k]*l.upper[k][i]*l.upper[k][j];l.upper[i][j]=x/d;}
    }return l;
}
inline std::vector<Vec> short_vectors(const Mat& g,const Frame& frame,Budget& budget,bool stop_on_violation=false) {
    LDL l=decompose(g);if(!l.positive)throw std::runtime_error("SVP requires positive definiteness");
    std::set<Vec> found;Vec x(frame.n);long visited=0,trials=0;bool stop=false;
    for(const auto& c:frame.cosets){
        std::function<void(int,Q)> go=[&](int i,Q remaining){
            ++budget.svp_nodes;if(++visited>budget.max_svp_nodes)throw Limit("short-vector node limit reached");
            if((visited&1023)==0)budget.check();
            if(i<0){if(std::any_of(x.begin(),x.end(),[](const Q& q){return q!=0;})){found.insert(normalize(x));if(stop_on_violation&&remaining>0)stop=true;}return;}
            Q center=0;for(int j=i+1;j<frame.n;++j)center-=l.upper[i][j]*x[j];
            Q offset=frac(c[i],frame.exponent), shift=center-offset;
            long radius=sqrt_ceil(remaining/l.diagonal[i]);
            long lo=ceilq(shift-radius),hi=floorq(shift+radius);
            long nearest=floorq(shift+frac(1,2));
            // Nearest-first enumeration avoids walking along a very long,
            // nearly null direction before finding a small violated vector.
            for(long step=0;!stop;++step){
                bool in_range=false;
                for(int sign:{1,-1}){if(step==0&&sign<0)continue;Z zz=Z(nearest)+Z(sign)*step;if(zz<lo||zz>hi)continue;in_range=true;
                    if(++trials>budget.max_svp_nodes)throw Limit("short-vector interval work limit reached");
                    long z=integer(zz);x[i]=z+offset;Q term=l.diagonal[i]*(x[i]-center)*(x[i]-center);if(term<=remaining)go(i-1,remaining-term);if(stop)break;
                }
                if(!in_range&&Z(nearest)+step>hi&&Z(nearest)-step<lo)break;
            }
        };go(frame.n-1,Q(1));
        if(stop)break;
    }return std::vector<Vec>(found.begin(),found.end());
}
inline Vec small_nonpositive_cut(const Mat& g,const Vec& direction) {
    // Clearing denominators in a nearly singular LP point can create enormous
    // coefficients. First round the same direction to small integer vectors;
    // a cut only needs G[v]<1, not G[v]<=0. Every proposed cut is checked exactly.
    Q largest=0;for(auto& x:direction)largest=std::max(largest,absq(x));
    for(long scale=1;scale<=64;++scale){Vec v;bool nonzero=false;
        for(auto& x:direction){long t=floorq(scale*x/largest+frac(1,2));v.emplace_back(t);nonzero|=t!=0;}
        if(nonzero&&value(g,v)<1)return normalize(v);
    }
    return normalize(direction); // decompose returns an integral exact witness
}
struct Witness {Mat g;std::vector<Vec> minimal;};
class ClosureOracle {
    const Frame& frame;
    Budget& budget;
    std::set<Vec> cuts;

    std::optional<Vec> linear_program(const std::vector<Vec>& equalities,const Vec* objective,bool margin);

    std::optional<Witness> optimize(const std::vector<Vec>& equalities,const Vec* objective=nullptr,bool margin=false) {
        for(;;){
            auto y=linear_program(equalities,objective,margin);if(!y)return {};
            Mat g=gram(frame.n,*y);auto l=decompose(g);
            if(!l.positive){if(!cuts.insert(small_nonpositive_cut(g,l.bad)).second)throw std::runtime_error("repeated nonpositive cut");continue;}
            if(budget.trace)std::cerr<<"SVP after LP "<<budget.lp<<" start\n";
            auto minimal=short_vectors(g,frame,budget,true);bool violated=false;
            if(budget.trace)std::cerr<<"SVP after LP "<<budget.lp<<" done, "<<minimal.size()<<" vectors\n";
            for(auto& v:minimal)if(value(g,v)<1){violated=true;if(!cuts.insert(v).second)throw std::runtime_error("repeated short-vector cut");}
            if(!violated)return Witness{std::move(g),std::move(minimal)};
        }
    }
public:
    ClosureOracle(const Frame& f,Budget& b):frame(f),budget(b){
        for(int i=0;i<frame.n;++i)for(int j=i+1;j<frame.n;++j)for(int s:{-1,1}){Vec v(frame.n);v[i]=1;v[j]=s;cuts.insert(normalize(v));}
        for(auto& c:frame.cosets){Vec v;for(long x:c)v.push_back(frac(2*x>frame.exponent?x-frame.exponent:x,frame.exponent));if(std::any_of(v.begin(),v.end(),[](const Q& q){return q!=0;}))cuts.insert(normalize(v));}
    }
    std::optional<Witness> close(const Config& family,int cutoff) {
        std::vector<Vec> forced;Span tensor_span,affine;
        for(auto& v:family){auto x=frame.coordinates(v);forced.push_back(x);tensor_span.add(tensor(x));affine.add(augmented_constraint(x));}
        if(tensor_span.rank()>cutoff)return {};
        auto result=optimize(forced,nullptr,true);if(!result)return {};
        // Remove accidental minima by averaging feasible optimizers. A minimum
        // survives this averaging exactly when it is minimal at both endpoints.
        for(;;){bool changed=false;
            for(auto v:result->minimal){
                if(affine.contains(augmented_constraint(v)))continue;
                auto other=optimize(forced,&v,false);
                if(!other)throw std::runtime_error("a known feasible face became infeasible");
                Q maximum=value(other->g,v);
                if(maximum==1){forced.push_back(v);affine.add(augmented_constraint(v));tensor_span.add(tensor(v));if(tensor_span.rank()>cutoff)return {};}
                else if(maximum>1){
                    for(int i=0;i<frame.n;++i)for(int j=0;j<frame.n;++j)result->g[i][j]=(result->g[i][j]+other->g[i][j])/2;
                    std::vector<Vec> intersection;std::set_intersection(result->minimal.begin(),result->minimal.end(),other->minimal.begin(),other->minimal.end(),std::back_inserter(intersection));result->minimal=std::move(intersection);changed=true;break;
                }else throw std::runtime_error("invalid upper bound for a known minimum");
            }
            if(!changed)break;
        }
        return result;
    }
};
} // namespace vor
