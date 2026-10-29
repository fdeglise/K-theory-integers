#pragma once
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace vor {
using Z = mpz_class;
using Q = mpq_class;
using Vec = std::vector<Q>;
using Mat = std::vector<Vec>;
using IVec = std::vector<long>;
using Config = std::vector<IVec>;

struct Limit : std::runtime_error { using std::runtime_error::runtime_error; };
struct Budget {
    long lp = 0, svp_nodes = 0, families = 0, states = 0;
    long max_lp = 2000000, max_svp_nodes = 20000000, max_families = 1000000;
    double seconds = 0;
    bool trace=false;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    void check() const {
        if (seconds && std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() > seconds)
            throw Limit("time limit reached");
    }
};
inline long integer(const Z& z) {
    if (!z.fits_slong_p() || z==std::numeric_limits<long>::min()) throw Limit("integer coordinate exceeds machine-long representation");
    return z.get_si();
}
inline long integer(const Q& q) {
    if (q.get_den()!=1) throw std::runtime_error("expected an integer");
    return integer(q.get_num());
}
inline Q frac(long a, long b) { Q q(a,b); q.canonicalize(); return q; }
inline Z absz(const Z& z) { return z<0 ? -z : z; }
inline Q absq(const Q& q) { return q<0 ? -q : q; }
inline long floorq(const Q& q) { Z z; mpz_fdiv_q(z.get_mpz_t(),q.get_num_mpz_t(),q.get_den_mpz_t()); return integer(z); }
inline long ceilq(const Q& q) { return -floorq(-q); }
inline long sqrt_ceil(const Q& q) {
    if (q<0) throw std::runtime_error("negative square-root bound");
    Z z, r; mpz_fdiv_q(z.get_mpz_t(),q.get_num_mpz_t(),q.get_den_mpz_t());
    mpz_sqrt(r.get_mpz_t(),z.get_mpz_t());
    if (Q(r*r)<q) ++r;
    return integer(r);
}
inline Vec rational(const IVec& v) { Vec r; for (long x:v) r.emplace_back(x); return r; }
inline Vec mv(const Mat& a,const Vec& x) {
    Vec y(a.size()); for(size_t i=0;i<a.size();++i) for(size_t j=0;j<x.size();++j) y[i]+=a[i][j]*x[j]; return y;
}
inline Q dot(const Vec& a,const Vec& b) { Q r; for(size_t i=0;i<a.size();++i) r+=a[i]*b[i]; return r; }
inline Q value(const Mat& g,const Vec& v) { return dot(v,mv(g,v)); }
inline Mat transpose(const Mat& a) {
    if(a.empty()) return {}; Mat t(a[0].size(),Vec(a.size()));
    for(size_t i=0;i<a.size();++i) for(size_t j=0;j<a[i].size();++j) t[j][i]=a[i][j]; return t;
}
inline Mat multiply(const Mat& a,const Mat& b) {
    Mat c(a.size(),Vec(b[0].size()));
    for(size_t i=0;i<a.size();++i) for(size_t k=0;k<b.size();++k) if(a[i][k]!=0)
        for(size_t j=0;j<b[0].size();++j) c[i][j]+=a[i][k]*b[k][j];
    return c;
}
inline Mat eye(int n) { Mat a(n,Vec(n)); for(int i=0;i<n;++i) a[i][i]=1; return a; }
inline Mat columns(const Config& c) { Mat a; for(const auto& v:c) a.push_back(rational(v)); return transpose(a); }
inline Q det(Mat a) {
    Q d=1; const int n=a.size();
    for(int j=0;j<n;++j) {
        int p=j; while(p<n && a[p][j]==0) ++p; if(p==n) return 0;
        if(p!=j) {std::swap(a[p],a[j]); d=-d;}
        Q pivot=a[j][j]; d*=pivot;
        for(int i=j+1;i<n;++i) {Q t=a[i][j]/pivot; for(int k=j+1;k<n;++k) a[i][k]-=t*a[j][k];}
    } return d;
}
inline Mat inverse(Mat a) {
    int n=a.size(); Mat b=eye(n);
    for(int j=0;j<n;++j) {
        int p=j; while(p<n && a[p][j]==0) ++p;
        if(p==n) throw std::runtime_error("singular matrix");
        std::swap(a[p],a[j]); std::swap(b[p],b[j]); Q t=a[j][j];
        for(int k=0;k<n;++k) {a[j][k]/=t; b[j][k]/=t;}
        for(int i=0;i<n;++i) if(i!=j) {t=a[i][j]; for(int k=0;k<n;++k) {a[i][k]-=t*a[j][k]; b[i][k]-=t*b[j][k];}}
    } return b;
}
struct Span {
    Mat rows; std::vector<int> pivots;
    bool add(Vec v) {
        for(size_t i=0;i<rows.size();++i) {Q t=v[pivots[i]]; if(t!=0) for(size_t j=0;j<v.size();++j) v[j]-=t*rows[i][j];}
        int p=0; while(p<int(v.size()) && v[p]==0) ++p;
        if(p==int(v.size())) return false;
        Q t=v[p]; for(auto& x:v) x/=t;
        // Sorted pivots allow arbitrary subsequent insertion order.
        for(auto& row:rows) {t=row[p]; if(t!=0) for(size_t j=0;j<v.size();++j) row[j]-=t*v[j];}
        auto it=std::lower_bound(pivots.begin(),pivots.end(),p); auto pos=it-pivots.begin();
        pivots.insert(it,p); rows.insert(rows.begin()+pos,std::move(v)); return true;
    }
    bool contains(Vec v) const {
        for(size_t i=0;i<rows.size();++i) {Q t=v[pivots[i]]; if(t!=0) for(size_t j=0;j<v.size();++j) v[j]-=t*rows[i][j];}
        return std::all_of(v.begin(),v.end(),[](const Q& x){return x==0;});
    }
    int rank() const {return rows.size();}
};
inline Vec tensor(const Vec& v) {Vec r; for(size_t i=0;i<v.size();++i) for(size_t j=i;j<v.size();++j) r.push_back(v[i]*v[j]); return r;}
inline int cone_rank(const Config& s) {Span p; for(auto& v:s) p.add(tensor(rational(v))); return p.rank();}
inline IVec normalize(IVec v) {for(long x:v) if(x) {if(x<0) for(auto& y:v) y=-y; break;} return v;}
inline Config normalize(Config s) {for(auto& v:s) v=normalize(v); std::sort(s.begin(),s.end()); s.erase(std::unique(s.begin(),s.end()),s.end()); return s;}
inline Vec normalize(Vec v) {for(const Q& x:v) if(x!=0) {if(x<0) for(auto& y:v) y=-y; break;} return v;}
inline std::string key(const Config& s) {std::ostringstream o; for(auto& v:s) {for(auto x:v) o<<x<<','; o<<';';} return o.str();}
inline void choose(int size,int count,const std::function<void(const std::vector<int>&)>& f) {
    std::vector<int> v;
    std::function<void(int)> go=[&](int first) {if(int(v.size())==count){f(v);return;} for(int j=first;j<=size-count+int(v.size());++j){v.push_back(j);go(j+1);v.pop_back();}};
    go(0);
}
inline Mat gram(int n,const Vec& y) {Mat g=eye(n);int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)g[i][j]=g[j][i]=y[k++];return g;}
inline Vec offdiag_constraint(const Vec& v) {Vec a;for(size_t i=0;i<v.size();++i)for(size_t j=i+1;j<v.size();++j)a.push_back(2*v[i]*v[j]);return a;}
inline Q constraint_rhs(const Vec& v) {Q b=1;for(auto& x:v)b-=x*x;return b;}
inline Vec augmented_constraint(const Vec& v) {auto a=offdiag_constraint(v);a.push_back(constraint_rhs(v));return a;}
} // namespace vor
