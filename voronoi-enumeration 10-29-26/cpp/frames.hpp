#pragma once
#include "exact.hpp"

namespace vor {
// Finite quotient codes: coordinate permutations, signs, and quotient-group
// automorphisms are removed before any expensive lattice realizability calls.
struct FiniteGroup {
    std::vector<int> factors;
    std::vector<std::vector<int>> elements;
    int order=1, exponent=1;
    explicit FiniteGroup(std::vector<int> d):factors(std::move(d)) {
        for(int x:factors) order*=x;
        if(!factors.empty()) exponent=factors.back();
        for(int a=0;a<order;++a) {int t=a;std::vector<int> x;for(int d:factors){x.push_back(t%d);t/=d;}elements.push_back(x);}
    }
    int encode(const std::vector<int>& v) const {int r=0,m=1;for(size_t j=0;j<factors.size();++j){r+=m*((v[j]%factors[j]+factors[j])%factors[j]);m*=factors[j];}return r;}
    int add(int a,int b) const {auto v=elements[a];for(size_t j=0;j<v.size();++j)v[j]+=elements[b][j];return encode(v);}
    int times(int a,int k) const {auto v=elements[a];for(auto& x:v)x*=k;return encode(v);}
    int pair(int a,int b) const {int t=0;for(size_t j=0;j<factors.size();++j)t+=elements[a][j]*elements[b][j]*(exponent/factors[j]);return t%exponent;}
    int sign_class(int a) const {return std::min(a,times(a,-1));}
    bool generates(const std::vector<int>& gs) const {
        std::set<int> found{0};std::vector<int> pending{0};
        for(size_t i=0;i<pending.size();++i) for(int g:gs){int x=add(pending[i],g);if(found.insert(x).second)pending.push_back(x);}
        return int(found.size())==order;
    }
    std::vector<std::vector<int>> automorphisms() const {
        std::vector<std::vector<int>> result; std::vector<int> images;
        std::function<void()> go=[&](){
            int j=images.size();
            if(j==int(factors.size())) {
                if(!generates(images))return;
                std::vector<int> map(order);
                for(int a=0;a<order;++a)for(size_t k=0;k<images.size();++k)map[a]=add(map[a],times(images[k],elements[a][k]));
                result.push_back(std::move(map));return;
            }
            for(int x=0;x<order;++x)if(times(x,factors[j])==0){images.push_back(x);go();images.pop_back();}
        };go();return result;
    }
};
inline std::vector<std::vector<int>> group_types(int order) {
    std::vector<std::vector<int>> result;std::vector<int> d;
    std::function<void(int,int)> go=[&](int left,int previous){
        if(left==1){result.push_back(d);return;}
        for(int x=2;x<=left;++x)if(left%x==0 && (previous==1 || x%previous==0)){d.push_back(x);go(left/x,x);d.pop_back();}
    };go(order,1);return result;
}
inline Mat row_lattice_basis(std::vector<std::vector<Z>> a,int n) {
    for(int j=0;j<n;++j) {
        int p=j;while(p<int(a.size())&&a[p][j]==0)++p;
        if(p==int(a.size()))throw std::runtime_error("lattice rank defect");
        std::swap(a[p],a[j]);
        for(int i=j+1;i<int(a.size());++i)while(a[i][j]!=0){
            Z q=a[i][j]/a[j][j];for(int k=0;k<n;++k)a[i][k]-=q*a[j][k];
            if(a[i][j]!=0)std::swap(a[i],a[j]);
        }
        if(a[j][j]<0)for(auto& x:a[j])x=-x;
        for(int i=0;i<j;++i){Z q;mpz_fdiv_q(q.get_mpz_t(),a[i][j].get_mpz_t(),a[j][j].get_mpz_t());for(int k=0;k<n;++k)a[i][k]-=q*a[j][k];}
    }
    Mat b(n,Vec(n));for(int i=0;i<n;++i)for(int j=0;j<n;++j)b[i][j]=a[i][j];return b;
}
struct Frame {
    int n,index,exponent;
    std::vector<IVec> cosets;
    Mat basis, b; // L=basis*Z^n, b=basis^-1; b has integral entries.
    Config axes;
    Frame(int dim,int idx,int den,std::vector<IVec> code):n(dim),index(idx),exponent(den),cosets(std::move(code)) {
        std::vector<std::vector<Z>> generators;
        for(int i=0;i<n;++i){std::vector<Z> row(n);row[i]=exponent;generators.push_back(row);}
        for(const auto& x:cosets){std::vector<Z> row;for(long y:x)row.emplace_back(y);generators.push_back(row);}
        basis=transpose(row_lattice_basis(generators,n));
        for(auto& row:basis)for(auto& x:row)x/=exponent;
        b=inverse(basis);
        if(det(b)!=index)throw std::runtime_error("incorrect quotient lattice index");
        for(auto& row:b)for(auto& x:row)if(x.get_den()!=1)throw std::runtime_error("nonintegral frame matrix");
        for(int i=0;i<n;++i){Vec e(n);e[i]=1;axes.push_back(ambient(e));}
        axes=normalize(axes);
    }
    IVec ambient(const Vec& x) const {auto y=mv(b,x);IVec v;for(auto& q:y)v.push_back(integer(q));return v;}
    Vec coordinates(const IVec& x) const {return mv(basis,rational(x));}
    Config cube(Budget& budget) const {
        std::set<IVec> result;
        for(const auto& c:cosets){Vec x(n);
            std::function<void(int)> go=[&](int j){
                if(j==n){auto a=normalize(ambient(x));long gcd=0;for(long t:a)gcd=std::gcd(gcd,std::abs(t));if(gcd==1)result.insert(a);return;}
                for(long t=c[j]-exponent;t<=exponent;t+=exponent){x[j]=frac(t,exponent);go(j+1);}
            };go(0);budget.check();
        }return Config(result.begin(),result.end());
    }
    Config watson_forced() const {
        Config result=axes;
        for(auto a:cosets){long weight=0;for(auto& x:a){if(2*x>exponent)x-=exponent;weight+=std::abs(x);}
            if(weight!=2*exponent)continue;
            Vec x;for(long t:a)x.push_back(frac(t,exponent));
            for(int i=0;i<n;++i)if(a[i]){Vec w=x;for(auto& z:w)z*=a[i]>0?-1:1;w[i]+=1;result.push_back(ambient(w));}
            // A binary coset of weight four has average squared length one
            // over all coordinate signs. Every sign choice must be minimal.
            bool binary=std::all_of(a.begin(),a.end(),[&](long t){return t==0||2*std::abs(t)==exponent;});
            if(binary){Vec w(n);std::function<void(int)> signs=[&](int j){
                if(j==n){result.push_back(ambient(w));return;}
                if(a[j]==0){w[j]=0;signs(j+1);}else for(int s:{-1,1}){w[j]=frac(s,2);signs(j+1);}
            };signs(0);}
        }return normalize(result);
    }
};
inline int certified_index_bound(int n) {
    static const int bounds[]={0,1,1,1,2,2,4,8};
    if(n<1||n>7)throw std::runtime_error("certified global enumeration currently supports ranks 1--7");
    return bounds[n];
}
inline std::vector<Frame> frame_catalogue(int n,int bound,Budget& budget) {
    std::vector<Frame> result;
    for(int index=1;index<=bound;++index)for(const auto& type:group_types(index)) {
        FiniteGroup g(type);auto aut=g.automorphisms();
        std::vector<int> classes;for(int a=0;a<g.order;++a)if(g.sign_class(a)==a)classes.push_back(a);
        std::vector<int> code;
        std::function<void(int)> go=[&](int first){
            if(int(code.size())!=n){for(int j=first;j<int(classes.size());++j){code.push_back(classes[j]);go(j);code.pop_back();}return;}
            budget.check();if(!g.generates(code))return;
            std::vector<IVec> cosets;
            for(int x=0;x<g.order;++x){IVec r;long weight=0;for(int a:code){int t=g.pair(x,a);r.push_back(t);weight+=std::min(t,g.exponent-t);}
                if(x && weight<2*g.exponent)return; // necessary Watson inequality
                cosets.push_back(std::move(r));
            }
            for(const auto& map:aut){std::vector<int> image;for(int a:code)image.push_back(g.sign_class(map[a]));std::sort(image.begin(),image.end());if(image<code)return;}
            std::sort(cosets.begin(),cosets.end());
            result.emplace_back(n,index,g.exponent,std::move(cosets));
        };go(0);
    }
    return result;
}
} // namespace vor
