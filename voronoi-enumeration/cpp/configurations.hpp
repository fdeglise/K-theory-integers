#pragma once
#include "exact.hpp"
#include <memory>
#include <unordered_map>

namespace vor {
struct Metric {
    Config vectors;
    Mat dual, base_inverse, basis;
    std::vector<int> base;
    std::vector<std::vector<Z>> pairings;
    std::vector<std::string> labels;
    std::string invariant;
    mutable std::vector<int> exchange_components;
    explicit Metric(Config s):vectors(normalize(std::move(s))) {
        int n=vectors[0].size(),size=vectors.size();Mat c(n,Vec(n));
        for(auto& v:vectors)for(int i=0;i<n;++i)for(int j=0;j<n;++j)c[i][j]+=Q(v[i])*v[j];
        Q determinant=det(c);dual=inverse(c);pairings.resize(size,std::vector<Z>(size));
        std::vector<Vec> d;for(auto& v:vectors)d.push_back(mv(dual,rational(v)));
        for(int i=0;i<size;++i)for(int j=0;j<size;++j){Q q=determinant*dot(rational(vectors[i]),d[j]);if(q.get_den()!=1)throw std::runtime_error("nonintegral intrinsic pairing");pairings[i][j]=q.get_num();}
        for(int i=0;i<size;++i){std::vector<Z> row;for(auto& x:pairings[i])row.push_back(absz(x));std::sort(row.begin(),row.end());std::ostringstream o;o<<pairings[i][i]<<':';for(auto& x:row)o<<x<<',';labels.push_back(o.str());}
        auto sorted=labels;std::sort(sorted.begin(),sorted.end());std::ostringstream o;o<<n<<'/'<<size<<'/'<<determinant<<'/';for(auto& x:sorted)o<<x<<';';invariant=o.str();
        Span span;Config chosen;for(int i=0;i<size&&int(base.size())<n;++i)if(span.add(rational(vectors[i]))){base.push_back(i);chosen.push_back(vectors[i]);}
        if(int(base.size())!=n)throw std::runtime_error("configuration is not full rank");basis=columns(chosen);base_inverse=inverse(basis);
    }
    const std::vector<int>& exchanges() const {
        // Certified source automorphisms permit symmetry breaking when only
        // existence of an isometry is wanted. Do not use it for group orders.
        if(!exchange_components.empty())return exchange_components;
        int n=base.size();exchange_components.resize(n);std::iota(exchange_components.begin(),exchange_components.end(),0);
        auto root=[&](int i){while(exchange_components[i]!=i)i=exchange_components[i];return i;};
        std::set<IVec> wanted(vectors.begin(),vectors.end());
        for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){
            Vec difference(n);for(int k=0;k<n;++k)difference[k]=basis[k][j]-basis[k][i];bool ok=true;
            for(int a=0;a<n&&ok;++a)for(int b=0;b<n;++b){Q entry=difference[a]*(base_inverse[i][b]-base_inverse[j][b]);if(entry.get_den()!=1){ok=false;break;}}
            if(!ok)continue;
            for(auto& v:vectors){Q delta=dot(base_inverse[i],rational(v))-dot(base_inverse[j],rational(v));IVec image;
                for(int k=0;k<n;++k)image.push_back(integer(Q(v[k])+delta*difference[k]));
                if(!wanted.count(normalize(image))){ok=false;break;}
            }
            if(ok)exchange_components[root(i)]=root(j);
        }
        for(int i=0;i<n;++i)exchange_components[i]=root(i);
        return exchange_components;
    }
};
inline bool isometries(const Metric& source,const Metric& target,const std::function<bool(const Mat&)>& callback,bool quotient_source=false) {
    if(source.invariant!=target.invariant)return false;
    int n=source.base.size(),size=target.vectors.size();std::vector<int> chosen,signs;bool found=false,stop=false;
    std::set<IVec> wanted(target.vectors.begin(),target.vectors.end());
    std::vector<int> components;if(quotient_source)components=source.exchanges();
    std::function<void()> go=[&](){
        if(stop)return;int j=chosen.size();
        if(j==n){
            Config image;for(int k=0;k<n;++k){auto v=target.vectors[chosen[k]];for(auto& x:v)x*=signs[k];image.push_back(v);}
            Mat u=multiply(columns(image),source.base_inverse);
            for(auto& row:u)for(auto& x:row)if(x.get_den()!=1)return;
            if(absq(det(u))!=1)return;
            std::set<IVec> actual;
            for(auto& v:source.vectors){IVec w;for(auto& x:mv(u,rational(v)))w.push_back(integer(x));actual.insert(normalize(w));}
            if(actual!=wanted)return;found=true;stop=callback(u);return;
        }
        int src=source.base[j];
        for(int i=0;i<size;++i)if(std::find(chosen.begin(),chosen.end(),i)==chosen.end() && source.labels[src]==target.labels[i])for(int s:{-1,1}){
            if(quotient_source){bool sorted=true;for(int k=0;k<j;++k)if(components[k]==components[j]&&chosen[k]>i){sorted=false;break;}if(!sorted)continue;}
            bool ok=true;for(int k=0;k<j;++k)if(s*signs[k]*target.pairings[i][chosen[k]]!=source.pairings[src][source.base[k]]){ok=false;break;}
            if(ok){chosen.push_back(i);signs.push_back(s);go();chosen.pop_back();signs.pop_back();if(stop)return;}
        }
    };go();return found;
}
struct Registry {
    std::vector<std::shared_ptr<Metric>> entries;
    std::unordered_map<std::string,std::vector<int>> buckets;
    std::unordered_map<std::string,int> raw;
    std::pair<int,bool> insert(const Config& s) {
        auto text=key(s);auto hit=raw.find(text);if(hit!=raw.end())return {hit->second,false};
        auto m=std::make_shared<Metric>(s);
        for(int i:buckets[m->invariant])if(isometries(*m,*entries[i],[](const Mat&){return true;},true)){raw[text]=i;return {i,false};}
        int id=entries.size();buckets[m->invariant].push_back(id);entries.push_back(std::move(m));raw[text]=id;return {id,true};
    }
};
struct Stabilizer {long order=0;bool reversing=false;};
inline Stabilizer stabilizer(const Metric& m,Budget& budget) {
    Span span;Config basis;
    for(auto& v:m.vectors)if(span.add(tensor(rational(v))))basis.push_back(v);
    Mat before;for(auto& v:basis){auto a=tensor(rational(v));Vec row;for(int p:span.pivots)row.push_back(a[p]);before.push_back(row);}
    Q initial=det(before);Stabilizer s;
    isometries(m,m,[&](const Mat& u){
        budget.check();++s.order;Mat after;
        for(auto& v:basis){auto a=tensor(mv(u,rational(v)));Vec row;for(int p:span.pivots)row.push_back(a[p]);after.push_back(row);}
        Q sign=det(after)/initial;if(absq(sign)!=1)throw std::runtime_error("invalid stabilizer orientation");s.reversing=s.reversing||(sign<0);return false;
    });return s;
}
inline bool maximal_frame(const Config& s,int index) {
    bool ok=true;int n=s[0].size();choose(s.size(),n,[&](const std::vector<int>& indices){if(!ok)return;Config rows;for(int i:indices)rows.push_back(s[i]);if(absq(det(columns(rows)))>index)ok=false;});return ok;
}
struct ExtensionFilter {
    std::vector<std::pair<IVec,long>> inequalities;
    explicit ExtensionFilter(const Config& s,int index) {
        int n=s[0].size();std::set<IVec> normals;
        choose(s.size(),n-1,[&](const std::vector<int>& indices){
            IVec normal(n);for(int row=0;row<n;++row){Mat minor;for(int i=0;i<n;++i)if(i!=row){Vec v;for(int j:indices)v.emplace_back(s[j][i]);minor.push_back(v);}normal[row]=integer(det(minor))*(row%2?-1:1);}
            normal=normalize(normal);if(std::any_of(normal.begin(),normal.end(),[](long x){return x!=0;}))normals.insert(normal);
        });
        for(auto v:normals){long gcd=0;for(long x:v)gcd=std::gcd(gcd,std::abs(x));for(auto& x:v)x/=gcd;inequalities.emplace_back(v,index/gcd);}
    }
    bool allows(const IVec& v) const {
        for(auto& [normal,bound]:inequalities){Z q=0;for(size_t j=0;j<v.size();++j)q+=Z(normal[j])*v[j];if(absz(q)>bound)return false;}return true;
    }
};
} // namespace vor
