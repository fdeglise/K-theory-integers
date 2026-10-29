#include "closure.hpp"
#include "configurations.hpp"
#include "orbits.hpp"
#include <fstream>
#include <iostream>
#include <optional>

namespace vor {
struct Cell {Config rays;Mat form;int dimension,index;std::optional<Stabilizer> group;};
struct Options {int n=0,excess=-1,index_bound=0;bool groups=false,verbose=false,self_test=false;std::string output;Budget budget;};
struct Run {std::vector<Cell> cells;int frames=0;};

inline Run enumerate(Options& options) {
    int n=options.n,bound=certified_index_bound(n);
    if(options.excess<0||options.excess>n*(n-1)/2)throw std::runtime_error("invalid cone-dimension excess");
    int cutoff=n+options.excess;
    if(options.index_bound){if(options.index_bound<bound)throw std::runtime_error("index bound is below the certified completeness bound");bound=options.index_bound;}
    auto& budget=options.budget;Run run;Registry states;
    // Trial families are searched under a fixed maximal-index restriction.
    // A forced root may already contain a larger-index frame; rejecting it in
    // one index must not suppress the equivalent family at its proper index.
    std::map<int,Registry> families;
    auto catalogue=frame_catalogue(n,bound,budget);run.frames=catalogue.size();
    if(options.verbose)std::cerr<<"Symmetry-reduced frame models surviving Watson: "<<run.frames<<'\n';
    for(size_t fi=0;fi<catalogue.size();++fi){
        auto& frame=catalogue[fi];auto forced=frame.watson_forced();
        if(cone_rank(forced)>cutoff)continue;
        if(options.verbose)std::cerr<<"Frame "<<fi+1<<'/'<<catalogue.size()<<", index "<<frame.index<<", exponent "<<frame.exponent<<'\n';
        ClosureOracle oracle(frame,budget);std::vector<Config> pending;
        std::function<void(Config)> visit=[&](Config family){
            budget.check();if(++budget.families>budget.max_families)throw Limit("family count limit reached");
            if(budget.trace)std::cerr<<"Family "<<budget.families<<": equivalence check start\n";
            if(!families[frame.index].insert(family).second)return;
            if(budget.trace)std::cerr<<"Family "<<budget.families<<": closure start\n";
            auto witness=oracle.close(family,cutoff);if(!witness)return;
            Config minimal;for(auto& v:witness->minimal)minimal.push_back(frame.ambient(v));minimal=normalize(minimal);
            int r=cone_rank(minimal);if(r>cutoff||!maximal_frame(minimal,frame.index))return;
            auto [id,fresh]=states.insert(minimal);if(!fresh)return;
            ++budget.states;Mat form=multiply(transpose(frame.basis),multiply(witness->g,frame.basis));
            Cell cell{minimal,form,r,frame.index,{}};
            if(options.groups)cell.group=stabilizer(*states.entries[id],budget);
            run.cells.push_back(std::move(cell));if(r<cutoff)pending.push_back(minimal);
            if(options.verbose)std::cerr<<"  orbit "<<run.cells.size()<<": r="<<r<<", rays="<<minimal.size()<<", max-index="<<frame.index<<", LP="<<budget.lp<<'\n';
        };
        visit(forced);if(pending.empty())continue;
        auto cube=coordinate_cube(frame,frame.cube(budget));
        while(!pending.empty()){
            Config s=std::move(pending.back());pending.pop_back();
            Span span;for(auto& v:s)span.add(tensor(rational(v)));ExtensionFilter filter(s,frame.index);
            Config candidates=elementary_orbits(frame,s,cube,budget);
            for(auto& v:candidates){
                budget.check();if(std::binary_search(s.begin(),s.end(),v)||!filter.allows(v)||span.contains(tensor(rational(v))))continue;
                Config family=s;family.push_back(v);visit(normalize(std::move(family)));
            }
        }
    }
    std::sort(run.cells.begin(),run.cells.end(),[](const Cell& a,const Cell& b){return std::tie(a.dimension,a.rays)<std::tie(b.dimension,b.rays);});
    return run;
}

inline void write_json(std::ostream& o,const Options& options,const Run& run) {
    o<<"{\n  \"n\": "<<options.n<<",\n  \"complete\": true,\n  \"max_cone_dimension\": "<<options.n+options.excess
     <<",\n  \"index_bound\": "<<(options.index_bound?options.index_bound:certified_index_bound(options.n))
     <<",\n  \"frame_models\": "<<run.frames<<",\n  \"lp_solves\": "<<options.budget.lp<<",\n  \"svp_nodes\": "<<options.budget.svp_nodes<<",\n  \"cells\": [\n";
    for(size_t k=0;k<run.cells.size();++k){auto& c=run.cells[k];if(k)o<<",\n";
        o<<"    {\"cone_dimension\": "<<c.dimension<<", \"projective_dimension\": "<<c.dimension-1<<", \"maximal_frame_index\": "<<c.index<<",\n     \"rays\": [";
        for(size_t i=0;i<c.rays.size();++i){if(i)o<<',';o<<'[';for(size_t j=0;j<c.rays[i].size();++j){if(j)o<<',';o<<c.rays[i][j];}o<<']';}
        o<<"],\n     \"form\": [";
        for(size_t i=0;i<c.form.size();++i){if(i)o<<',';o<<'[';for(size_t j=0;j<c.form[i].size();++j){if(j)o<<',';o<<'"'<<c.form[i][j].get_str()<<'"';}o<<']';}
        o<<"], \"stabilizer_order\": ";if(c.group)o<<c.group->order;else o<<"null";
        o<<", \"orientation_reversing\": ";if(c.group)o<<(c.group->reversing?"true":"false");else o<<"null";o<<'}';
    }o<<"\n  ]\n}\n";
}
inline void require(bool ok,const std::string& what){if(!ok)throw std::runtime_error("self-test failed: "+what);}
inline void self_test() {
    Budget b;auto f2=frame_catalogue(2,1,b);require(f2.size()==1,"rank-two catalogue");
    Mat g={{Q(1),Q(-5)},{Q(-5),Q(26)}};auto shortv=short_vectors(g,f2[0],b);
    require(std::set<Vec>(shortv.begin(),shortv.end())==std::set<Vec>{{Q(1),Q(0)},{Q(5),Q(1)}},"exact ellipsoid enumeration");
    for(Mat bad:{Mat{{1,1},{1,1}},Mat{{1,2},{2,1}}}){auto d=decompose(bad);require(!d.positive&&value(bad,d.bad)<=0,"nonpositive direction");}
    auto f4=frame_catalogue(4,2,b);require(f4.size()==2,"rank-four quotient catalogue");
    for(auto& frame:f4)if(frame.index==2){ClosureOracle oracle(frame,b);auto result=oracle.close(frame.watson_forced(),10);require(result&&result->minimal.size()==12,"forced D4 minima");Config s;for(auto& v:result->minimal)s.push_back(frame.ambient(v));Metric m(normalize(s));require(stabilizer(m,b).order==1152,"D4 stabilizer");}
    bool tested_e7=false;
    for(auto& frame:frame_catalogue(7,8,b))if(frame.index==8&&frame.exponent==2){
        auto forced=frame.watson_forced();require(cone_rank(forced)==28,"noncyclic E7 frame rank");
        ClosureOracle oracle(frame,b);auto result=oracle.close(forced,28);
        require(result&&result->minimal.size()==63,"noncyclic E7 frame closure");tested_e7=true;
    }
    require(tested_e7,"index-eight elementary-abelian quotient is present");
    auto z7=frame_catalogue(7,1,b);auto cube=coordinate_cube(z7[0],z7[0].cube(b));
    require(elementary_orbits(z7[0],z7[0].axes,cube,b).size()==7,"hyperoctahedral cube orbits");
    for(int n=1;n<=3;++n){Options o;o.n=n;o.excess=n*(n+1)/2-n;auto r=enumerate(o);require(r.cells.size()==size_t(n==1?1:n==2?2:5),"small-rank orbit counts");}
    std::cout<<"Exact arithmetic, quotient frames, SVP, LP closure, D4, isometry and small-rank tests passed.\n";
}
inline long number(const std::string& s) {size_t end=0;long v=std::stol(s,&end);if(end!=s.size())throw std::runtime_error("invalid integer: "+s);return v;}
inline int small_number(const std::string& s) {long v=number(s);if(v<std::numeric_limits<int>::min()||v>std::numeric_limits<int>::max())throw std::runtime_error("integer argument is out of range");return int(v);}
inline Options parse(int argc,char** argv) {
    Options o;for(int i=1;i<argc;++i){std::string a=argv[i];auto next=[&](){if(i+1>=argc)throw std::runtime_error("missing value for "+a);return std::string(argv[++i]);};
        if(a=="--self-test")o.self_test=true;
        else if(a=="--excess"){o.excess=small_number(next());if(o.excess<0)throw std::runtime_error("excess must be nonnegative");}
        else if(a=="--index-bound"){o.index_bound=small_number(next());if(o.index_bound<1||o.index_bound>8)throw std::runtime_error("this quotient-code implementation supports index bounds 1--8");}
        else if(a=="--json")o.output=next();
        else if(a=="--stabilizers")o.groups=true;
        else if(a=="--verbose")o.verbose=true;
        else if(a=="--trace")o.budget.trace=true;
        else if(a=="--max-lp")o.budget.max_lp=number(next());
        else if(a=="--max-families")o.budget.max_families=number(next());
        else if(a=="--max-svp-nodes")o.budget.max_svp_nodes=number(next());
        else if(a=="--seconds")o.budget.seconds=number(next());
        else if(a=="--help"){std::cout<<"Usage: voronoi_cpp N [--excess K] [--json FILE] [--stabilizers] [--verbose]\n"
            <<"N=1..7; cone dimensions N..N+K; default K=1.\n"
            <<"Limits: --max-lp N --max-families N --max-svp-nodes N --seconds N.\n"
            <<"--index-bound may increase, but never decrease, the certified bound.\n"
            <<"Run --self-test for exact arithmetic and geometric regression checks.\n";std::exit(0);}
        else if(!a.empty()&&a[0]!='-'&&o.n==0)o.n=small_number(a);
        else throw std::runtime_error("unknown argument: "+a);
    }
    if(!o.self_test){certified_index_bound(o.n);if(o.excess==-1)o.excess=std::min(1,o.n*(o.n+1)/2-o.n);}
    if(o.index_bound<0||o.budget.max_lp<1||o.budget.max_families<1||o.budget.max_svp_nodes<1||o.budget.seconds<0)throw std::runtime_error("invalid resource bound");
    return o;
}
} // namespace vor

int main(int argc,char** argv) {
    try {
        auto options=vor::parse(argc,argv);if(options.self_test){vor::self_test();return 0;}
        auto run=vor::enumerate(options);std::map<int,int> counts;for(auto& c:run.cells)++counts[c.dimension];
        // Only open output after successful exhaustive completion. A failed or
        // capped search leaves any existing result file untouched.
        if(!options.output.empty()){std::ofstream f(options.output);if(!f)throw std::runtime_error("cannot open output");vor::write_json(f,options,run);if(!f)throw std::runtime_error("output write failed");}
        std::cout<<"Complete enumeration for n="<<options.n<<", cone dimensions "<<options.n<<".."<<options.n+options.excess<<"\ncone dimension | cell orbits\n";
        for(auto [r,count]:counts)std::cout<<r<<" | "<<count<<'\n';
        std::cout<<"Frames="<<run.frames<<", LP solves="<<options.budget.lp<<", SVP nodes="<<options.budget.svp_nodes<<"\n";
    } catch(const std::exception& e){std::cerr<<"Enumeration failed; no complete result: "<<e.what()<<'\n';return 2;}
}
