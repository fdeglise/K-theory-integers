#pragma once
#include "frames.hpp"

namespace vor {
using Cube = std::map<IVec,IVec>; // frame-coordinate numerator -> ambient ray
inline IVec scaled_coordinates(const Frame& frame,const IVec& v) {
    IVec a;for(auto& x:frame.coordinates(v))a.push_back(integer(frame.exponent*x));return normalize(a);
}
inline Cube coordinate_cube(const Frame& frame,const Config& pool) {
    Cube result;for(auto& v:pool)result.emplace(scaled_coordinates(frame,v),v);return result;
}
inline Config elementary_orbits(const Frame& frame,const Config& s,const Cube& cube,Budget& budget) {
    // A deliberately small, certified subgroup of the cell stabilizer. Each
    // coordinate flip/swap must preserve both the quotient code and all rays.
    // Orbit BFS is over the finite candidate cube, never over a huge group.
    struct Generator {int i,j;};
    std::vector<Generator> generators;
    std::set<IVec> code(frame.cosets.begin(),frame.cosets.end()),minimal;
    for(auto& v:s)minimal.insert(scaled_coordinates(frame,v));
    auto transform=[](IVec v,Generator g){if(g.j<0)v[g.i]=-v[g.i];else std::swap(v[g.i],v[g.j]);return v;};
    auto consider=[&](Generator g){
        for(auto c:frame.cosets){c=transform(c,g);for(auto& x:c)x=(x%frame.exponent+frame.exponent)%frame.exponent;if(!code.count(c))return;}
        for(auto v:minimal)if(!minimal.count(normalize(transform(v,g))))return;
        generators.push_back(g);
    };
    for(int i=0;i<frame.n;++i){consider({i,-1});for(int j=i+1;j<frame.n;++j)consider({i,j});}
    Config representatives;std::set<IVec> seen;
    for(auto& [a,v]:cube){
        if(!seen.insert(a).second)continue;
        representatives.push_back(v);std::vector<IVec> pending{a};
        for(size_t i=0;i<pending.size();++i)for(auto g:generators){
            IVec image=normalize(transform(pending[i],g));
            if(!cube.count(image))throw std::runtime_error("certified symmetry escaped the candidate cube");
            if(seen.insert(image).second)pending.push_back(std::move(image));
        }
        budget.check();
    }
    return representatives;
}
} // namespace vor
