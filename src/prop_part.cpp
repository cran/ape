/* prop_part.cpp    2017-07-26 */

/* Copyright 2017 Klaus Schliep */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::vector<int> > bipartition2(IntegerMatrix orig, int nTips) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent), j=0;
    int nnode = m - nTips;
    // create list for results
    std::vector< std::vector<int> > out(nnode);
    std::vector<int> y;
    for(int i = 0; i<parent.size(); i++){
        j = parent[i] - nTips - 1L;
        if(children[i] > nTips){
            y = out[children[i] - nTips -1L];
            out[j].insert( out[j].end(), y.begin(), y.end() );
        }
        else out[j].push_back(children[i]);
    }
    for(int i=0; i<nnode; ++i){
        sort(out[i].begin(), out[i].end());
    }
    return out;
}


// [[Rcpp::export]]
List prop_part2(SEXP trees, int nTips){
    List tr(trees);
    int nbtree = tr.size();//, KeepPartition=1; // unused (EP 2020-05-02)
    List M = tr(0);
    IntegerMatrix E = M["edge"];
    std::vector< std::vector<int> > ans = bipartition2(E, nTips);
    std::vector<int> no;
    for(unsigned int i=0; i<ans.size();++i) no.push_back(1);
    no[0] = nbtree;
    for(int k=1; k<nbtree; ++k){
        List tmpTree = tr(k);
        IntegerMatrix tmpE = tmpTree["edge"];
        std::vector< std::vector<int> > bp = bipartition2(tmpE, nTips);
        for (unsigned int i = 1; i < bp.size(); i++) {
            unsigned int j = 1;
            next_j:
                if (bp[i] == ans[j]) {
                    no[j]++;
                    continue;
                }
                j++;
                if (j < ans.size()) goto next_j;
                else  {   //if(KeepPartition)
                    ans.push_back(bp[i]);
                    no.push_back(1);
                }
        }
    }
    List output = wrap(ans);
    output.attr("number") = no;
    output.attr("class") = "prop.part";
    return output;
}
