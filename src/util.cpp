//
// Created by H. Shinaoka on 2015/10/15.
//

#include <vector>
#include "util.h"

/*
 */


double permutation(size_t N, size_t k) {
    assert(k>0);
    double r=1.0;
    for(size_t i=N-k+1; i<=N; ++i) {
        r *= static_cast<double>(i);
    }
    return r;
}


