#include <algorithm>

#include "gtest.h"
#include "common.hpp"

#include "../src/submatrix.hpp"

TEST(FastUpdate, BlockMatrixReplaceRowsColsSingular) {
    typedef alps::numeric::matrix<double> matrix_t;

    for (int i=0; i<20; ++i) {
        double alpha = pow(0.1, 20*i);
        double a=4., b=8.;

        matrix_t G2(2,2), G2p(2,2), M2p(2,2);
        std::vector<std::pair<int, int> > swap_list;

        G2(0,0) = alpha;
        G2(1,1) = alpha;
        G2(1,0) = a;
        G2(0,1) = a;

        M2p(0,0) = alpha;
        M2p(1,0) = -b;
        M2p(0,1) = -b;
        M2p(1,1) = alpha;
        M2p.block() /= (alpha+b)*(alpha-b);

        matrix_t Q(1,1), R(1,1), S(1,1);
        Q(0,0) = b;
        R(0,0) = b;
        S(0,0) = alpha;

        matrix_t M2p_fast = inverse(G2);
        matrix_t Mmat, inv_tSp;
        matrix_t tPp, tQp, tRp, tSp;

        double det_rat_fast = compute_det_ratio_replace_rows_cols_safe(M2p_fast, Q, R, S, Mmat, inv_tSp);
        double det_rat = (alpha-b)*(alpha+b)/((alpha-a)*(alpha+a));
    }
}

TEST(FastUpdate, ReplaceDiagonalElements) {
    typedef double T;
    typedef alps::numeric::matrix<T> matrix_t;

    const int N=10, m=2, offset=2;
    //const int N=2, m=1, offset=0;
    assert(m+offset<=N);

    matrix_t A_old(N,N), A_new(N,N), new_elems(m,1);
    std::vector<T> elems_diff(m);
    std::vector<int> pos(m);

    randomize_matrix(A_old, 100);
    randomize_matrix(new_elems, 200);
    A_new = A_old;
    matrix_t invA_old = inverse(A_old);
    for (int i=0; i<m; ++i) {
        pos[i] = i+offset;
    }
    for (int i=0; i<m; ++i) {
        elems_diff[i] = new_elems(i,0)-A_old(pos[i],pos[i]);
        A_new(pos[i],pos[i]) = new_elems(i,0);
    }

    const T det_rat = determinant(A_new)/determinant(A_old);
    const T det_rat_fast = compute_det_ratio_replace_diaognal_elements(invA_old, m, pos, elems_diff, true);
    ASSERT_TRUE(std::abs((det_rat-det_rat_fast)/det_rat)<1E-8);

    /* inverse matrix update */
    matrix_t invA_new = inverse(A_new);
    matrix_t invA_new_fast = invA_old;
    compute_det_ratio_replace_diaognal_elements(invA_new_fast, m, pos, elems_diff, false);

    ASSERT_TRUE(std::abs(alps::fastupdate::norm_square(invA_new-invA_new_fast))<1E-5);
}

