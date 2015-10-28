//
// Created by H. Shinaoka on 2015/10/28.
//

#ifndef IMPSOLVER_FASTUPDATE_FORMULA_H
#define IMPSOLVER_FASTUPDATE_FORMULA_H

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

//Implementing Appendix B.1.1 of Luitz's thesis
template<class T>
T
compute_det_ratio_up(
        const alps::numeric::matrix<T> &B, const alps::numeric::matrix<T> &C, const alps::numeric::matrix<T> &D,
        const alps::numeric::matrix<T> &invA) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const size_t N = num_rows(invA);
    const size_t M = num_rows(D);

    assert(num_rows(invA)==num_cols(invA));

    //compute H
    matrix_t C_invA(M,N,0.0), C_invA_B(M,M,0.0);
    gemm(C, invA, C_invA);
    gemm(C_invA, B, C_invA_B);

    return determinant(D-C_invA_B);
}

template<class T>
T
compute_inverse_matrix_up(
        const alps::numeric::matrix<T> &B, const alps::numeric::matrix<T> &C, const alps::numeric::matrix<T> &D,
        const alps::numeric::matrix<T> &invA,
        alps::numeric::matrix<T> &E, alps::numeric::matrix<T> &F, alps::numeric::matrix<T> &G,
        alps::numeric::matrix<T> &H) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const size_t N = num_rows(invA);
    const size_t M = num_rows(D);

    assert(num_rows(invA)==num_cols(invA));

    //fill E, F, G by zero for safety
    std::fill(E.get_values().begin(), E.get_values().end(),0);
    std::fill(F.get_values().begin(), F.get_values().end(),0);
    std::fill(G.get_values().begin(), G.get_values().end(),0);

    //compute H
    matrix_t C_invA(M,N,0.0), C_invA_B(M,M,0.0);
    gemm(C, invA, C_invA);
    gemm(C_invA, B, C_invA_B);
    H = inverse(D-C_invA_B);

    //compute G
    gemm(H,C_invA,G);
    G *= -1.;

    //compute F
    matrix_t invA_B(N, M, 0);
    gemm(invA, B, invA_B);
    gemm(invA_B, H, F);
    F *= -1.0;

    //compute E
    gemm(invA_B, G, E);
    E *= -1;
    E += invA;

    return 1./determinant(H);
}

//Note: invA and invBigMat can point to the same matrix object.
// invBigMat is resized automatically.
template<class T>
T
compute_inverse_matrix_up2(
        const alps::numeric::matrix<T> &B, const alps::numeric::matrix<T> &C, const alps::numeric::matrix<T> &D,
        const alps::numeric::matrix<T> &invA,
        alps::numeric::matrix<T> &invBigMat) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const size_t N = num_rows(invA);
    const size_t M = num_rows(D);

    assert(num_rows(invA)==num_cols(invA));

    matrix_t E(N,N,0), F(N,M,0), G(M,N,0), H(M,M,0);

    //compute H
    matrix_t C_invA(M,N,0.0), C_invA_B(M,M,0.0);
    gemm(C, invA, C_invA);
    gemm(C_invA, B, C_invA_B);
    H = inverse(D-C_invA_B);

    //compute G
    gemm(H,C_invA,G);
    G *= -1.;

    //compute F
    matrix_t invA_B(N, M, 0);
    gemm(invA, B, invA_B);
    gemm(invA_B, H, F);
    F *= -1.0;

    //compute E
    gemm(invA_B, G, E);
    E *= -1;
    E += invA;

    resize(invBigMat, N+M, N+M);
    copy_block(E,0,0,invBigMat,0,0,N,N);
    copy_block(F,0,0,invBigMat,0,N,N,M);
    copy_block(G,0,0,invBigMat,N,0,M,N);
    copy_block(H,0,0,invBigMat,N,N,M,M);

    return 1./determinant(H);
}

#endif //IMPSOLVER_FASTUPDATE_FORMULA_H
