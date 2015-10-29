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

    assert(M>0);

    assert(num_rows(invA)==num_cols(invA));
    assert(num_rows(B)==N && num_cols(B)==M);
    assert(num_rows(C)==M && num_cols(C)==N);
    assert(num_rows(D)==M && num_cols(D)==M);

    if (N==0) {
        return determinant(D);
    } else {
        //compute H
        matrix_t C_invA(M,N,0.0), C_invA_B(M,M,0.0);
        gemm(C, invA, C_invA);
        gemm(C_invA, B, C_invA_B);
        return determinant(D-C_invA_B);
    }
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

    assert(M>0);

    assert(num_rows(invA)==num_cols(invA));
    assert(num_rows(B)==N && num_cols(B)==M);
    assert(num_rows(C)==M && num_cols(C)==N);
    assert(num_rows(D)==M && num_cols(D)==M);

    //compute H
    if (N==0) {
        H = inverse(D);
        E.resize(0,0);
        F.resize(0,M);
        G.resize(M,0);

        return 1./determinant(H);
    } else {
        E.resize(N,N);
        F.resize(N,M);
        G.resize(M,N);

        //fill E, F, G by zero for safety
        std::fill(E.get_values().begin(), E.get_values().end(),0);
        std::fill(F.get_values().begin(), F.get_values().end(),0);
        std::fill(G.get_values().begin(), G.get_values().end(),0);

        matrix_t C_invA(M, N, 0.0), C_invA_B(M, M, 0.0);
        gemm(C, invA, C_invA);
        gemm(C_invA, B, C_invA_B);
        H = inverse(D - C_invA_B);

        //compute G
        gemm(H, C_invA, G);
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

    assert(M>0);

    assert(num_rows(invA)==num_cols(invA));
    assert(num_rows(B)==N && num_cols(B)==M);
    assert(num_rows(C)==M && num_cols(C)==N);
    assert(num_rows(D)==M && num_cols(D)==M);

    if (N==0) {
        invBigMat = inverse(D);
        return determinant(D);
    } else {
        matrix_t E(N, N, 0), F(N, M, 0), G(M, N, 0), H(M, M, 0);

        //compute H
        matrix_t C_invA(M, N, 0.0), C_invA_B(M, M, 0.0);
        gemm(C, invA, C_invA);
        gemm(C_invA, B, C_invA_B);
        H = inverse(D - C_invA_B);

        //compute G
        gemm(H, C_invA, G);
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

        resize(invBigMat, N + M, N + M);
        copy_block(E, 0, 0, invBigMat, 0, 0, N, N);
        copy_block(F, 0, 0, invBigMat, 0, N, N, M);
        copy_block(G, 0, 0, invBigMat, N, 0, M, N);
        copy_block(H, 0, 0, invBigMat, N, N, M, M);

        return 1. / determinant(H);
    }
}

template<class T>
T
compute_det_ratio_down(
        const size_t num_rows_cols_removed,
        const std::vector<size_t>& rows_removed,
        const std::vector<size_t>& cols_removed,
        const alps::numeric::matrix<T>& invBigMat) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const size_t NpM = num_rows(invBigMat);
    const size_t M = num_rows_cols_removed;
    assert(num_cols(invBigMat)==NpM);
    assert(rows_removed.size()>=M);
    assert(cols_removed.size()>=M);

    matrix_t H(M,M);
    for (size_t j=0; j<M; ++j) {
        for (size_t i=0; i<M; ++i) {
            H(i,j) = invBigMat(rows_removed[i], cols_removed[j]);
        }
    }
    return determinant(H);
}

//Note: invBigMat will be shrinked.
template<class T>
void
compute_inverse_matrix_down(
    const size_t num_rows_cols_removed,
    const std::vector<size_t>& rows_removed,
    const std::vector<size_t>& cols_removed,
    alps::numeric::matrix<T>& BigMat,
    alps::numeric::matrix<T>& invBigMat
    ) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const size_t NpM = num_rows(invBigMat);
    const size_t M = num_rows_cols_removed;
    const size_t N = NpM-M;
    assert(num_cols(invBigMat)==NpM);
    assert(rows_removed.size()>=M);
    assert(cols_removed.size()>=M);
    assert(M>0);

    if (M==0) {
        throw std::logic_error("M should be larger than 0!");
    }

#ifndef NDEBUG
    //make sure the indices are in ascending order.
    for (size_t i=0; i<M-1; ++i) {
        assert(cols_removed[i]<cols_removed[i+1]);
        assert(rows_removed[i]<rows_removed[i+1]);
    }
#endif

    //move rows and cols to be removed to the end.
    for (size_t i=0; i<M; ++i) {
        if(cols_removed[M-1-i]!=NpM-1-i) {
           BigMat.swap_cols(cols_removed[M-1-i], NpM-1-i);
           invBigMat.swap_cols(cols_removed[M-1-i], NpM-1-i);
        }
    }
    for (size_t i=0; i<M; ++i) {
        if (rows_removed[M-1-i]!=NpM-1-i) {
            BigMat.swap_rows(rows_removed[M-1-i], NpM-1-i);
            invBigMat.swap_rows(rows_removed[M-1-i], NpM-1-i);
        }
    }

    matrix_t E(N,N), F(N,M), G(M,N), H(M,M);
    copy_block(invBigMat,0,0,E,0,0,N,N);
    copy_block(invBigMat,0,N,F,0,0,N,M);
    copy_block(invBigMat,N,0,G,0,0,M,N);
    copy_block(invBigMat,N,N,H,0,0,M,M);

    matrix_t invH_G(M,N,0), F_invH_G(N,N,0);//one might reuse memories...
    gemm(inverse(H), G, invH_G);
    gemm(F, invH_G, F_invH_G);

    invBigMat.resize(N,N);
    invBigMat = E-F_invH_G;
}
#endif //IMPSOLVER_FASTUPDATE_FORMULA_H
