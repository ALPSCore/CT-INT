//
// Created by H. Shinaoka on 2015/10/28.
//

#ifndef IMPSOLVER_FASTUPDATE_FORMULA_H
#define IMPSOLVER_FASTUPDATE_FORMULA_H

#include <boost/timer/timer.hpp>

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include "util.h"

//Implementing equations in Appendix B.1.1 of Luitz's thesis
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
        return safe_determinant(D);
    } else {
        //compute H
        matrix_t C_invA(M,N,0.0), C_invA_B(M,M,0.0);
        gemm(C, invA, C_invA);
        gemm(C_invA, B, C_invA_B);
        return safe_determinant(D-C_invA_B);
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
        H = safe_inverse(D);
        E.resize(0,0);
        F.resize(0,M);
        G.resize(M,0);

        return 1./safe_determinant(H);
    } else {
        E.resize(N,N);
        F.resize(N,M);
        G.resize(M,N);

        //fill E, F, G by zero for safety
        //std::fill(E.get_values().begin(), E.get_values().end(),0);
        //std::fill(F.get_values().begin(), F.get_values().end(),0);
        //std::fill(G.get_values().begin(), G.get_values().end(),0);

        matrix_t C_invA(M, N, 0.0), C_invA_B(M, M, 0.0);
        gemm(C, invA, C_invA);
        gemm(C_invA, B, C_invA_B);
        H = safe_inverse(D - C_invA_B);

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

        return 1./safe_determinant(H);
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
        invBigMat = safe_inverse(D);
        return safe_determinant(D);
    } else {
        static matrix_t H, C_invA, C_invA_B, invA_B, F;
        H.resize_values_not_retained(M, M);
        C_invA.resize_values_not_retained(M, N);
        C_invA_B.resize_values_not_retained(M, M);
        invA_B.resize_values_not_retained(N, M);
        F.resize_values_not_retained(N,M);

        //compute H
        gemm(C, invA, C_invA);
        gemm(C_invA, B, C_invA_B);
        gemm(C_invA, B, C_invA_B);
        H = safe_inverse(D - C_invA_B);

        //compute F
        gemm(invA, B, invA_B);
        mygemm((T)-1.0, invA_B, H, (T)0.0, F);

        //submatrix views to invBigMat
        invBigMat.resize_values_not_retained(N + M, N + M);//this may destroy the contents of invA as well
        submatrix_view<T> E_view(invBigMat, 0, 0, N, N);
        submatrix_view<T> G_view(invBigMat, N, 0, M, N);

        //compute G
        mygemm((T)-1.0, H, C_invA, (T)0.0, G_view);

        //compute E
        my_copy_block(invA,0,0,E_view,0,0,N,N);
        mygemm(static_cast<T>(-1.0), invA_B, G_view, static_cast<T>(1.0), E_view);

        copy_block(H, 0, 0, invBigMat, N, N, M, M);
        copy_block(F, 0, 0, invBigMat, 0, N, N, M);

        T r = 1. / safe_determinant(H);
        return r;
    }
}

template<class T>
T
compute_det_ratio_down(
        const size_t num_rows_cols_removed,
        const std::vector<size_t>& rows_cols_removed,
        const alps::numeric::matrix<T>& invBigMat) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const size_t NpM = num_rows(invBigMat);
    const size_t M = num_rows_cols_removed;
    assert(num_cols(invBigMat)==NpM);
    assert(rows_cols_removed.size()>=M);
    assert(M>0);

    matrix_t H(M,M);
    for (size_t j=0; j<M; ++j) {
        for (size_t i=0; i<M; ++i) {
            //std::cout << " Debug "  << rows_cols_removed[i] << " " << M << " " << NpM << std::endl;
            assert(rows_cols_removed[i]<NpM);
            assert(rows_cols_removed[j]<NpM);
            H(i,j) = invBigMat(rows_cols_removed[i], rows_cols_removed[j]);
        }
    }
    return safe_determinant(H);
}

//Note: invBigMat will be shrinked.
template<class T>
T
compute_inverse_matrix_down(
    const size_t num_rows_cols_removed,
    const std::vector<size_t>& rows_cols_removed,
    alps::numeric::matrix<T>& invBigMat,
    std::vector<std::pair<size_t,size_t> >& swap_list
    ) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    static matrix_t H, invH_G, F_invH_G;

    const size_t NpM = num_rows(invBigMat);
    const size_t M = num_rows_cols_removed;
    const size_t N = NpM-M;
    assert(num_cols(invBigMat)==NpM);
    assert(rows_cols_removed.size()>=M);
    assert(M>0);
    assert(NpM>=M);

    if (NpM<M) {
        throw std::logic_error("N should not be negative!");
    }

    if (M==0) {
        throw std::logic_error("M should be larger than 0!");
    }

#ifndef NDEBUG
    //make sure the indices are in ascending order.
    for (size_t i=0; i<M-1; ++i) {
        assert(rows_cols_removed[i]<rows_cols_removed[i+1]);
    }
#endif

    //move rows and cols to be removed to the end.
    swap_list.resize(M);
    for (size_t i=0; i<M; ++i) {
        invBigMat.swap_cols(rows_cols_removed[M-1-i], NpM-1-i);
        invBigMat.swap_rows(rows_cols_removed[M-1-i], NpM-1-i);
        swap_list[i] = std::pair<size_t,size_t>(rows_cols_removed[M-1-i], NpM-1-i);
    }

    if (N==0) {
        matrix_t H(invBigMat);
        invBigMat.resize(0,0);
        return safe_determinant(H);
    } else {
        H.resize_values_not_retained(M, M);
        invH_G.resize_values_not_retained(M, N);
        F_invH_G.resize_values_not_retained(N, N);

        submatrix_view<T> E_view(invBigMat, 0, 0, N, N);
        submatrix_view<T> F_view(invBigMat, 0, N, N, M);
        submatrix_view<T> G_view(invBigMat, N, 0, M, N);
        copy_block(invBigMat, N, N, H, 0, 0, M, M);//we need to copy a submatrix to H because save_determinant() does not support a matrix view.

        gemm(safe_inverse(H), G_view, invH_G);
        mygemm((T)-1.0, F_view, invH_G, (T)1.0, E_view);

        invBigMat.resize(N, N);
        return safe_determinant(H);
    }
}

//Implementing Eq. (80) in Otsuki and Kusunose
template<class T>
T
compute_det_ratio_replace_row(const alps::numeric::matrix<T>& M, const alps::numeric::matrix<T> Dr, int m) {
    assert(num_cols(M)==num_rows(M));
    assert(num_rows(Dr)==1);
    assert(num_cols(Dr)==num_rows(M));
    const int N = num_cols(M);

    T lambda = 0;
    for (int i=0; i<N; ++i) {
        lambda += Dr(0,i)*M(i,m);
    }
    return lambda;
}

//Implementing Eq. (81) in Otsuki and Kusunose
template<class T>
T
compute_imverse_matrix_replace_row(alps::numeric::matrix<T>& M, const alps::numeric::matrix<T> Dr, int m) {
    assert(num_cols(M) == num_rows(M));
    assert(num_rows(Dr) == 1);
    assert(num_cols(Dr) == num_rows(M));
    const int N = num_cols(M);

    T lambda = 0;
    for (int i=0; i<N; ++i) {
        lambda += Dr(0,i)*M(i,m);
    }
    const T inv_lambda = 1/lambda;

    alps::numeric::matrix<T> R(1,N,0);
    gemm(Dr,M,R);
    R *= -inv_lambda;
    R(0,m) = -1;

    alps::numeric::matrix<T> new_M(N,N);
    for (int j=0; j<N; ++j) {
        for (int i=0; i<N; ++i) {
            new_M(i,j) = M(i,m)*R(0,j)-M(i,j)*R(0,m);
        }
    }
    for (int i=0; i<N; ++i) {
        new_M(i,m) += M(i,m)*inv_lambda;
    }
    std::swap(M, new_M);
}

//Implementing Eq. (84) in Otsuki and Kusunose
// Not optimized yet.
/*
template<class T>
T
compute_det_ratio_replace_row_col(const alps::numeric::matrix<T>& M, const alps::numeric::matrix<T> Dr, const alps::numeric::matrix<T> Dc, int m) {
    typedef alps::numeric::matrix<T> matrix_t;
    const int N = num_cols(M);
    assert(num_rows(M)==N);
    assert(num_rows(Dr)==1);
    assert(num_cols(Dr)==N);
    assert(num_rows(Dc)==N);
    assert(num_cols(Dc)==1);
    assert(Dr(0,m)==Dc(m,0));

    T lambda = 0;
    {
        T sum_r = 0, sum_c = 0;
        for (int i=0; i<N; ++i) {
            sum_r += Dr(0,i)*M(i,m);
        }
        for (int i=0; i<N; ++i) {
            sum_c += Dc(i,0)*M(m,i);
        }
        lambda += sum_r*sum_c;
    }
    lambda += -M(m,m)*mygemm(mygemm(Dr,M),Dc)(0,0) + Dr(0,m)*M(m,m);

    return lambda;
}
*/

template<class T>
T
compute_inverse_matrix_replace_row_col(alps::numeric::matrix<T>& invG, const alps::numeric::matrix<T> Dr, const alps::numeric::matrix<T> Dc, int m,
    bool assume_intermediate_state_is_singular=false) {
    typedef alps::numeric::matrix<T> matrix_t;
    const double eps = 1E-10;

    const int N = num_cols(invG);
    const int Nm1 = N-1;
    assert(num_rows(invG)==N);
    assert(num_rows(Dr)==1);
    assert(num_cols(Dr)==N);
    assert(num_rows(Dc)==N);
    assert(num_cols(Dc)==1);
    assert(Dr(0,m)==Dc(m,0));

    //original G^{-1}
    matrix_t tP(Nm1, Nm1), tQ(Nm1, 1), tR(1, Nm1);
    invG.swap_cols(m,N-1);
    invG.swap_rows(m,N-1);
    copy_block(invG, 0, 0, tP, 0, 0, Nm1, Nm1);
    copy_block(invG, 0, Nm1, tQ, 0, 0, Nm1, 1);
    copy_block(invG, Nm1, 0, tR, 0, 0, 1, Nm1);
    const T tS = invG(N-1,N-1);

    //new G
    matrix_t Q(Dc), R(Dr);
    const T S = Dc(m,0);
    Q.swap_rows(m,N-1); Q.resize(Nm1,1);
    R.swap_cols(m,N-1); R.resize(1,Nm1);

    //compute lambda
    T lambda = tS*(S-mygemm(R,mygemm(tP,Q))(0,0)) +  mygemm(mygemm(R,tQ),mygemm(tR,Q))(0,0);
    const T tSp = tS/lambda;
    //std::cout << "lamba split " << tS*(S-mygemm(R,mygemm(tP,Q))(0,0)) << " " <<  mygemm(mygemm(R,tQ),mygemm(tR,Q))(0,0) << std::endl;
    //std::cout << "debug tS/tSp " << tS << " " << tSp << std::endl;
    if (assume_intermediate_state_is_singular || (std::abs(tS)<eps && std::abs(tSp)<eps) ) {
        const double R_tQ = mygemm(R, tQ)(0,0);
        const double tR_Q = mygemm(tR, Q)(0,0);
        matrix_t tQp(tQ); tQp /= R_tQ;
        matrix_t tRp(tR); tRp /= tR_Q;

        matrix_t tPp(tP);
        tPp += -mygemm(tQp,mygemm(R,tP))-mygemm(mygemm(tP,Q),tRp)+mygemm(mygemm(tQp,mygemm(R,mygemm(tP,Q))),tRp);

        copy_block(tPp, 0, 0, invG, 0, 0, Nm1, Nm1);
        copy_block(tQp, 0, 0, invG, 0, Nm1, Nm1, 1);
        copy_block(tRp, 0, 0, invG, Nm1, 0, 1, Nm1);
        invG(N-1,N-1) = tSp;

        invG.swap_cols(m,N-1);
        invG.swap_rows(m,N-1);

        return lambda;
    } else {
        matrix_t Mmat = mygemm(tQ,tR);
        Mmat *= -1/tS;
        Mmat += tP;

        const matrix_t MQ = mygemm(Mmat, Q);
        const matrix_t RM = mygemm(R, Mmat);
        const matrix_t tQp = (-tS/lambda)*MQ;
        const matrix_t tRp = (-tS/lambda)*RM;
        const matrix_t tPp = Mmat+(tS/lambda)*mygemm(MQ, RM);

        copy_block(tPp, 0, 0, invG, 0, 0, Nm1, Nm1);
        copy_block(tQp, 0, 0, invG, 0, Nm1, Nm1, 1);
        copy_block(tRp, 0, 0, invG, Nm1, 0, 1, Nm1);
        invG(N-1,N-1) = tSp;

        invG.swap_cols(m,N-1);
        invG.swap_rows(m,N-1);

        return lambda;
    }

/*
    matrix_t L(N,1,0.), R(1,N,0.);
    alps::numeric::gemm(M, Dc, L);
    alps::numeric::gemm(Dr, M, R);

    matrix_t M_new(N,N);
    const double coeff = -1/lambda;
    const double Mmm = M(m,m);
    for(int j=0; j<N; ++j) {
        for(int i=0; i<N; ++i) {
            M_new(i,j) = coeff*(Mmm*L(i,0)-L(m,0)*M(i,m) + Mmm*R(0,j)-M(m,j)*R(0,m));
        }
    }
    M_new(m,m) += Mmm/lambda;
    M.swap(M_new);*/


}

//Assuming the intermediate state is singular, one replaces a row and a column.
template<class T>
T
compute_inverse_matrix_replace_row_col2(alps::numeric::matrix<T>& invG, const alps::numeric::matrix<T>& Dr, const alps::numeric::matrix<T>& Dc, int m,
                                       bool compute_only_det_rat) {
    typedef alps::numeric::matrix<T> matrix_t;
    const double eps = 1E-10;

    const int N = num_cols(invG);
    const int Nm1 = N-1;
    assert(num_rows(invG)==N);
    assert(num_rows(Dr)==1);
    assert(num_cols(Dr)==N);
    assert(num_rows(Dc)==N);
    assert(num_cols(Dc)==1);
    assert(Dr(0,m)==Dc(m,0));

    //original G^{-1}
    matrix_t tQ(Nm1, 1), tR(1, Nm1);
    invG.swap_cols(m,N-1);
    invG.swap_rows(m,N-1);
    copy_block(invG, 0, Nm1, tQ, 0, 0, Nm1, 1);
    copy_block(invG, Nm1, 0, tR, 0, 0, 1, Nm1);
    const T tS = invG(N-1,N-1);

    //new G
    matrix_t Q(Dc), R(Dr);
    Q.swap_rows(m,N-1); Q.resize(Nm1,1);
    R.swap_cols(m,N-1); R.resize(1,Nm1);

    const T R_tQ = mygemm(R, tQ)(0,0);
    const T tR_Q = mygemm(tR, Q)(0,0);
    const T lambda = R_tQ*tR_Q;

    if (compute_only_det_rat) {
        invG.swap_cols(m,N-1);
        invG.swap_rows(m,N-1);
        return lambda;
    }

    matrix_t tP(Nm1, Nm1);
    copy_block(invG, 0, 0, tP, 0, 0, Nm1, Nm1);

    matrix_t tQp(tQ); tQp /= R_tQ;
    matrix_t tRp(tR); tRp /= tR_Q;

    matrix_t tPp(tP);
    tPp += -mygemm(tQp,mygemm(R,tP))-mygemm(mygemm(tP,Q),tRp)+mygemm(mygemm(tQp,mygemm(R,mygemm(tP,Q))),tRp);

    copy_block(tPp, 0, 0, invG, 0, 0, Nm1, Nm1);
    copy_block(tQp, 0, 0, invG, 0, Nm1, Nm1, 1);
    copy_block(tRp, 0, 0, invG, Nm1, 0, 1, Nm1);
    const T tSp = tS/lambda;
    invG(N-1,N-1) = tSp;

    invG.swap_cols(m,N-1);
    invG.swap_rows(m,N-1);

    return lambda;
}


template<class T>
T
compute_inverse_matrix_replace_rows_cols_succesive(alps::numeric::matrix<T>& invG, const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                                                    const std::vector<int>& rows_cols, bool assume_intermediate_state_is_singular=false) {
    typedef alps::numeric::matrix<T> matrix_t;

    const int N = num_cols(invG);
    const int M = rows_cols.size();
    //const int Nm1 = N - 1;
    assert(num_rows(invG)==N);
    assert(num_rows(Q)==N-M && num_cols(Q)==M);
    assert(num_rows(R)==M && num_cols(R)==N-M);
    assert(num_rows(S)==M && num_cols(S)==M);
#ifndef NDEBUG
    for(int im=0; im<rows_cols.size()-1; ++im) {
        assert(rows_cols[im]<rows_cols[im+1]);
    }
#endif

    matrix_t Dr(1,N), Dc(N,1);
    //std::vector<int> idx_in_R(M);
    std::vector<bool> is_included(N,false);
    for (int im=0; im<M; ++im) {
        is_included[rows_cols[im]] = true;
    }
    //std::vector<int> rows_cols_rest(N-M);
    //int idx = 0;
    //for (int i=0; i<N; ++i) {
        //if (!is_included[i]) {
            //rows_cols_rest[idx] = i;
            //++idx;
        //}
    //}

    T lambda = 1.;
    for (int im=0; im<M; ++im) {
        int idx=0, idx2=0;
        for (int i=0; i<N; ++i) {
            if (!is_included[i]) {
                Dr(0,i) = R(im,idx);
                Dc(i,0) = Q(idx,im);
                ++idx;
            } else {
                Dr(0,i) = S(im,idx2);
                Dc(i,0) = S(idx2,im);
                ++idx2;
            }
        }
        assert(idx==N-M);
        assert(idx2==M);

        double rtmp = compute_inverse_matrix_replace_row_col(invG, Dr, Dc, rows_cols[im], assume_intermediate_state_is_singular);
        //std::cout << " im = " << im << " " << rtmp << std::endl;
        lambda *= rtmp;
    }
    return lambda;
}

template<class T>
T
compute_inverse_matrix_replace_single_row_col(alps::numeric::matrix<T>& invG, const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                                                   const std::vector<int>& rows_cols, bool compute_only_det) {
    typedef alps::numeric::matrix<T> matrix_t;

    const int N = num_cols(invG);
    const int M = rows_cols.size();
    if (M!=1)
        throw std::logic_error("Error in compute_inverse_matrix_replace_single_row_col");

    assert(num_rows(invG)==N);
    assert(num_rows(Q)==N-M && num_cols(Q)==M);
    assert(num_rows(R)==M && num_cols(R)==N-M);
    assert(num_rows(S)==M && num_cols(S)==M);

    matrix_t Dr(1,N), Dc(N,1);
    std::vector<bool> is_included(N,false);
    for (int im=0; im<M; ++im) {
        is_included[rows_cols[im]] = true;
    }

    T lambda = 1.;
    for (int im=0; im<M; ++im) {
        int idx=0, idx2=0;
        for (int i=0; i<N; ++i) {
            if (!is_included[i]) {
                Dr(0,i) = R(im,idx);
                Dc(i,0) = Q(idx,im);
                ++idx;
            } else {
                Dr(0,i) = S(im,idx2);
                Dc(i,0) = S(idx2,im);
                ++idx2;
            }
        }
        assert(idx==N-M);
        assert(idx2==M);

        T rtmp = compute_inverse_matrix_replace_row_col2(invG, Dr, Dc, rows_cols[im], compute_only_det);
        lambda *= rtmp;
    }
    return lambda;
}


template<class T>
void
replace_rows_cols(alps::numeric::matrix<T>& A,
                  const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                  const std::vector<int>& rows_cols) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const int NpM = num_cols(A);
    const int M = rows_cols.size();
    const int N = NpM-M;

    std::vector<std::pair<int,int> > swap_list(M);
    for (int i=0; i<M; ++i) {
        A.swap_cols(rows_cols[M-1-i], NpM-1-i);
        A.swap_rows(rows_cols[M-1-i], NpM-1-i);
        swap_list[i] = std::pair<int,int>(rows_cols[M-1-i], NpM-1-i);
    }

    copy_block(Q, 0, 0, A, 0, N, N, M);
    copy_block(R, 0, 0, A, N, 0, M, N);
    copy_block(S, 0, 0, A, N, N, M, M);

    for(std::vector<std::pair<int,int> >::reverse_iterator it=swap_list.rbegin(); it!=swap_list.rend(); ++it) {
        A.swap_cols(it->first, it->second);
        A.swap_rows(it->first, it->second);
    }
}

template<class T>
void generate_indices(const std::vector<T>& rows_cols, int N, int M, std::vector<T>& rows_cols_rest) {
    assert(rows_cols.size()==M);
    const int NpM = N+M;

    std::vector<bool> is_included(N+M,false);
    for (int im=0; im<M; ++im) {
        is_included[rows_cols[im]] = true;
    }

    rows_cols_rest.resize(N);
    int idx = 0;
    for (int i=0; i<NpM; ++i) {
        if (!is_included[i]) {
            rows_cols_rest[idx]  = i;
            ++idx;
        }
    }
    assert(idx==N);
}

template<class T>
void
replace_rows_cols_respect_ordering(alps::numeric::matrix<T>& A,
                  const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                  const std::vector<int>& rows_cols) {

    const int NpM = num_cols(A);
    const int M = rows_cols.size();
    const int N = NpM-M;

    std::vector<bool> is_included(N+M,false);
    for (int im=0; im<M; ++im) {
        is_included[rows_cols[im]] = true;
    }
    std::vector<int> pos_N(N), pos_M(M);
    int idx_N=0, idx_M=0;
    for (int i=0; i<NpM; ++i) {
        if (is_included[i]) {
            pos_M[idx_M] = i;
            ++idx_M;
        } else {
            pos_N[idx_N] = i;
            ++idx_N;
        }
    }
    assert(idx_N==N && idx_M==M);

    for (int j=0; j<M; ++j) {
        for (int i=0; i<M; ++i) {
            A(pos_M[i],pos_M[j]) =  S(i,j);
        }
    }
    for (int j=0; j<M; ++j) {
        for (int i=0; i<N; ++i) {
            A(pos_N[i],pos_M[j]) =  Q(i,j);
            A(pos_M[j],pos_N[i]) =  R(j,i);
        }
    }
}

//Implementing Ye-Hua Lie and Lei Wang (2015): Eqs. (17)-(26) before taking the limit of tS->0
//T = double or complex<double>
template<class T>
T
compute_det_ratio_replace_rows_cols(const alps::numeric::matrix<T>& invBigMat,
                               const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                               alps::numeric::matrix<T>& Mmat, alps::numeric::matrix<T>& inv_tSp) {
    //const std::vector<int>& rows_cols, const std::vector<std::pair<int,int> >& swap_list, alps::numeric::matrix<T>& Mmat, alps::numeric::matrix<T>& inv_tSp) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const int NpM = num_cols(invBigMat);
    const int M = num_rows(R);
    const int N = NpM-M;
    assert(N>0);

    assert(num_cols(invBigMat)==num_rows(invBigMat));
    assert(num_rows(R)==M && num_cols(R)==N);
    assert(num_rows(Q)==N && num_cols(Q)==M);
    assert(num_rows(S)==M && num_cols(S)==M);

    static matrix_t invtS_tR, inv_tS, MQ, RMQ, ws1, ws2;
    inv_tS.resize_values_not_retained(M,M);
    invtS_tR.resize_values_not_retained(M,N);
    MQ.resize_values_not_retained(N,M);
    RMQ.resize_values_not_retained(M,M);
    ws1.resize_values_not_retained(M,M);
    ws2.resize_values_not_retained(M,M);

    submatrix_view<T> tQ_view(invBigMat, 0, N, N, M);
    submatrix_view<T> tR_view(invBigMat, N, 0, M, N);
    submatrix_view<T> tS_view(invBigMat, N, N, M, M);

    //compute inv_tS
    copy_block(invBigMat, N, N, inv_tS, 0, 0, M, M);
    inverse_in_place(inv_tS);

    gemm(inv_tS, tR_view, invtS_tR);

    Mmat.resize_values_not_retained(N,N);
    copy_block(invBigMat, 0, 0, Mmat, 0, 0, N, N);
    mygemm((T)-1.0, tQ_view, invtS_tR, (T) 1.0, Mmat);

    inv_tSp.resize_values_not_retained(M,M);
    gemm(Mmat, Q, MQ);
    copy_block(S, 0, 0, inv_tSp, 0, 0, M, M);
    mygemm((T) -1.0, R, MQ, (T) 1.0, inv_tSp);
    return determinant_no_copy(tS_view, ws1)*determinant_no_copy(inv_tSp, ws2);
}

//Implementing Ye-Hua Lie and Lei Wang (2015): Eqs. (17)-(26) before taking the limit of tS->0
//T = double or complex<double>
template<class T>
T
compute_inverse_matrix_replace_rows_cols(alps::numeric::matrix<T>& invBigMat,
                                    const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                                    const alps::numeric::matrix<T>& Mmat, const alps::numeric::matrix<T>& inv_tSp) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const int NpM = num_cols(invBigMat);
    const int M = num_rows(R);
    const int N = NpM-M;

    assert(N>0);
    assert(num_cols(invBigMat)==num_rows(invBigMat));
    assert(num_rows(R)==M && num_cols(R)==N);
    assert(num_rows(Q)==N && num_cols(Q)==M);
    assert(num_rows(S)==M && num_cols(S)==M);

    static matrix_t tmp_NM, tmp_MN;
    tmp_NM.resize_values_not_retained(N,M);
    tmp_MN.resize_values_not_retained(M,N);
    submatrix_view<T> tPp_view(invBigMat,0,0,N,N);
    submatrix_view<T> tQp_view(invBigMat,0,N,N,M);
    submatrix_view<T> tRp_view(invBigMat,N,0,M,N);
    submatrix_view<T> tSp_view(invBigMat,N,N,M,M);

    //tSp
    my_copy_block(inv_tSp, 0, 0, tSp_view, 0, 0, M, M);
    inverse_in_place(tSp_view);

    //tQp
    gemm(Q,tSp_view,tmp_NM);
    mygemm((T)-1.0, Mmat, tmp_NM, (T) 0.0, tQp_view);

    //tRp
    gemm(tSp_view,R,tmp_MN);
    mygemm((T)-1.0, tmp_MN, Mmat, (T) 0.0, tRp_view);

    //tPp
    gemm(Mmat, Q, tmp_NM);
    my_copy_block(Mmat, 0, 0, tPp_view, 0, 0, N, N);
    mygemm((T)-1.0, tmp_NM, tRp_view, (T) 1.0, tPp_view);
}

template<class T>
T
compute_det_ratio_replace_rows_cols_safe(const alps::numeric::matrix<T>& invBigMat,
                                    const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                                    alps::numeric::matrix<T>& Mmat, alps::numeric::matrix<T>& inv_tSp) {
    using namespace alps::numeric;
    typedef matrix<T> matrix_t;

    const int NpM = num_cols(invBigMat);
    const int M = num_rows(R);
    const int N = NpM-M;
    assert(N>0);
    assert(M==1);

    assert(num_cols(invBigMat)==num_rows(invBigMat));
    assert(num_rows(R)==M && num_cols(R)==N);
    assert(num_rows(Q)==N && num_cols(Q)==M);
    assert(num_rows(S)==M && num_cols(S)==M);

    matrix_t tP(N, N), tQ(N, M), tR(M, N), tS(M, M), invtS_tR(M,N,0.), tQ_invtS_tR(N,N,0.);

    copy_block(invBigMat, 0, 0, tP, 0, 0, N, N);
    copy_block(invBigMat, 0, N, tQ, 0, 0, N, M);
    copy_block(invBigMat, N, 0, tR, 0, 0, M, N);
    copy_block(invBigMat, N, N, tS, 0, 0, M, M);

    return (mygemm(tS, S-mygemm(R,mygemm(tP,Q)))+mygemm(mygemm(R,tQ),mygemm(tR,Q)))(0,0);

    //gemm(inverse(tS), tR, invtS_tR);
    //gemm(tQ, invtS_tR, tQ_invtS_tR);
    //Mmat = tP-tQ_invtS_tR;

    //matrix_t MQ(N,M,0.), RMQ(M,M,0.);//, inv_tSp(M,M);
    //gemm(Mmat, Q, MQ);
    //gemm(R, MQ, RMQ);
    //inv_tSp = S-RMQ;
    //return determinant(tS)*determinant(inv_tSp);
}

template<typename T>
T
compute_det_ratio_replace_diaognal_elements(alps::numeric::matrix<T>& invBigMat, int num_elem_updated, const std::vector<int>& pos, const std::vector<T>& elems_diff, bool compute_only_det_rat) {
    using namespace alps::numeric;

    assert(invBigMat.num_cols()==invBigMat.num_rows());

    //work space (static!)
    static matrix<T> x, invC_plus_x, invC_plus_x_times_zx, yx;
    static std::vector<std::pair<int,int> > swap_list;

    const int N = invBigMat.num_cols();

    x.resize_values_not_retained(num_elem_updated, num_elem_updated);
    for (int j=0; j<num_elem_updated; ++j) {
        for (int i=0; i<num_elem_updated; ++i) {
            assert(pos[i]>=0 && pos[i]<invBigMat.num_cols());
            assert(pos[j]>=0 && pos[j]<invBigMat.num_cols());
            x(i,j) = invBigMat(pos[i],pos[j])*elems_diff[j];
        }
        x(j,j) += (T) 1;
    }

    const T det_rat = determinant(x);
    if (compute_only_det_rat) {
        return det_rat;
    }

    //swap rows and cols to move all diagonal elements to be updated to the end rows and cols
    swap_list.resize(num_elem_updated);
    for (int i=0; i<num_elem_updated; ++i) {
        invBigMat.swap_cols(pos[num_elem_updated-1-i], N-1-i);
        invBigMat.swap_rows(pos[num_elem_updated-1-i], N-1-i);
        swap_list[i] = std::pair<int,int>(pos[num_elem_updated-1-i], N-1-i);
    }
    yx.resize_values_not_retained(N, num_elem_updated);
    copy_block(invBigMat, 0, N-num_elem_updated, yx, 0, 0, N, num_elem_updated);

    submatrix_view<T> zx_view(invBigMat, N-num_elem_updated, 0,     num_elem_updated, N);
    invC_plus_x.resize_values_not_retained(num_elem_updated, num_elem_updated);
    invC_plus_x_times_zx.resize_values_not_retained(num_elem_updated, N);

    copy_block(invBigMat, N-num_elem_updated, N-num_elem_updated, invC_plus_x, 0, 0, num_elem_updated, num_elem_updated);
    for (int i=0; i<num_elem_updated; ++i) {
        invC_plus_x(i,i) += 1.0/elems_diff[i];
    }
    inverse_in_place(invC_plus_x);
    mygemm((T) 1.0, invC_plus_x, zx_view, (T) 0.0, invC_plus_x_times_zx);
    mygemm((T) -1.0, yx, invC_plus_x_times_zx, (T) 1.0, invBigMat);

    for(std::vector<std::pair<int,int> >::reverse_iterator it=swap_list.rbegin(); it!=swap_list.rend(); ++it) {
        invBigMat.swap_cols(it->first, it->second);
        invBigMat.swap_rows(it->first, it->second);
    }

    return det_rat;
}
#endif //IMPSOLVER_FASTUPDATE_FORMULA_H

