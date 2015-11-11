//
// Created by H. Shinaoka on 2015/10/28.
//

#ifndef IMPSOLVER_FASTUPDATE_FORMULA_H
#define IMPSOLVER_FASTUPDATE_FORMULA_H

//#include <boost/multiprecision/float128.hpp>

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

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
    return determinant(H);
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
    /*
    std::vector<size_t> new_index(NpM);
    {
        std::vector<int> mark(NpM, 0);
        //put 1 on rows to be removed
        for (size_t i=0; i<M; ++i) {
            mark[rows_cols_removed[i]] = 1;
        }
        int pos_remain = 0;
        int pos_removed = N;
        for (size_t i=0; i<NpM; ++i) {

        }
    }
     */

    swap_list.resize(M);
    for (size_t i=0; i<M; ++i) {
        //if(rows_cols_removed[M-1-i]!=NpM-1-i) {
        invBigMat.swap_cols(rows_cols_removed[M-1-i], NpM-1-i);
        invBigMat.swap_rows(rows_cols_removed[M-1-i], NpM-1-i);
        swap_list[i] = std::pair<size_t,size_t>(rows_cols_removed[M-1-i], NpM-1-i);
        //}
    }

    if (N==0) {
        matrix_t H(invBigMat);
        invBigMat.resize(0,0);
        return determinant(H);
    } else {
        matrix_t E(N, N), F(N, M), G(M, N), H(M, M);
        copy_block(invBigMat, 0, 0, E, 0, 0, N, N);
        copy_block(invBigMat, 0, N, F, 0, 0, N, M);
        copy_block(invBigMat, N, 0, G, 0, 0, M, N);
        copy_block(invBigMat, N, N, H, 0, 0, M, M);

        matrix_t invH_G(M, N, 0), F_invH_G(N, N, 0);//one might reuse memories...
        gemm(inverse(H), G, invH_G);
        gemm(F, invH_G, F_invH_G);

        invBigMat.resize(N, N);
        invBigMat = E - F_invH_G;
        return determinant(H);
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
compute_inverse_matrix_replace_row_col(alps::numeric::matrix<T>& invG, const alps::numeric::matrix<T> Dr, const alps::numeric::matrix<T> Dc, int m) {
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
    if (std::abs(tS)<eps && std::abs(tSp)<eps) {
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

template<class T>
T
compute_inverse_matrix_replace_rows_cols_succesive(alps::numeric::matrix<T>& invG, const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                                                    const std::vector<int>& rows_cols) {
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

        double rtmp = compute_inverse_matrix_replace_row_col(invG, Dr, Dc, rows_cols[im]);
        //std::cout << " im = " << im << " " << rtmp << std::endl;
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

//Implementing Ye-Hua Lie and Lei Wang (2015): Eqs. (17)-(26).
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

    matrix_t tP(N, N), tQ(N, M), tR(M, N), tS(M, M), invtS_tR(M,N,0.), tQ_invtS_tR(N,N,0.);

    copy_block(invBigMat, 0, 0, tP, 0, 0, N, N);
    copy_block(invBigMat, 0, N, tQ, 0, 0, N, M);
    copy_block(invBigMat, N, 0, tR, 0, 0, M, N);
    copy_block(invBigMat, N, N, tS, 0, 0, M, M);

    const double norm = std::sqrt(alps::numeric::norm_square(tS)/(M*M));
    matrix_t tS_norm(tS);
    tS_norm /= norm;

    gemm(inverse(tS_norm), tR, invtS_tR);
    gemm(tQ, invtS_tR, tQ_invtS_tR);
    Mmat = tP-(1/norm)*tQ_invtS_tR;

    matrix_t MQ(N,M,0.), RMQ(M,M,0.);//, inv_tSp(M,M);
    gemm(Mmat, Q, MQ);
    gemm(R, MQ, RMQ);
    inv_tSp = S-RMQ;
    const double norm2 = std::sqrt(alps::numeric::norm_square(inv_tSp)/(M*M));
    matrix_t inv_tSp_norm(inv_tSp);
    inv_tSp_norm /= norm2;
    return determinant(tS_norm)*determinant(inv_tSp_norm)*pow(norm*norm2,(double)M);
}

template<class T>
T
compute_inverse_matrix_replace_rows_cols(alps::numeric::matrix<T>& invBigMat,
                                    const alps::numeric::matrix<T>& Q, const alps::numeric::matrix<T>& R, const alps::numeric::matrix<T>& S,
                                    const alps::numeric::matrix<T>& Mmat, const alps::numeric::matrix<T>& inv_tSp,
                                    //const std::vector<int>& rows_cols, const std::vector<std::pair<int,int> >& swap_list, const alps::numeric::matrix<T>& Mmat, const alps::numeric::matrix<T>& inv_tSp,
                                    alps::numeric::matrix<T>& tPp, alps::numeric::matrix<T>& tQp, alps::numeric::matrix<T>& tRp, alps::numeric::matrix<T>& tSp) {
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

    matrix_t tmp_NM(N,M,0.), tmp_MN(M,N,0.);
    tPp.resize(N,N);
    tQp.resize(N,M);
    tRp.resize(M,N);
    //tSp.resize(M,M);

    //tSp
    tSp = safe_inverse(inv_tSp);

    //tQp
    gemm(Q,tSp,tmp_NM);
    std::fill(tQp.get_values().begin(), tQp.get_values().end(), 0.0);
    gemm(Mmat,tmp_NM,tQp);
    tQp *= -1;

    //tRp
    gemm(tSp,R,tmp_MN);
    std::fill(tRp.get_values().begin(), tRp.get_values().end(), 0.0);
    gemm(tmp_MN,Mmat,tRp);
    tRp *= -1;

    //tPp
    std::fill(tPp.get_values().begin(), tPp.get_values().end(), 0.0);
    std::fill(tmp_NM.get_values().begin(), tmp_NM.get_values().end(), 0.0);
    gemm(Mmat, Q, tmp_NM);
    gemm(tmp_NM, tRp, tPp);
    tPp = Mmat-tPp;

    //write back
    copy_block(tPp, 0, 0, invBigMat, 0, 0, N, N);
    copy_block(tQp, 0, 0, invBigMat, 0, N, N, M);
    copy_block(tRp, 0, 0, invBigMat, N, 0, M, N);
    copy_block(tSp, 0, 0, invBigMat, N, N, M, M);

    ////swap rows and cols
    //for (std::vector<std::pair<int,int> >::const_reverse_iterator it=swap_list.rbegin(); it!=swap_list.rend(); ++it) {
        //invBigMat.swap_cols(it->first, it->second);
        //invBigMat.swap_rows(it->first, it->second);
    //}
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
#endif //IMPSOLVER_FASTUPDATE_FORMULA_H
