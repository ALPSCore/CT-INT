#include <algorithm>
#include "gtest/gtest.h"

#include <boost/random.hpp>
#include <boost/multi_array.hpp>

#include "gtest.hpp"

TEST(LegendreMeasurement, Tnl)
{
    const int n_matsubara = 1;
    const int n_legendre = 1000;
    LegendreTransformer trans(n_matsubara, n_legendre);

    //Check if norm is conserved
    for (std::size_t im=0; im<n_matsubara; ++im) {
        std::complex<double> tmp;
        for (std::size_t il=0; il<n_legendre; ++il) {
            tmp += trans.Tnl()(im,il)*std::conj(trans.Tnl()(im,il));
        }
        ASSERT_TRUE(std::abs(tmp-1.)<1E-5) << "norm is not conserved = " << tmp;
    }

    std::vector<double> legendre_values(n_legendre);
    const double xval = 0.2;
    trans.compute_legendre(xval, legendre_values);
    ASSERT_TRUE(std::abs(legendre_values[1]-xval)<1E-8);
    for (std::size_t il=1; il<n_legendre-1; ++il) {
        double left_val = (il+1)*legendre_values[il+1];
        double right_val = (2*il+1)*xval*legendre_values[il]-il*legendre_values[il-1];
        ASSERT_TRUE(std::abs(left_val-right_val)<1E-8);
    }
}

TEST(FastUpdate, BlockMatrixAdd)
{
    std::vector<size_t> N_list, M_list;
    N_list.push_back(0);
    N_list.push_back(10);
    M_list.push_back(10);
    M_list.push_back(20);

    for (size_t n=0; n<N_list.size(); ++n) {
        for (size_t m=0; m<M_list.size(); ++m) {
            const size_t N = N_list[n];
            const size_t M = M_list[m];

            typedef alps::numeric::matrix<double> matrix_t;

            matrix_t A(N,N,0), B(N,M,0), C(M,N,0), D(M,M,0), invA(N,N,0);
            matrix_t E(N,N,0), F(N,M,0), G(M,N,0), H(M,M,0);
            matrix_t BigMatrix(N+M, N+M, 0);
            matrix_t invBigMatrix2(N+M, N+M, 0);

            randomize_matrix(A, 100);//100 is a seed
            randomize_matrix(B, 200);
            randomize_matrix(C, 300);
            randomize_matrix(D, 400);
            if (N>0) {
                invA = alps::numeric::inverse(A);
            } else {
                invA.resize(0,0);
            }

            alps::numeric::copy_block(A,0,0,BigMatrix,0,0,N,N);
            alps::numeric::copy_block(B,0,0,BigMatrix,0,N,N,M);
            alps::numeric::copy_block(C,0,0,BigMatrix,N,0,M,N);
            alps::numeric::copy_block(D,0,0,BigMatrix,N,N,M,M);

            matrix_t invBigMatrix(alps::numeric::inverse(BigMatrix));

            double det_rat = compute_inverse_matrix_up(B, C, D, invA, E, F, G, H);
            ASSERT_TRUE(std::abs(det_rat-alps::numeric::determinant(BigMatrix)/alps::numeric::determinant(A))<1E-8) << "N=" << N << " M=" << M;

            alps::numeric::copy_block(E,0,0,invBigMatrix2,0,0,N,N);
            alps::numeric::copy_block(F,0,0,invBigMatrix2,0,N,N,M);
            alps::numeric::copy_block(G,0,0,invBigMatrix2,N,0,M,N);
            alps::numeric::copy_block(H,0,0,invBigMatrix2,N,N,M,M);

            ASSERT_TRUE(std::abs(det_rat-alps::numeric::determinant(BigMatrix)/alps::numeric::determinant(A))) << "N=" << N << " M=" << M;
            ASSERT_TRUE(alps::numeric::norm_square(invBigMatrix-invBigMatrix2)<1E-8);

            det_rat = compute_det_ratio_up(B, C, D, invA);
            ASSERT_TRUE(std::abs(det_rat-alps::numeric::determinant(BigMatrix)/alps::numeric::determinant(A))<1E-8) << "N=" << N << " M=" << M;

            det_rat = compute_inverse_matrix_up2(B, C, D, invA, invBigMatrix2);
            ASSERT_TRUE(alps::numeric::norm_square(invBigMatrix-invBigMatrix2)<1E-8) << "N=" << N << " M=" << M;
        }
    }
}

TEST(FastUpdate, BlockMatrixRemove)
{
    std::vector<size_t> N_list, M_list;
    //N_list.push_back(0);
    N_list.push_back(10);
    //M_list.push_back(10);
    M_list.push_back(10);

    for (size_t n=0; n<N_list.size(); ++n) {
        for (size_t m=0; m<M_list.size(); ++m) {
            const size_t N = N_list[n];
            const size_t M = M_list[m];

            typedef alps::numeric::matrix<double> matrix_t;

            matrix_t BigMatrix(N+M, N+M, 0), invBigMatrix(N+M, N+M, 0);
            matrix_t SmallMatrix(N,N,0), invSmallMatrix(N,N,0);
            std::vector<std::pair<size_t,size_t> > swap_list;

            randomize_matrix(BigMatrix, 100);//100 is a seed
            invBigMatrix = inverse(BigMatrix);

            //which rows and cols are to be removed
            std::vector<size_t> rows_removed(N+M);
            std::vector<size_t> rows_remain(N);
            for (size_t i=0; i<N+M; ++i) {
                rows_removed[i] = i;
            }
            std::random_shuffle(rows_removed.begin(), rows_removed.end());
            for (size_t i=0; i<N; ++i) {
                rows_remain[i] = rows_removed[i+M];
            }
            rows_removed.resize(M);
            std::sort(rows_removed.begin(), rows_removed.end());
            std::sort(rows_remain.begin(), rows_remain.end());

            for (size_t j=0; j<N; ++j) {
                for (size_t i=0; i<N; ++i) {
                    SmallMatrix(i,j) = BigMatrix(rows_remain[i], rows_remain[j]);
                }
            }

            //testing compute_det_ratio_down
            double det_rat = compute_det_ratio_down(M,rows_removed,invBigMatrix);
            ASSERT_TRUE(std::abs(det_rat-alps::numeric::determinant(SmallMatrix)/alps::numeric::determinant(BigMatrix))<1E-8) << "N=" << N << " M=" << M;

            matrix_t invSmallMatrix2(invBigMatrix);
            double det_rat2 = compute_inverse_matrix_down(M,rows_removed,invSmallMatrix2, swap_list);
            ASSERT_TRUE(std::abs(det_rat-det_rat2)<1E-8) << "N=" << N << " M=" << M;

            matrix_t SmallMatrix3(BigMatrix);
            for (size_t s=0; s<swap_list.size(); ++s) {
                SmallMatrix3.swap_cols(swap_list[s].first, swap_list[s].second);
                SmallMatrix3.swap_rows(swap_list[s].first, swap_list[s].second);
            }
            SmallMatrix3.resize(N,N);
            ASSERT_TRUE(alps::numeric::norm_square(inverse(SmallMatrix3)-invSmallMatrix2)<1E-8) << "N=" << N << " M=" << M;
        }
    }
}
