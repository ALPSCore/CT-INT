#include <algorithm>
#include "gtest/gtest.h"


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
    M_list.push_back(20);
    M_list.push_back(30);

    for (size_t n=0; n<N_list.size(); ++n) {
        for (size_t m=0; m<M_list.size(); ++m) {
            const size_t N = N_list[n];
            const size_t M = M_list[m];

            typedef alps::numeric::matrix<double> matrix_t;

            matrix_t BigMatrix(N+M, N+M, 0), invBigMatrix(N+M, N+M, 0);
            matrix_t SmallMatrix(N,N,0);
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

TEST(Boost, Binomial) {
    const size_t k = 2;
    for (size_t N=k; N<10; ++N) {
        const double tmp = boost::math::binomial_coefficient<double>(N,k);
        ASSERT_TRUE(std::abs(tmp-0.5*N*(N-1.0))<1E-8);
    }
}

TEST(MyUtil, permutation) {
    assert(permutation(3,1)==3);
    assert(permutation(3,2)==6);
}

TEST(QuantumNumber, diagonal_GF) {
    size_t n_site = 4;
    size_t n_flavors = 2;
    size_t N=1000;

    const size_t n_rank=2;
    const size_t n_af=2;
    const size_t Nv=1;
    const double eps=0;

    green_function<double> gf(N, n_site, n_flavors);

    std::vector<vertex_definition<double> > vertices;
    boost::multi_array<double,2> alpha(boost::extents[n_af][n_rank]);
    std::fill(alpha.origin(), alpha.origin()+alpha.num_elements(), 0);

    //for (size_t iv=0; iv<Nv; ++iv) {
    std::vector<spin_t> flavors(n_rank);
    std::vector<size_t> sites(2*n_rank);
    flavors[0] = 0;
    flavors[1] = 1;
    sites[0] = 0;
    sites[1] = 1;
    sites[2] = 2;
    sites[3] = 3;
    vertices.push_back(vertex_definition<double>(2,2,flavors,sites,0.0,alpha));

    int qs[] = {1, -1, 0, 0, 0, 0, 1, -1};
    std::vector<std::valarray<int> > quantum_number_vertices;
    quantum_number_vertices = make_quantum_numbers(gf, vertices, eps);
    std::valarray<int> qs2(qs,n_site*n_flavors);
    //std::valarray<int> qs3(qs,n_site*n_flavors);
    bool flag = true;
    for (int i=0; i<qs2.size(); ++i) {
        //std::cout << i << " " << quantum_number_vertices[0][i] << std::endl;
        if (qs2[i]!=quantum_number_vertices[0][i]) {
            flag = false;
        }
    }
    ASSERT_TRUE(flag);
    //ASSERT_TRUE(qs2==std::valarray<int>(quantum_number_vertices[0]));
    //ASSERT_TRUE(1==1);
    //ASSERT_TRUE(qs2==qs3);
}

TEST(UpdateStatistics, EstimateSpread) {
    const double beta = 100;

    std::vector<simple_vertex> vertices;
    vertices.push_back(simple_vertex(0.0));
    vertices.push_back(simple_vertex(0.3*beta));
    ASSERT_TRUE(std::abs(compute_spread(vertices, beta)/beta-0.3)<1E-5);

    //vertices are distributed on [0,beta] uniformly.
    for (int Nv=2; Nv<10; ++Nv) {
        vertices.clear();
        for (int iv=0; iv<Nv; ++iv) {
            vertices.push_back(simple_vertex((beta*iv)/Nv));
        }
        ASSERT_TRUE(std::abs(compute_spread(vertices, beta)/beta-(1-1./Nv))<1E-5);
    }
}
