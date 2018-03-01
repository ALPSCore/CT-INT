#include <algorithm>
#include "unittests.hpp"

#include "gtest.h"

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
    for (int t=0; t<N; ++t) {
        for (int flavor=0; flavor<n_flavors; ++flavor) {
            for (int i=0; i<n_site; ++i) {
                for (int j=0; j<n_site; ++j) {
                    gf(t, i, j, flavor) = 0.0;
                }
                gf(t, i, i, flavor) = -0.5;
            }
        }
    }

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
    vertices.push_back(vertex_definition<double>(2,2,flavors,sites,0.0,alpha,0));

    int qs[] = {1, -1, 0, 0, 0, 0, 1, -1};
    std::vector<std::vector<std::valarray<int> > > quantum_number_vertices;
    std::vector<std::vector<std::vector<size_t> > > groups(n_flavors);
    std::vector<std::vector<int> > group_map;
    quantum_number_vertices = make_quantum_numbers<double,double,green_function<double> >(gf, vertices, groups, group_map, eps);
    std::valarray<int> qs2(qs,n_site*n_flavors);
    bool flag = true;
    int i_af = 0;
    for (int i=0; i<qs2.size(); ++i) {
        if (qs2[i]!=quantum_number_vertices[0][i_af][i]) {
            flag = false;
        }
    }
    ASSERT_TRUE(flag);

    //vertices[0].make_quantum_numbers(group_map, quantum_number_vertices[0].size()/n_flavors);

    //ASSERT_TRUE(qs2==std::valarray<int>(quantum_number_vertices[0]));
    //ASSERT_TRUE(1==1);
    //ASSERT_TRUE(qs2==qs3);
}

TEST(UpdateStatistics, EstimateSpread) {
    const double beta = 100;

    std::vector<simple_vertex> vertices;
    vertices.push_back(simple_vertex(0.0));
    vertices.push_back(simple_vertex(0.3*beta));
    ASSERT_TRUE(std::abs(compute_spread<std::vector<simple_vertex> >(vertices.begin(), vertices.end(), beta)/beta-0.3)<1E-5);

    //vertices are distributed on [0,beta] uniformly.
    for (int Nv=2; Nv<10; ++Nv) {
        vertices.clear();
        for (int iv=0; iv<Nv; ++iv) {
            vertices.push_back(simple_vertex((beta*iv)/Nv));
        }
        ASSERT_TRUE(std::abs(compute_spread<std::vector<simple_vertex> >(vertices.begin(), vertices.end(), beta)/beta-(1-1./Nv))<1E-5);
    }
}

TEST(TypeRange, Integer) {
    ASSERT_TRUE(std::numeric_limits<int>::max()>=2147483647);
}

TEST(Util, Mymod) {
    const double tau = 1.0, beta = 15.0;
    for (int i=-10; i<=10; ++i) {
        double x = tau+beta*i;
        ASSERT_TRUE(std::abs(mymod(x,beta)-tau)<1E-10);
    }
}

TEST(SubmatrixUpdate, single_vertex_insertion_spin_flip)
{
    typedef std::complex<double> T;
    const int n_sites = 3;
    const double U = 2.0;
    const double alpha = 1E-2;
    const double beta = 200.0;
    const int Nv_max = 2;
    const int n_flavors = 2;
    const int k_ins_max = 32;
    const int n_update = 5;
    const int seed = 100;

    std::vector<double> E(n_sites);
    boost::multi_array<T,2> phase(boost::extents[n_sites][n_sites]);

    for (int i=0; i<n_sites; ++i) {
        E[i] = (double) i;
        //std::cout << phase[i] << std::endl;
    }
    for (int i=0; i<n_sites; ++i) {
        for (int j=i; j<n_sites; ++j) {
            phase[i][j] = std::exp(std::complex<double>(0.0, 1.*i*(2*j+1.0)));
            phase[j][i] = myconj(phase[i][j]);
        }
    }

    general_U_matrix<T> Uijkl(n_sites, U, alpha);

    itime_vertex_container itime_vertices_init;
    itime_vertices_init.push_back(itime_vertex(0, 0, 0.5*beta, 2, true));

    /* initialize submatrix_update */
    //SubmatrixUpdate<T> submatrix_update(k_ins_max, n_flavors, DiagonalG0<T>(beta), &Uijkl, beta, itime_vertices_init);
    SubmatrixUpdate<T> submatrix_update(k_ins_max, n_flavors, OffDiagonalG0<T>(beta, n_sites, E, phase), &Uijkl, beta, itime_vertices_init);

    submatrix_update.sanity_check();

    /* init udpate_manager */
    alps::params params;
    params["BETA"] = beta;
    params["FLAVORS"] = n_flavors;
    params["N_MULTI_VERTEX_UPDATE"] = Nv_max;
    params["DOUBLE_VERTEX_UPDATE_A"] = 1.0/beta;
    params["DOUBLE_VERTEX_UPDATE_B"] = 1.0e-2;
    VertexUpdateManager<T> manager(params, Uijkl, OffDiagonalG0<T>(beta, n_sites, E, phase), false);

    /* initialize RND generator */
    //std::vector<double> probs(Nv_max, 1.0);
    boost::random::uniform_smallint<> dist(1,Nv_max);
    boost::random::uniform_01<> dist01;
    boost::random::mt19937 gen(seed);
    boost::random::variate_generator<boost::random::mt19937&, boost::random::uniform_smallint<> > Nv_prob(gen, dist);
    boost::random::variate_generator<boost::random::mt19937&, boost::random::uniform_01<> > random01(gen, dist01);

    std::vector<alps::numeric::matrix<T> > M(n_flavors), M_scratch(n_flavors);


    for (int i_update=0; i_update<n_update; ++i_update) {
        T sign_from_M0, weight_from_M0;
        boost::tie(sign_from_M0,weight_from_M0) = submatrix_update.compute_M_from_scratch(M_scratch);

        //const T weight_rat = submatrix_update.vertex_insertion_removal_update(manager, random01);
        const T weight_rat = manager.do_ins_rem_update(submatrix_update, Uijkl, random01, 1.0);
        const T sign_bak = submatrix_update.sign();

        ASSERT_TRUE(submatrix_update.sanity_check());
        submatrix_update.recompute_matrix(true);
        submatrix_update.compute_M(M);
        T sign_from_M, weight_from_M;
        boost::tie(sign_from_M,weight_from_M) = submatrix_update.compute_M_from_scratch(M_scratch);

        ASSERT_TRUE(my_equal(weight_from_M/weight_from_M0, weight_rat, 1E-5));

        //std::cout << "sign " << sign_bak << " " << submatrix_update.sign() << std::endl;
        //std::cout << "sign_from_M " << sign_from_M << std::endl;
        ASSERT_TRUE(std::abs(sign_bak-submatrix_update.sign())<1.0e-5);
        ASSERT_TRUE(std::abs(sign_from_M-submatrix_update.sign())<1.0e-5);
        //std::cout << " Nv " << submatrix_update.itime_vertices().num_interacting() << std::endl;
        for (int flavor=0; flavor<n_flavors; ++flavor) {
            if (M[flavor].size2()>0) {
                ASSERT_TRUE(alps::fastupdate::norm_square(M[flavor]-M_scratch[flavor])/alps::fastupdate::norm_square(M[flavor])<1E-8);
            }
        }

        const T weight_rat2 = manager.do_spin_flip_update(submatrix_update, Uijkl, random01);

        T sign_from_M2, weight_from_M2;
        boost::tie(sign_from_M2,weight_from_M2) = submatrix_update.compute_M_from_scratch(M_scratch);
        ASSERT_TRUE(my_equal(weight_from_M2/weight_from_M, weight_rat2, 1E-5));

        const T weight_rat3 = manager.do_shift_update(submatrix_update, Uijkl, random01, false);
    }

    //std::cout << DiagonalG0<T>(beta)(0.0) << std::endl;
    //std::cout << DiagonalG0<T>(beta)(0.5*beta) << std::endl;
    //std::cout << DiagonalG0<T>(beta)(0.9999*beta) << std::endl;
}
