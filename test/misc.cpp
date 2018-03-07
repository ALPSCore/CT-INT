#include <algorithm>

#include <complex>
#include <limits>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <boost/range/irange.hpp>

#include "../src/legendre.h"
#include "../src/util.h"
#include "../src/green_function.h"

#include "gtest.h"
#include "common.hpp"

using namespace alps::ctint;

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

/*
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
*/

TEST(Util, Mymod) {
    const double tau = 1.0, beta = 15.0;
    for (int i=-10; i<=10; ++i) {
        double x = tau+beta*i;
        ASSERT_TRUE(std::abs(mymod(x,beta)-tau)<1E-10);
    }
}


