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
#include "../src/spline.h"

#include "gtest.h"
#include "common.hpp"

using namespace alps::ctint;

TEST(Spline, uniform_mesh) {
    tk::spline spline;

    std::size_t N = 100;

    std::vector<double> x(N), y(N);

    for (int i=0; i<N; ++i) {
        x[i] = static_cast<double>(i)/(N-1);
        y[i] = x[i] * x[i];
    }

    spline.set_points(x, y);

  //do some test

}

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

TEST(Util, Mymod) {
    const double tau = 1.0, beta = 15.0;
    for (int i=-10; i<=10; ++i) {
        double x = tau+beta*i;
        ASSERT_TRUE(std::abs(mymod(x,beta)-tau)<1E-10);
    }
}


