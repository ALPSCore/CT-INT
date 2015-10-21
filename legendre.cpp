//
// Created by H. Shinaoka on 2015/10/21.
//

#include "legendre.h"


LegendreTransformer::LegendreTransformer(int n_matsubara, int n_legendre)
        : n_matsubara_(n_matsubara), n_legendre_(n_legendre), Tnl_(n_matsubara,n_legendre) {

    double sign_tmp = 1.0;
    for (int im=0; im<n_matsubara_; ++im) {
        std::complex<double> ztmp(0.0,1.0);
        for (int il=0; il<n_legendre_; ++ il) {
            Tnl_(im,il) = sign_tmp*ztmp*std::sqrt(2*il+1.0)*boost::math::sph_bessel(il,0.5*(2*im+1)*M_PI);
            ztmp *= std::complex<double>(0.0,1.0);
        }
        sign_tmp *= -1;
    }
    //do something
};

void
LegendreTransformer::compute_legendre(double x, std::vector<double> &val) const {
    assert(val.size()>=n_legendre_);
    for(int l=0; l<n_legendre_; l++) {
        if (l == 0) {
            val[l] = 1;
        } else if (l == 1) {
            val[l] = x;
        } else {
            val[l] = ((2 * l - 1) * x * val[l-1] - (l - 1) * val[l-2]) / static_cast<double>(l);//l
        }
    }
}

const alps::numeric::matrix<std::complex<double> > &LegendreTransformer::Tnl() const {
    return Tnl_;
}
