#pragma once

#include <complex>
#include "types.h" // for the multiple_vector_type
#include "green_function.h" // for the multiple_vector_type

#include <alps/params.hpp>

namespace alps {
    namespace ctint {


        template<typename T>
        std::complex<double> f_omega(std::complex<double> iw, T c1, T c2, T c3) {
          std::complex<double> iwsq = iw * iw;
          return c1 / iw + c2 / (iwsq) + c3 / (iw * iwsq);
        }


        template<typename T>
        double f_tau(double tau, double beta, T c1, T c2, T c3) {
          return -0.5 * c1 + (c2 * 0.25) * (-beta + 2. * tau) + (c3 * 0.25) * (beta * tau - tau * tau);
        }


        class FourierTransformer {
        public:
            typedef green_function<std::complex<double> > matsubara_green_function_t;
            //typedef green_function<double> itime_green_function_t;
            typedef std::complex<double> GTAU_TYPE;
            typedef green_function<GTAU_TYPE> itime_green_function_t;

            FourierTransformer(const double beta, const int n_flavor, const int n_site) {
              beta_ = beta;
              //mu_=mu;
              c1_.resize(n_flavor);
              c2_.resize(n_flavor);
              c3_.resize(n_flavor);
              Sc0_.resize(n_flavor);
              Sc1_.resize(n_flavor);
              Sc2_.resize(n_flavor);
              for (int f = 0; f < n_flavor; ++f) {
                c1_[f].resize(n_site);
                c2_[f].resize(n_site);
                c3_[f].resize(n_site);
                Sc0_[f].resize(n_site);
                Sc1_[f].resize(n_site);
                Sc2_[f].resize(n_site);
                for (int i = 0; i < n_site; ++i) {
                  c1_[f][i].resize(n_site);
                  c2_[f][i].resize(n_site);
                  c3_[f][i].resize(n_site);
                  Sc0_[f][i].resize(n_site);
                  Sc1_[f][i].resize(n_site);
                  Sc2_[f][i].resize(n_site);
                  for (int j = 0; j < n_site; ++j) {
                    c1_[f][i][j] = (i == j) ? 1. : 0.;
                    c2_[f][i][j] = 0.;
                    c3_[f][i][j] = 0.;
                    Sc0_[f][i][j] = 0.;
                    Sc1_[f][i][j] = 0.;
                    Sc2_[f][i][j] = 0.;
                  }
                }
              }
            }


            virtual ~FourierTransformer() {}

            virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const;

            virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const;

            //virtual void backward_ft_fullG(itime_full_green_function_t &G_tau, const matsubara_full_green_function_t &G_omega) const;
            virtual void append_tail(matsubara_green_function_t &G_omega, const matsubara_green_function_t &G0_omega,
                                     const int nfreq_measured) const;

            static void generate_transformer_lowest_order(const alps::params &parms,
                                                          boost::shared_ptr<FourierTransformer> &fourier_ptr);

        protected:

            double beta_;
            std::vector<std::vector<std::vector<double> > > c1_;
            std::vector<std::vector<std::vector<double> > > c2_;
            std::vector<std::vector<std::vector<double> > > c3_;
            std::vector<std::vector<std::vector<double> > > Sc0_; //coefficients for the self-energy
            std::vector<std::vector<std::vector<double> > > Sc1_;
            std::vector<std::vector<std::vector<double> > > Sc2_;


        };

        class SimpleG0FourierTransformer : public FourierTransformer {
        public:
            SimpleG0FourierTransformer(const double beta, const double mu, const double h, const int n_flavor,
                                       const std::vector<double> &eps, const std::vector<double> &epssq)
              : FourierTransformer(beta, n_flavor, 1) {
              for (int f = 0; f < n_flavor; f++) {
                double hsign = f % 2 ? h : -h;
                double mub = mu + hsign;
                c1_[f][0][0] = 1.;
                c2_[f][0][0] = eps[f] - mub;
                c3_[f][0][0] = epssq[f] - 2 * mub * eps[f] + mub * mub;
                // Note: the coefficients are the same if h is uniform field as if it is staggered
              }
            }
        };
    }
}
