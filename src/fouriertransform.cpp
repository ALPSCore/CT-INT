#include "fouriertransform.h"
#include <valarray>

#include <math.h>
#include "util.h"

namespace alps {
    namespace ctint {


        void FourierTransformer::backward_ft(itime_green_function_t &G_tau,
                                             const matsubara_green_function_t &G_omega) const {
          /*
          assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
          unsigned int N_tau = G_tau.ntime()-1;
          unsigned int N_omega = G_omega.nfreq();
          unsigned int N_site = G_omega.nsite();
          matsubara_green_function_t G_omega_no_model(G_omega);
          itime_green_function_t G_tau_no_model(G_tau);
          const double dt = beta_/N_tau;

          for(unsigned int f=0;f<G_omega.nflavor();++f){
            for (unsigned int s1=0; s1<N_site; ++s1){
              for (unsigned int s2=0; s2<N_site; ++s2) {
                for (unsigned int k=0; k<N_omega; k++) {
                  std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
                  G_omega_no_model(k,s1,s2,f) -= f_omega(iw, c1_[f][s1][s2], c2_[f][s1][s2], c3_[f][s1][s2]);
                }
                for (unsigned int i=0; i<N_tau+1; i++) {
                  G_tau_no_model(i,s1,s2,f) = 0.0;
                  for (unsigned int k=0; k<N_omega; k++) {
                    const double wt((2*k+1)*i*M_PI/N_tau);
                    G_tau_no_model(i,s1,s2,f) += 1/beta_*G_omega_no_model(k,s1,s2,f)*std::complex<double>(cos(wt),-sin(wt));
                  }
                }
              }
            }
          }
          for(unsigned int f=0;f<G_omega.nflavor();++f){
            for (unsigned int i=0; i<N_tau+1; i++) {
              for (unsigned int s1 = 0; s1 < N_site; ++s1) {
                for (unsigned int s2 = 0; s2 < N_site; ++s2) {
                  G_tau(i, s1, s2, f) = f_tau(i * dt, beta_, c1_[f][s1][s2], c2_[f][s1][s2], c3_[f][s1][s2])
                                        +G_tau_no_model(i, s1, s2, f)+std::conj(G_tau_no_model(i, s2, s1, f));
                }
              }
            }
          }
           */
        }

        void FourierTransformer::forward_ft(const itime_green_function_t & gtau, matsubara_green_function_t & gomega) const
        {
          throw std::runtime_error("Complex version not implemented yet!");
        }

        void FourierTransformer::append_tail(matsubara_green_function_t& G_omega,
                                             const matsubara_green_function_t& G0_omega,
                                             const int nfreq_measured) const
        {
          /*
          for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
            for (site_t k=0; k<G0_omega.nsite(); ++k) {
              std::cout << "append tail to self-energy with coefficients: "
                        << " " << Sc0_[flavor][k][k]
                        << " " << Sc1_[flavor][k][k]
                        << " " << Sc2_[flavor][k][k] << std::endl;
              for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
                std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
                std::complex<double> Sigma = Sc0_[flavor][k][k] + Sc1_[flavor][k][k]/iw + Sc2_[flavor][k][k]/(iw*iw);
                G_omega(freq, k, k, flavor) = 1./(1./G0_omega(freq, k, k, flavor) - Sigma);
              }
            }
          }
           */
        }

        void FourierTransformer::generate_transformer_lowest_order(const alps::params &parms,
                                                                   boost::shared_ptr<FourierTransformer> &fourier_ptr)
        {
          /*
          int n_flavors = parms["FLAVORS"];
          int n_site = parms["SITES"];
          std::cout << "using general fourier transformer (up to lowest order)" << "\n";
          fourier_ptr.reset(new FourierTransformer((double)parms["BETA"], n_flavors, n_site));
           */
        }
    }
}
