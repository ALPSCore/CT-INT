//
// Created by H. Shinaoka on 2015/12/23.
//

#ifndef IMPSOLVER_MEASUREMENTS_H
#define IMPSOLVER_MEASUREMENTS_H

#include <complex>

#include "interaction_expansion.hpp"
#include "legendre.h"

template<class SOLVER_TYPE>
void evaluate_selfenergy_measurement_matsubara(const alps::accumulators::result_set &results,
                                               typename SOLVER_TYPE::matsubara_green_function_t &green_matsubara_measured,
                                               const typename SOLVER_TYPE::matsubara_green_function_t &bare_green_matsubara,
                                               std::vector<double>& densities,
                                               const double &beta, std::size_t n_site,
                                               std::size_t n_flavors, std::size_t n_matsubara, std::size_t n_matsubara_meas) {
  green_matsubara_measured.clear();
  std::cout << "evaluating self energy measurement: matsubara, reciprocal space" << std::endl;
  double sign = results["Sign"].template mean<double>();
  for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
    for (std::size_t site1 = 0; site1 < n_site; ++site1) {
      for (std::size_t site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Wk_real_name, Wk_imag_name;
        Wk_real_name << "Wk_real_" << flavor1 << "_" << site1 << "_" << site2;
        Wk_imag_name << "Wk_imag_" << flavor1 << "_" << site1 << "_" << site2;
        std::vector<double> mean_real = results[Wk_real_name.str().c_str()].template mean<std::vector<double> >();
        std::vector<double> mean_imag = results[Wk_imag_name.str().c_str()].template mean<std::vector<double> >();
        assert(mean_real.size()==n_matsubara_meas && mean_imag.size()==n_matsubara_meas);
        if (!(mean_real.size()==n_matsubara_meas && mean_imag.size()==n_matsubara_meas)) {
          throw std::logic_error("Logic error in evaluate_selfenergy_measurement_matsubara.");
        }
        for (unsigned int w = 0; w < n_matsubara_meas; ++w) {
          green_matsubara_measured(w, site1, site2, flavor1) = -std::complex<double>(mean_real[w], mean_imag[w]) / (sign*beta);
        }
      }
    }
  }
  for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
    for (std::size_t site1 = 0; site1 < n_site; ++site1) {
      for (std::size_t site2 = 0; site2 < n_site; ++site2) {
        for (std::size_t w = 0; w < n_matsubara; ++w) {
          green_matsubara_measured(w, site1, site2, flavor1) += bare_green_matsubara(w, site1, site2, flavor1);
        }
      }
    }
  }
  std::vector<double> dens = results["densities"].template mean<std::vector<double> >();
  for (std::size_t z=0; z<n_flavors; ++z)
    densities[z] = dens[z];

  //std::vector<double> tau = results["PerturbationOrderVertex"].template tau<std::vector<double> >();
  //std::cout << "Autocorrelation times of PerturbationOrderVertex : " << std::endl;
  //for (int iv=0; iv<tau.size(); ++iv) {
    //std::cout << " iv = " << iv << " tau = " << tau[iv] << std::endl;
  //}
}

template<class SOLVER_TYPE>
void evaluate_selfenergy_measurement_legendre(const typename alps::accumulators::result_set &results,
                                              typename SOLVER_TYPE::matsubara_green_function_t &green_matsubara_measured,
                                              typename SOLVER_TYPE::matsubara_green_function_t &sigma_green_matsubara_measured,
                                              typename SOLVER_TYPE::itime_green_function_t &green_itime_measured,
                                              const typename SOLVER_TYPE::matsubara_green_function_t &bare_green_matsubara,
                                              std::vector<double>& densities,
                                              const double &beta, std::size_t n_site,
                                              std::size_t n_flavors, std::size_t n_matsubara, std::size_t n_tau, std::size_t n_legendre) {
  green_matsubara_measured.clear();
  green_itime_measured.clear();
  std::cout << "evaluating self energy measurement: lengendre, real space" << std::endl;
  double sign = results["Sign"].template mean<double>();

  //Legendre expansion utils
  LegendreTransformer legendre_transformer(n_matsubara,n_legendre);
  auto& Tnl = legendre_transformer.Tnl();

  boost::multi_array<std::complex<double>,4> Sl(boost::extents[n_legendre][n_site][n_site][n_flavors]);
  boost::multi_array<std::complex<double>,4> Sw(boost::extents[n_matsubara][n_site][n_site][n_flavors]);

  //load S_l
  for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
    for (std::size_t site1 = 0; site1 < n_site; ++site1) {
      for (std::size_t site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Sl_real_name, Sl_imag_name;
        Sl_real_name << "Sl_real_" << flavor1 << "_" << site1 << "_" << site2;
        Sl_imag_name << "Sl_imag_" << flavor1 << "_" << site1 << "_" << site2;
        std::vector<double> mean_real = results[Sl_real_name.str().c_str()].template mean<std::vector<double> >();
        std::vector<double> mean_imag = results[Sl_imag_name.str().c_str()].template mean<std::vector<double> >();
        for (unsigned int i_l = 0; i_l < n_legendre; ++i_l) {
          Sl[i_l][site1][site2][flavor1] = std::complex<double>(mean_real[i_l], mean_imag[i_l])/sign;
        }
      }
    }
  }

  //compute S(iomega_n)
  {
    Eigen::MatrixXcd Sl_vec(n_legendre,1);
    for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
      for (std::size_t site1 = 0; site1 < n_site; ++site1) {
        for (std::size_t site2 = 0; site2 < n_site; ++site2) {
          for(std::size_t i_l = 0; i_l < n_legendre;  ++ i_l) {
            Sl_vec(i_l,0) = Sl[i_l][site1][site2][flavor1];
          }
          Eigen::MatrixXcd Sw_vec {Tnl * Sl_vec};
          for (std::size_t w = 0; w < n_matsubara; ++w) {
            Sw[w][site1][site2][flavor1] = Sw_vec(w,0);
          }
        }
      }
    }
  }

  //compute G(iomega_n) by Dyson eq.
  Eigen::MatrixXcd g0_mat(n_site, n_site, 0.0), s_mat(n_site, n_site, 0.0);//tmp arrays for gemm
  Eigen::MatrixXcd g_mat(n_site, n_site, 0.0);
  for (std::size_t w = 0; w < n_matsubara; ++w) {
    for (std::size_t z = 0; z < n_flavors; ++z) {

      //prepare a matrix for S and bare G
      for (std::size_t site1 = 0; site1 < n_site; ++site1) {
        for (std::size_t site2 = 0; site2 < n_site; ++site2) {
          s_mat(site1, site2) = Sw[w][site1][site2][z];
          g0_mat(site1, site2) = bare_green_matsubara(w, site1, site2, z);
          g_mat(site1, site2) = 0.0;//just for safety
        }
      }

      //solve Dyson equation
      g_mat = g0_mat + g0_mat * s_mat;

      //write back full G
      for (std::size_t site1 = 0; site1 < n_site; ++site1) {
        for (std::size_t site2 = 0; site2 < n_site; ++site2) {
          sigma_green_matsubara_measured(w, site1, site2, z) = s_mat(site1, site2);
          green_matsubara_measured(w, site1, site2, z) = g_mat(site1, site2);
        }
      }
    }
  }
}

template<class SOLVER_TYPE>
void compute_greens_functions(const typename alps::accumulators::result_set &results, const typename alps::params& parms, const std::string &output_file)
{
  std::cout<<"getting result!"<<std::endl;
  unsigned int n_matsubara = parms["N_MATSUBARA"];
  unsigned int n_matsubara_measurements = parms["NMATSUBARA_MEASUREMENTS"];
  unsigned int n_tau = parms["N_TAU"];
  spin_t n_flavors = parms["FLAVORS"];
  unsigned int n_site = parms["SITES"];
  double beta = parms["BETA"];
  int n_legendre = parms["N_LEGENDRE"];
  typename SOLVER_TYPE::itime_green_function_t green_itime_measured(n_tau+1, n_site, n_flavors);
  typename SOLVER_TYPE::matsubara_green_function_t green_matsubara_measured(n_matsubara, n_site, n_flavors);
  typename SOLVER_TYPE::matsubara_green_function_t sigma_green_matsubara_measured(n_matsubara, n_site, n_flavors);

  boost::shared_ptr<FourierTransformer> fourier_ptr;
  boost::shared_ptr<FourierTransformer> fourier_ptr_g0;
  FourierTransformer::generate_transformer_lowest_order(parms, fourier_ptr_g0);

  //find whether our data is in imaginary time or frequency:
  bool measure_in_matsubara=true;
  if(parms["HISTOGRAM_MEASUREMENT"]) {
    measure_in_matsubara=false;
  }

  std::vector<double> mean_order=results["PertOrder"].template mean<std::vector<double> >();

  std::cout<<"average matrix size was: "<<std::endl;
  std::ofstream matrix_size("matrix_size", std::ios::app);
  for(unsigned int i=0;i<n_flavors;++i){
    std::cout<<mean_order[i]<<"\t";
    matrix_size<<mean_order[i]<<"\t";
  }
  std::cout<<std::endl;
  matrix_size<<std::endl;
  std::cout<<"average sign was: "<<results["Sign"].template mean<double>()<<" error: "<<results["Sign"].template error<double>()<<std::endl;

  //Load G0 in Matsubara frequency
  typename SOLVER_TYPE::matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, n_flavors);
  typename SOLVER_TYPE::itime_green_function_t bare_green_itime(n_tau+1, n_site, n_flavors);//this is not used. just dummy
  std::vector<double> densities(n_flavors);
  boost::tie(bare_green_matsubara,bare_green_itime) =
      read_bare_green_functions<typename SOLVER_TYPE::COMPLEX_TYPE>(parms);

  //Matsubara freq. measurement
  if (n_matsubara_measurements>0) {
    evaluate_selfenergy_measurement_matsubara<SOLVER_TYPE>(results, green_matsubara_measured,
                                            bare_green_matsubara, densities,
                                            beta, n_site, n_flavors, n_matsubara, n_matsubara_measurements);

    //Fourier transformations
    FourierTransformer::generate_transformer_lowest_order(parms, fourier_ptr);
    fourier_ptr->append_tail(green_matsubara_measured, bare_green_matsubara, n_matsubara_measurements);
    fourier_ptr->backward_ft(green_itime_measured, green_matsubara_measured);

    alps::hdf5::archive ar(output_file, "w");
    green_matsubara_measured.write_hdf5(ar, "/G_omega");
    green_itime_measured.write_hdf5(ar, "/G_tau");
    bare_green_matsubara.write_hdf5(ar, "/G0_omega");
    bare_green_itime.write_hdf5(ar, "/G0_tau");
  }

  //Legendre measurement
  if (n_legendre>0) {
    typename SOLVER_TYPE::itime_green_function_t green_itime_measured_l(n_tau+1, n_site, n_flavors);
    typename SOLVER_TYPE::matsubara_green_function_t green_matsubara_measured_l(n_matsubara, n_site, n_flavors);
    typename SOLVER_TYPE::matsubara_green_function_t sigma_green_matsubara_measured_l(n_matsubara, n_site, n_flavors);

    evaluate_selfenergy_measurement_legendre<SOLVER_TYPE>(results, green_matsubara_measured_l, sigma_green_matsubara_measured_l, green_itime_measured_l,
                                             bare_green_matsubara, densities,
                                             beta, n_site, n_flavors, n_matsubara, n_tau, n_legendre);

    //just for test use: G(tau) should be computed directly from S_l to prevent truncation errors.
    //FourierTransformer::generate_transformer_lowest_order(parms, fourier_ptr);
    //fourier_ptr->append_tail(green_matsubara_measured_l, bare_green_matsubara, n_matsubara);
    //fourier_ptr->backward_ft(green_itime_measured_l, green_matsubara_measured_l);

    alps::hdf5::archive ar(output_file, "a");
    green_matsubara_measured_l.write_hdf5(ar, "/G_omega_legendre");
    sigma_green_matsubara_measured_l.write_hdf5(ar, "/SigmaG_omega_legendre");
    //green_itime_measured_l.write_hdf5(ar, "/G_tau_legendre");
  }
}

#endif //IMPSOLVER_MEASUREMENTS_H
