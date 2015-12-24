//
// Created by H. Shinaoka on 2015/12/23.
//

#ifndef IMPSOLVER_MEASUREMENTS_H
#define IMPSOLVER_MEASUREMENTS_H

#include "interaction_expansion.hpp"
#include <complex>
#include <alps/alea.h>
#include "alps/ngs/make_deprecated_parameters.hpp"
#include "legendre.h"

template<class SOLVER_TYPE>
void evaluate_selfenergy_measurement_matsubara(const typename alps::results_type<SOLVER_TYPE>::type &results,
                                               matsubara_green_function_t &green_matsubara_measured,
                                               const matsubara_green_function_t &bare_green_matsubara,
                                               std::vector<double>& densities,
                                               const double &beta, std::size_t n_site,
                                               std::size_t n_flavors, std::size_t n_matsubara, std::size_t n_matsubara_measurements);

template<class SOLVER_TYPE>
void evaluate_selfenergy_measurement_itime_rs(const typename alps::results_type<SOLVER_TYPE>::type &results,
                                              itime_green_function_t &green_result,
                                              const itime_green_function_t &green0,
                                              const double &beta, const int n_site,
                                              const int n_flavors, const int n_tau, const int n_self);

template<class SOLVER_TYPE>
void evaluate_selfenergy_measurement_legendre(const typename alps::results_type<SOLVER_TYPE>::type &results,
                                              matsubara_green_function_t &green_matsubara_measured,
                                              matsubara_green_function_t &sigma_green_matsubara_measured,
                                              itime_green_function_t &green_itime_measured,
                                              const matsubara_green_function_t &bare_green_matsubara,
                                              std::vector<double>& densities,
                                              const double &beta, std::size_t n_site,
                                              std::size_t n_flavors, std::size_t n_matsubara, std::size_t n_tau, std::size_t n_legendre);


template<class SOLVER_TYPE>
void compute_greens_functions(const typename alps::results_type<SOLVER_TYPE>::type &results, const typename alps::parameters_type<SOLVER_TYPE>::type& parms, const std::string &output_file)
{
  std::cout<<"getting result!"<<std::endl;
  unsigned int n_matsubara = parms["NMATSUBARA"]|parms["N_MATSUBARA"];
  unsigned int n_matsubara_measurements=parms["NMATSUBARA_MEASUREMENTS"] | (int)n_matsubara;
  unsigned int n_tau=parms["N"]|parms["N_TAU"];
  unsigned int n_self=parms["NSELF"] | (int)(10*n_tau);
  spin_t n_flavors(parms["FLAVORS"] | (parms["N_ORBITALS"]| 2));
  unsigned int n_site(parms["SITES"] | 1);
  double beta(parms["BETA"]);
  const int n_legendre(parms["N_LEGENDRE"] | 0);
  itime_green_function_t green_itime_measured(n_tau+1, n_site, n_flavors);
  matsubara_green_function_t green_matsubara_measured(n_matsubara, n_site, n_flavors);
  matsubara_green_function_t sigma_green_matsubara_measured(n_matsubara, n_site, n_flavors);

  boost::shared_ptr<FourierTransformer> fourier_ptr;
  boost::shared_ptr<FourierTransformer> fourier_ptr_g0;
  //FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr_g0);
  FourierTransformer::generate_transformer_lowest_order(alps::make_deprecated_parameters(parms), fourier_ptr_g0);
  //find whether our data is in imaginary time or frequency:
  bool measure_in_matsubara=true;
  if(parms["HISTOGRAM_MEASUREMENT"] | false) {
    measure_in_matsubara=false;
  }
  std::cout << "debug: measure_in_matsubara " << measure_in_matsubara << std::endl;
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
  //single particle Green function measurements
  matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, n_flavors);
  itime_green_function_t bare_green_itime(n_tau+1, n_site, n_flavors);//this is not used. just dummy
  std::vector<double> densities(n_flavors);

  //Load G0 in Matsubara frequency
  boost::tie(bare_green_matsubara,bare_green_itime) = read_bare_green_functions<double>(parms);


  if(measure_in_matsubara) {
    evaluate_selfenergy_measurement_matsubara<SOLVER_TYPE>(results, green_matsubara_measured,
                                              bare_green_matsubara, densities,
                                              beta, n_site, n_flavors, n_matsubara, n_matsubara_measurements);
  }
  else {
    throw std::runtime_error("Please make sure that sign is correcly treated!");
    //fourier_ptr_g0->backward_ft(bare_green_itime, bare_green_matsubara);
    //evaluate_selfenergy_measurement_itime_rs(results, green_itime_measured, bare_green_itime,
    //beta, n_site, n_flavors, n_tau, n_self);
  }
  //Fourier transformations
  FourierTransformer::generate_transformer_lowest_order(alps::make_deprecated_parameters(parms), fourier_ptr);
  if (measure_in_matsubara) {
    fourier_ptr->append_tail(green_matsubara_measured, bare_green_matsubara, n_matsubara_measurements);
    fourier_ptr->backward_ft(green_itime_measured, green_matsubara_measured);
  } else {
    throw std::runtime_error("Please make sure that sign is correcly treated!");
    //fourier_ptr->forward_ft(green_itime_measured, green_matsubara_measured);
  }


  {
    //alps::hdf5::archive ar(output_file, "a");
    alps::hdf5::archive ar(output_file, "w");
    green_matsubara_measured.write_hdf5(ar, "/G_omega");
    green_itime_measured.write_hdf5(ar, "/G_tau");
    bare_green_matsubara.write_hdf5(ar, "/G0_omega");
    bare_green_itime.write_hdf5(ar, "/G0_tau");
  }

  //Legendre measurement
  if (n_legendre>0) {
    itime_green_function_t green_itime_measured_l(n_tau+1, n_site, n_flavors);
    matsubara_green_function_t green_matsubara_measured_l(n_matsubara, n_site, n_flavors);
    matsubara_green_function_t sigma_green_matsubara_measured_l(n_matsubara, n_site, n_flavors);

    evaluate_selfenergy_measurement_legendre<SOLVER_TYPE>(results, green_matsubara_measured_l, sigma_green_matsubara_measured_l, green_itime_measured_l,
                                             bare_green_matsubara, densities,
                                             beta, n_site, n_flavors, n_matsubara, n_tau, n_legendre);

    //just for test use: G(tau) should be computed directly from S_l to prevent truncation errors.
    FourierTransformer::generate_transformer_lowest_order(alps::make_deprecated_parameters(parms), fourier_ptr);
    fourier_ptr->append_tail(green_matsubara_measured_l, bare_green_matsubara, n_matsubara);
    fourier_ptr->backward_ft(green_itime_measured_l, green_matsubara_measured_l);

    alps::hdf5::archive ar(output_file, "a");
    green_matsubara_measured_l.write_hdf5(ar, "/G_omega_legendre");
    sigma_green_matsubara_measured_l.write_hdf5(ar, "/SigmaG_omega_legendre");
    green_itime_measured_l.write_hdf5(ar, "/G_tau_legendre");
  }
}

#endif //IMPSOLVER_MEASUREMENTS_H
