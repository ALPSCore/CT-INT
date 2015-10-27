/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "interaction_expansion.hpp"
#include <ctime>
#include "xml.h"
#include "alps/ngs/make_deprecated_parameters.hpp"

#include <boost/lexical_cast.hpp>
#include "boost/tuple/tuple.hpp"

//global variables

frequency_t c_or_cdagger::nm_;
bool c_or_cdagger::use_static_exp_;
unsigned int c_or_cdagger::ntau_;
double c_or_cdagger::beta_;
double *c_or_cdagger::omegan_;
std::complex<double> *c_or_cdagger::exp_iomegan_tau_;


InteractionExpansion::InteractionExpansion(const alps::params &parms, int node)
: alps::mcbase(parms,node),
max_order(parms["MAX_ORDER"] | 2048),
n_flavors(parms["FLAVORS"] | (parms["N_ORBITALS"] | 2)),
n_site(parms["SITES"] | 1),
n_matsubara((int)(parms["NMATSUBARA"]|parms["N_MATSUBARA"])),
n_matsubara_measurements(parms["NMATSUBARA_MEASUREMENTS"] | (int)n_matsubara),
n_tau((int)(parms["N"]|parms["N_TAU"])),
n_tau_inv(1./n_tau),
n_self(parms["NSELF"] | (int)(10*n_tau)),
n_legendre(parms["N_LEGENDRE"] | 0),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned int)parms["THERMALIZATION"]),        
max_time_in_seconds(parms["MAX_TIME"] | 86400),
beta((double)parms["BETA"]),                        
temperature(1./beta),
onsite_U((double)parms["U"]),                        
alpha((double)parms["ALPHA"]),
U(alps::make_deprecated_parameters(parms)),
Uijkl(alps::make_deprecated_parameters(parms)),
recalc_period(parms["RECALC_PERIOD"] | 5000),
//measurement_period(parms["MEASUREMENT_PERIOD"] | (parms["N_MEAS"] | 200)),
measurement_period(parms["MEASUREMENT_PERIOD"] | 500*n_flavors*n_site),
convergence_check_period(parms["CONVERGENCE_CHECK_PERIOD"] | (int)recalc_period),
almost_zero(parms["ALMOSTZERO"] | 1.e-16),
seed(parms["SEED"] | 0),
green_matsubara(n_matsubara, n_site, n_flavors),
bare_green_matsubara(n_matsubara,n_site, n_flavors), 
bare_green_itime(n_tau+1, n_site, n_flavors),
green_itime(n_tau+1, n_site, n_flavors),
pert_hist(max_order),
legendre_transformer(n_matsubara,n_legendre)
{
  //initialize measurement method
  if (parms["HISTOGRAM_MEASUREMENT"] | false) {
    measurement_method=selfenergy_measurement_itime_rs;
    //std::cout << "debug: measure_in_itime_rs " << std::endl;
  } else {
    measurement_method=selfenergy_measurement_matsubara;
    //std::cout << "debug: measure_matsubara " << std::endl;
  }
  for(unsigned int i=0;i<n_flavors;++i)
    g0.push_back(green_matrix(n_tau, 20));
  //other parameters
  weight=0;
  sign=1;
  step=0;
  start_time=time(NULL);
  measurement_time=0;
  update_time=0;
  thermalized=therm_steps==0?true:false;


  //load bare Green's function
  boost::tie(bare_green_matsubara,bare_green_itime) = read_bare_green_functions<double>(parms);//G(tau) is assume to be real.

  //FourierTransformer::generate_transformer(alps::make_deprecated_parameters(parms), fourier_ptr);
  //fourier_ptr->backward_ft(bare_green_itime, bare_green_matsubara);

  //initialize the simulation variables
  initialize_simulation(parms);
  if(node==0) {print(std::cout);}
  vertex_histograms=new simple_hist *[n_flavors*n_flavors];
  vertex_histogram_size=100;
  for(unsigned int i=0;i<n_flavors*n_flavors;++i){
    vertex_histograms[i]=new simple_hist(vertex_histogram_size);
  }
  c_or_cdagger::initialize_simulation(parms);
  
  //if(n_site !=1) throw std::invalid_argument("you're trying to run this code for more than one site. Do you know what you're doing?!?");
}



void InteractionExpansion::update()
{
  boost::timer::cpu_timer timer;
  //std::cout << "step " << step << std::endl;
  for(std::size_t i=0;i<measurement_period;++i){
    step++;
    interaction_expansion_step();                
    if(vertices.size()<max_order)
      pert_hist[vertices.size()]++;
    if(step % recalc_period ==0) {
      //just for debug
      sanity_check();
      reset_perturbation_series();
    }
  }
  measurements["UpdateTimeMsec"] << timer.elapsed().wall*1E-6;
}

void InteractionExpansion::measure(){
  std::valarray<double> timings(2);
  measure_observables(timings);
  measurements["MeasurementTimeMsec"] << timings;
}



double InteractionExpansion::fraction_completed() const{
  //check for error convergence
  //std::cout << "fraction " << ((step-therm_steps) / (double) mc_steps) << std::endl;
  //std::cout << "debug fraction " << "step=" << step << " therm_steps " << therm_steps << " mc_steps= " << mc_steps << std::endl;
  if (!thermalized) {
    return 0.;
  }
  if(time(NULL)-start_time> max_time_in_seconds){
    std::cout<<"we ran out of time!"<<std::endl;
    return 1;
  }
  assert(step>=therm_steps);
  return ((step-therm_steps) / (double) mc_steps);
}



///do all the setup that has to be done before running the simulation.
void InteractionExpansion::initialize_simulation(const alps::params &parms)
{
  weight=0;
  sign=1;
  //set the right dimensions:
  for(spin_t flavor=0;flavor<n_flavors;++flavor)
    M.push_back(inverse_m_matrix());
  vertices.clear();
  pert_hist.clear();
  //initialize ALPS observables
  initialize_observables();
  green_matsubara=bare_green_matsubara;
  green_itime=bare_green_itime;
}


void InteractionExpansion::sanity_check() {
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    alps::numeric::matrix<GTYPE> tmp(M[flavor].matrix()), G0(M[flavor].matrix());
    std::fill(tmp.get_values().begin(), tmp.get_values().end(), 0);
    std::fill(G0.get_values().begin(), G0.get_values().end(), 0);

    const size_t Nv = M[flavor].matrix().num_rows();
    for (size_t q=0; q<Nv; ++q) {
      for (size_t p=0; p<Nv; ++p) {
        G0(p, q) = green0_spline_new(M[flavor].annihilators()[p], M[flavor].creators()[q]);
      }
    }
    for (size_t p=0; p<Nv; ++p) {
      G0(p, p) += M[flavor].alpha()[p];
    }

    gemm(G0, M[flavor].matrix(), tmp);
    bool OK = true;
    for (size_t q=0; q<Nv; ++q) {
      for (size_t p=0; p<Nv; ++p) {
        if (p==q) {
          OK = OK && (std::abs(tmp(p,q)-1.)<1E-8);
        } else {
          OK = OK && (std::abs(tmp(p,q))<1E-8);
        }
      }
    }
    if (!OK) {
      throw std::runtime_error("There is something wrong: G^{-1} != M.");
    }
    //std::cout << " OK " << OK << std::endl;
  }
}


void c_or_cdagger::initialize_simulation(const alps::params &p)
{
  beta_=p["BETA"];
  nm_=p["NMATSUBARA_MEASUREMENTS"] | (p["NMATSUBARA"]|p["N_MATSUBARA"]);
  omegan_ = new double[nm_];
  for(unsigned int i=0;i<nm_;++i) {
    omegan_[i]=(2.*i+1.)*M_PI/beta_;
  }
  if(p.defined("TAU_DISCRETIZATION_FOR_EXP")) {
    ntau_=p["TAU_DISCRETIZATION_FOR_EXP"];
    use_static_exp_=true;
    exp_iomegan_tau_=new std::complex<double> [2*nm_*ntau_];
    if(exp_iomegan_tau_==0){throw std::runtime_error("not enough memory for computing exp!"); }
    std::cout<<"starting computation of exp values for measurement"<<std::endl;
    for(unsigned int i=0;i<ntau_;++i){
      double tau=i*beta_/(double)ntau_;
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[2*nm_*i + o] = std::complex<double>(cos(omegan_[o]*tau), sin(omegan_[o]*tau));
      for(unsigned int o=0;o<nm_;++o)
        exp_iomegan_tau_[2*nm_*i + nm_ + o] = std::complex<double>(cos(omegan_[o]*tau), -sin(omegan_[o]*tau));
    }
    std::cout<<"done exp computation."<<std::endl;
  } else {
    use_static_exp_=false;
  }
}

