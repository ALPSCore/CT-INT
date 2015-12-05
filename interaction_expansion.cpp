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
node(node),
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
//onsite_U((double)parms["U"]),
//alpha((double)parms["ALPHA"]),
//U(alps::make_deprecated_parameters(parms)),
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
legendre_transformer(n_matsubara,n_legendre),
n_multi_vertex_update(parms["N_MULTI_VERTEX_UPDATE"] | 1),
statistics_ins((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, n_multi_vertex_update-1),
statistics_rem((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, n_multi_vertex_update-1),
statistics_shift((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, 1),
statistics_dv_rem(0, 0, 0),
statistics_dv_ins(0, 0, 0),
simple_statistics_ins(n_multi_vertex_update),
simple_statistics_rem(n_multi_vertex_update),
is_thermalized_in_previous_step_(false),
window_width(parms.defined("WINDOW_WIDTH") ? beta*static_cast<double>(parms["WINDOW_WIDTH"]) : 1000.0*beta),
window_dist(boost::random::exponential_distribution<>(1/window_width)),
add_helper(n_flavors),
remove_helper(n_flavors),
shift_helper(n_flavors, parms.defined("SHIFT_WINDOW_WIDTH") ? beta*static_cast<double>(parms["SHIFT_WINDOW_WIDTH"]) : 1000.0*beta),
n_ins_rem(parms["N_INS_REM_VERTEX"] | 1),
n_shift(parms["N_SHIFT_VERTEX"] | 0),
force_quantum_number_conservation(parms.defined("FORCE_QUANTUM_NUMBER_CONSERVATION") ? parms["FORCE_QUANTUM_NUMBER_CONSERVATION"] : false),
//force_quantum_number_within_range(parms.defined("FORCE_QUANTUM_NUMBER_WITHIN_RANGE") ? parms["FORCE_QUANTUM_NUMBER_WITHIN_RANGE"] : false),
alpha_scale_min(1),
alpha_scale_max(parms["ALPHA_SCALE_MAX"] | 1),
alpha_scale(1.),
alpha_scale_max_meas(parms["ALPHA_SCALE_MEASURE_MAX"] | 1),
alpha_scale_update_period(parms["ALPHA_SCALE_UPDATE_PERIOD"] | -1),
flat_histogram_alpha(std::log(alpha_scale_max), std::log(alpha_scale_min), 100),
single_vertex_update_non_density_type(parms.defined("SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE") ? parms["SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE"] : true)
{
  //initialize measurement method
  if (parms["HISTOGRAM_MEASUREMENT"] | false) {
    measurement_method=selfenergy_measurement_itime_rs;
  } else {
    measurement_method=selfenergy_measurement_matsubara;
  }
  for(unsigned int i=0;i<n_flavors;++i)
    g0.push_back(green_matrix(n_tau, 20));

  //other parameters
  weight=0;
  sign=1;
  det=1;
  step=0;
  start_time=time(NULL);
  measurement_time=0;
  update_time=0;

  //load bare Green's function
  boost::tie(bare_green_matsubara,bare_green_itime) = read_bare_green_functions<double>(parms);//G(tau) is assume to be real.

  //make quantum numbers
  if (n_multi_vertex_update>1) {
    quantum_number_vertices = make_quantum_numbers(bare_green_itime, Uijkl.get_vertices(), groups, group_map, almost_zero);
    qn_dim = quantum_number_vertices[0][0].size();
    group_dim.clear(); group_dim.resize(qn_dim, 0);
    const int qn_dim_f = qn_dim/n_flavors;
    for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
      for (int g=0; g<groups[flavor].size(); ++g) {
        group_dim[g+flavor*qn_dim_f] = groups[flavor][g].size();
      }
    }

    //for double vertex update
    find_valid_pair_multi_vertex_update(Uijkl.get_vertices(), quantum_number_vertices, mv_update_valid_pair, mv_update_valid_pair_flag);
    if (mv_update_valid_pair.size()==0)
      throw std::runtime_error("No valid vertex pair for double vertex update. Please deactivate double vertex update.");
  }

  is_density_density_type.resize(Uijkl.n_vertex_type());
  for (int iv=0; iv<Uijkl.n_vertex_type(); ++iv) {
    is_density_density_type[iv] = Uijkl.get_vertex(iv).is_density_type();
  }

  //occ changes
  if (n_multi_vertex_update>1) {
    for (int iv=0; iv<Uijkl.n_vertex_type(); ++iv) {
      Uijkl.get_vertex(iv).make_quantum_numbers(group_map, qn_dim/n_flavors);
    }
  }

  if (n_multi_vertex_update>1) {
    symm_exp_dist = SymmExpDist(parms["DOUBLE_VERTEX_UPDATE_A"], parms["DOUBLE_VERTEX_UPDATE_B"], beta);
    statistics_dv_ins = scalar_histogram_flavors((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, mv_update_valid_pair.size());
    statistics_dv_rem = scalar_histogram_flavors((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, mv_update_valid_pair.size());
  }

  //set up parameters for updates
  std::vector<double> proposal_prob(n_multi_vertex_update, 1.0);
  if (params.defined("MULTI_VERTEX_UPDATE_PROPOSAL_RATE")) {
    proposal_prob.resize(0);
    std::stringstream ss(params["MULTI_VERTEX_UPDATE_PROPOSAL_RATE"].cast<std::string>());
    double rtmp;
    while (ss >> rtmp) {
      proposal_prob.push_back(rtmp);
    }
    if (proposal_prob.size()!=n_multi_vertex_update)
      throw std::runtime_error("The number of elements in MULTI_VERTEX_UPDATE_PROPOSAL_RATE is different from N_MULTI_VERTEX_UPDATE");
  }
  update_prop = update_proposer(n_multi_vertex_update, proposal_prob);

  //initialize the simulation variables
  initialize_simulation(parms);

  if(node==0 && n_multi_vertex_update>1) {
    print(std::cout);

    std::cout << std::endl << "Analysis of quantum numbers"  << std::endl;
    for (int flavor=0; flavor<n_flavors; ++flavor) {
      std::cout << "  Flavor " << flavor << " has " << groups[flavor].size() << " group(s)." << std::endl;
      print_group(groups[flavor]);
    }
    std::cout << std::endl;
    std::cout << std::endl;

    if (mv_update_valid_pair.size()>0) {
      std::cout << std::endl << "Vertex pairs for double vertex update." << std::endl;
      for (int i=0; i<mv_update_valid_pair.size(); ++i)
        std::cout << " type " << mv_update_valid_pair[i].first << ", type " << mv_update_valid_pair[i].second << std::endl;
    } else {
      std::cout << std::endl << "No vertex pairs for double vertex update." << std::endl;
    }
  }
  vertex_histograms=new simple_hist *[n_flavors];
  vertex_histogram_size=100;
  for(unsigned int i=0;i<n_flavors;++i){
    vertex_histograms[i]=new simple_hist(vertex_histogram_size);
  }
  c_or_cdagger::initialize_simulation(parms);
}



void InteractionExpansion::update()
{
  //std::cout << " update called  node " << node << std::endl;
  std::valarray<double> t_meas(0.0, 2);
  const bool is_thermalized_now = is_thermalized();
  if (!is_thermalized_in_previous_step_ && is_thermalized_now) {
    prepare_for_measurement();
  }
  is_thermalized_in_previous_step_ = is_thermalized_now;

  //if (node==0)
  //std::cout << " node = " << node << " step = " << step << " random= " << random() << std::endl;
  //std::cout << " node = " << node << " step = " << step << std::endl;

  itime_vertex_container itime_vertices_bak = itime_vertices;
  big_inverse_m_matrix M_bak(M);

  for(std::size_t i=0;i<measurement_period;++i){
    step++;
    boost::timer::cpu_timer timer;

    //alpha
    flat_histogram_alpha.measure(std::log(alpha_scale));
    if (flat_histogram_alpha.count()>1000) {
      flat_histogram_alpha.update_dos();
    }

    try{
      for (int i_ins_rem=0; i_ins_rem<n_ins_rem; ++i_ins_rem)
        removal_insertion_update();

      double t_m = timer.elapsed().wall;
  
      for (int i_shift=0; i_shift<n_shift; ++i_shift)
        shift_update();
  
      t_meas[0] += t_m;
      t_meas[1] += (timer.elapsed().wall-t_m);
      if(itime_vertices.size()<max_order)
        pert_hist[itime_vertices.size()]++;

      //update alpha_scale
      if(alpha_scale_update_period>0 && step%alpha_scale_update_period==0) {
        alpha_update();
      }
    } catch (std::runtime_error e) {
      std::cout << " Runtime error at rank = " << node << " step = " << step << " : " << e.what() << ". This may be because we encountered a singular matrix." << std::endl;
      itime_vertices = itime_vertices_bak;
      M = M_bak;
      reset_perturbation_series(false);
    }

    if(step % recalc_period ==0) {
      reset_perturbation_series(true);
    }
  }


  t_meas *= 1E-6/measurement_period;
  measurements["UpdateTimeMsec"] << t_meas;
}

void InteractionExpansion::measure(){
  if (alpha_scale>alpha_scale_max_meas)
    return;

  //In the below, real physical quantities are measured.
  //if (n_multi_vertex_update>1) {
    //bool flag = is_quantum_number_conserved(itime_vertices);
    //measurements["QuantumNumberConserved"] << (int) flag;
    //if(!flag) return;
  //}

  std::valarray<double> timings(2);
  measure_observables(timings);
  measurements["MeasurementTimeMsec"] << timings;
}



double InteractionExpansion::fraction_completed() const{
  if (!is_thermalized()) {
    return 0.;
  } else {
    return ((step - therm_steps) / (double) mc_steps);
  }
}



///do all the setup that has to be done before running the simulation.
void InteractionExpansion::initialize_simulation(const alps::params &parms)
{
  weight=0;
  sign=1;
  //set the right dimensions:
  for(spin_t flavor=0;flavor<n_flavors;++flavor)
    M.push_back(inverse_m_matrix());
  pert_hist.clear();
  //initialize ALPS observables
  initialize_observables();
  green_matsubara=bare_green_matsubara;
  green_itime=bare_green_itime;
}

bool InteractionExpansion::is_quantum_number_conserved(const itime_vertex_container& vertices) {
  //using namespace boost::lambda;

  const int Nv = vertices.size();

  if (Nv==0)
    return true;

  std::valarray<int> qn_t(0, qn_dim), qn_max(0, qn_dim), qn_min(0, qn_dim);
  itime_vertex_container vertices_sorted(vertices);//sort vertices in decreasing order (in time)
  std::sort(vertices_sorted.begin(), vertices_sorted.end());

  for (int iv=0; iv<Nv; ++iv) {
    const vertex_definition<GTYPE> vd = Uijkl.get_vertex(vertices_sorted[iv].type());
    vd.apply_occ_change(vertices_sorted[iv].af_state(), qn_t, qn_max, qn_min);
  }

  //check if the quantum number is conserved
  for (int i=0; i<qn_t.size(); ++i) {
    if (qn_t[i]!=0) {
      return false;
    }
  }

  return true;
}


//This makes sence in the absence of a bath
bool InteractionExpansion::is_quantum_number_within_range(const itime_vertex_container& vertices) {
  //using namespace boost::lambda;

  const int Nv = vertices.size();

  if (Nv==0)
    return true;

  std::valarray<int> qn_t(0, qn_dim), qn_max(0, qn_dim), qn_min(0, qn_dim);
  itime_vertex_container vertices_sorted(vertices);//sort vertices in decreasing order (in time)
  std::sort(vertices_sorted.begin(), vertices_sorted.end());

  for (int iv=0; iv<Nv; ++iv) {
    const vertex_definition<GTYPE> vd = Uijkl.get_vertex(vertices_sorted[iv].type());
    vd.apply_occ_change(vertices_sorted[iv].af_state(), qn_t, qn_max, qn_min);
  }

  //check if the quantum number is within range
  for (int iq=0; iq<qn_dim; ++iq) {
    //note: group_dim[iq] is zero for non-existing group
    if (qn_max[iq]-qn_min[iq]>group_dim[iq]) {
      return false;
    }
  }

  return true;
}

void InteractionExpansion::sanity_check() {
#ifndef NDEBUG
  M.sanity_check(itime_vertices);

  //recompute M
  double sign_exact = 1;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (num_rows(M[flavor].matrix())==0) {
      continue;
    }
    alps::numeric::matrix<GTYPE> G0(M[flavor].matrix());
    std::fill(G0.get_values().begin(), G0.get_values().end(), 0);

    const size_t Nv = M[flavor].matrix().num_rows();
    for (size_t q=0; q<Nv; ++q) {
      for (size_t p=0; p<Nv; ++p) {
        G0(p, q) = green0_spline_for_M(flavor, p, q);
      }
    }
    for (size_t p=0; p<Nv; ++p) {
      G0(p, p) -= M[flavor].alpha_at(p);
    }

    sign_exact *= boost::math::sign(alps::numeric::determinant(G0));

    alps::numeric::matrix<GTYPE> tmp = mygemm(G0, M[flavor].matrix());
    bool OK = true;
    double max_diff = 0;
    for (size_t q=0; q<Nv; ++q) {
      for (size_t p=0; p<Nv; ++p) {
        bool flag;
        if (p==q) {
          max_diff = std::max(max_diff, std::abs(tmp(p,q)-1.));
          flag = (std::abs(tmp(p,q)-1.)<1E-5);
        } else {
          max_diff = std::max(max_diff, std::abs(tmp(p,q)));
          flag = (std::abs(tmp(p,q))<1E-5);
        }
        OK = OK && flag;
        if(!flag) {
          std::cout << " p, q = " << p << " " << q << " " << tmp(p,q) << std::endl;
        }
      }
    }
    std::cout << "sanity check max_diff " << max_diff << std::endl;
  }

  //contribution from (-U)^n
  {
    const std::vector<vertex_definition<GTYPE> >& vd = Uijkl.get_vertices();
    for (int iv=0; iv<itime_vertices.size(); ++iv)
      sign_exact *= boost::math::sign(-vd[itime_vertices[iv].type()].Uval());
  }
  assert(sign==sign_exact);

  //check itime_vertex list
  for (itime_vertex_container::iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    int type = it->type();
    assert(Uijkl.get_vertices()[type].is_density_type()==it->is_density_type());
  }
#endif
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

void HubbardInteractionExpansion::prepare_for_measurement()
{
  //update_prop.finish_learning((node==0));
  //update_prop.finish_learning(true);//REMOVE AFTER DEBUG
  statistics_ins.reset();
  statistics_rem.reset();
  statistics_shift.reset();
  simple_statistics_ins.reset();
  simple_statistics_rem.reset();
  statistics_dv_ins.reset();
  statistics_dv_rem.reset();
}
