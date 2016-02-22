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
 *
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

#include <ctime>

#include <boost/lexical_cast.hpp>
#include "boost/tuple/tuple.hpp"

#include "alps/ngs/make_deprecated_parameters.hpp"

#include "interaction_expansion.hpp"

template<typename  T>
BareGreenInterpolate<T>::BareGreenInterpolate(const alps::params &p) :
    beta_((double)p["BETA"]),
    temp_(1.0/beta_),
    ntau_((int)(p["N"]|p["N_TAU"])),
    n_flavors_(p["FLAVORS"] | (p["N_ORBITALS"] | 2)),
    n_sites_(p["SITES"] | 1),
    dbeta_(beta_/ntau_)
{
  green_function<std::complex<double> > bare_green_matsubara(ntau_,  n_sites_, n_flavors_),
      bare_green_itime(ntau_+1, n_sites_, n_flavors_);

  boost::tie(bare_green_matsubara,bare_green_itime) =
      read_bare_green_functions<std::complex<double> >(p);

  assert(ntau_==bare_green_itime.ntime()-1);
  AB_.resize(boost::extents[n_flavors_][n_sites_][n_sites_][ntau_+1]);

  for (int flavor=0; flavor<n_flavors_; ++flavor) {
    for (int site1=0; site1<n_sites_; ++site1) {
      for (int site2=0; site2<n_sites_; ++site2) {
        for (int tau=0; tau<ntau_; ++tau) {
          const T a =
              mycast<T>(
                  (bare_green_itime(tau+1,site1,site2,flavor)-bare_green_itime(tau,site1,site2,flavor))/dbeta_
              );
          const T b = mycast<T>(bare_green_itime(tau,site1,site2,flavor));

          AB_[flavor][site1][site2][tau] = std::make_pair(a,b);
        }
        AB_[flavor][site1][site2][ntau_] = std::make_pair(0.0,mycast<T>(bare_green_itime(ntau_,site1,site2,flavor)));
      }
    }
  }
}

/*
 * -beta <= d tau < beta
 */
template<typename  T>
T BareGreenInterpolate<T>::operator()(const annihilator& c, const creator& cdagger) const {
  assert(c.flavor()==cdagger.flavor());

  const int flavor = c.flavor();

  const int site1=c.s(), site2=cdagger.s();
  double dt = c.t().time()-cdagger.t().time();
  if (dt==0.0) {
    if (c.t().small_index() > cdagger.t().small_index()) { //G(+delta)
      return AB_[flavor][site1][site2][0].second;
    } else { //G(-delta)
      return -AB_[flavor][site1][site2][ntau_].second;
    }
  } else {
    T coeff = 1.0;
    while (dt>=beta_) {
      dt -= beta_;
      coeff *= -1.0;
    }
    while (dt<0.0) {
      dt += beta_;
      coeff *= -1.0;
    }

    assert(dt>=0 && dt<=beta_);
    const int time_index_1 = (int)(dt*ntau_*temp_);
    return coeff*(AB_[flavor][site1][site2][time_index_1].first*(dt-time_index_1*dbeta_) + AB_[flavor][site1][site2][time_index_1].second);
  }
}

//if delta_t==0, we assume delta_t = +0
template<typename  T>
T BareGreenInterpolate<T>::operator()(double delta_t, int flavor, int site1, int site2) const {
  double dt = delta_t;
  T coeff = 1.0;
  while (dt>=beta_) {
    dt -= beta_;
    coeff *= -1.0;
  }
  while (dt<0.0) {
    dt += beta_;
    coeff *= -1.0;
  }

  if (dt==0.0) dt += 1E-8;

  const int time_index_1 = (int)(dt*ntau_*temp_);
  return coeff*(AB_[flavor][site1][site2][time_index_1].first*(dt-time_index_1*dbeta_) + AB_[flavor][site1][site2][time_index_1].second);
}

template<typename  T>
bool BareGreenInterpolate<T>::is_zero(int site1, int site2, int flavor, double eps) const {
  return std::abs(operator()(beta_*1E-5, flavor, site1, site2))<eps && std::abs(operator()(beta_*(1-1E-5), flavor, site1, site2))<eps;
}

template<class TYPES>
InteractionExpansion<TYPES>::InteractionExpansion(const alps::params &parms, int node, const boost::mpi::communicator& communicator)
: InteractionExpansionBase(parms,node,communicator),
node(node),
max_order(parms["MAX_ORDER"] | 2048),
n_flavors(parms["FLAVORS"]),
n_site(parms["SITES"]),
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
Uijkl(alps::make_deprecated_parameters(parms)),
num_U_scale(parms["NUM_U_SCALE"] | 1),
min_U_scale(parms["MIN_U_SCALE"] | 1.0),
U_scale_vals(num_U_scale),
U_scale_index(num_U_scale),
M_flavors(n_flavors),
recalc_period(parms["RECALC_PERIOD"] | 5000),
measurement_period(parms["MEASUREMENT_PERIOD"] | 500*n_flavors*n_site),
convergence_check_period(parms["CONVERGENCE_CHECK_PERIOD"] | (int)recalc_period),
almost_zero(parms["ALMOSTZERO"] | 1.e-16),
seed(parms["SEED"] | 0),
bare_green_matsubara(n_matsubara,n_site, n_flavors),
bare_green_itime(n_tau+1, n_site, n_flavors),
pert_hist(max_order),
legendre_transformer(n_matsubara,n_legendre),
//n_multi_vertex_update(parms["N_MULTI_VERTEX_UPDATE"] | 1),
//statistics_ins((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, n_multi_vertex_update-1),
//statistics_rem((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, n_multi_vertex_update-1),
//statistics_shift((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, 1),
//statistics_dv_rem(0, 0, 0),
//statistics_dv_ins(0, 0, 0),
//simple_statistics_ins(n_multi_vertex_update),
//simple_statistics_rem(n_multi_vertex_update),
is_thermalized_in_previous_step_(false),
//window_width(parms.defined("WINDOW_WIDTH") ? beta*static_cast<double>(parms["WINDOW_WIDTH"]) : 1000.0*beta),
//window_dist(boost::random::exponential_distribution<>(1/window_width)),
//shift_helper(n_flavors, parms.defined("SHIFT_WINDOW_WIDTH") ? beta*static_cast<double>(parms["SHIFT_WINDOW_WIDTH"]) : 1000.0*beta),
n_ins_rem(parms["N_INS_REM_VERTEX"] | 1),
n_shift(parms["N_VERTEX_SHIFT"] | 1),
n_spin_flip(parms["N_SPIN_FLIP"] | 1),
force_quantum_number_conservation(parms.defined("FORCE_QUANTUM_NUMBER_CONSERVATION") ? parms["FORCE_QUANTUM_NUMBER_CONSERVATION"] : false),
single_vertex_update_non_density_type(parms.defined("SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE") ? parms["SINGLE_VERTEX_UPDATE_FOR_NON_DENSITY_TYPE"] : true),
pert_order_hist(max_order+1),
comm(communicator),
g0_intpl(parms),
update_manager(parms, Uijkl, g0_intpl, node==0)
{
  //other parameters
  step=0;
  start_time=time(NULL);
  measurement_time=0;
  update_time=0;

  //initialize the simulation variables
  initialize_simulation(parms);

  //submatrix update
  itime_vertex_container itime_vertices_init;
  if (params.defined("PREFIX_LOAD_CONFIG")) {
    std::ifstream is((params["PREFIX_LOAD_CONFIG"].template cast<std::string>()
                      +std::string("-node")+boost::lexical_cast<std::string>(node)+std::string(".txt")).c_str());
    load_config<typename TYPES::M_TYPE>(is, Uijkl, itime_vertices_init);
  }

  if (min_U_scale<1.0 && num_U_scale==1) {
    throw std::runtime_error("Invalid min_U_scale and num_U_scale!");
  }

  //submatrix_update =
    //new SubmatrixUpdate<M_TYPE>(
      //(parms["K_INS_MAX"] | 32), n_flavors,
      //g0_intpl, &Uijkl, beta, itime_vertices_init);
  for (int iu=0; iu<num_U_scale; ++iu) {
    U_scale_vals[iu] = num_U_scale!=1 ? iu*(1.0-min_U_scale)/(num_U_scale-1.0)+min_U_scale : 1.0;
    U_scale_index[iu] = iu;
    walkers.push_back(
      WALKER_P_TYPE(
        new SubmatrixUpdate<M_TYPE>(
        (parms["K_INS_MAX"] | 32), n_flavors,
        g0_intpl, &Uijkl, beta, itime_vertices_init)
      )
    );
  }
  submatrix_update = walkers[walkers.size()-1];

  vertex_histograms=new simple_hist *[n_flavors];
  vertex_histogram_size=100;
  for(unsigned int i=0;i<n_flavors;++i) {
    vertex_histograms[i]=new simple_hist(vertex_histogram_size);
  }

}

template<class TYPES>
InteractionExpansion<TYPES>::~InteractionExpansion()
{
  //delete submatrix_update;
}

template<class TYPES>
void InteractionExpansion<TYPES>::update()
{
  std::valarray<double> t_meas(0.0, 3*num_U_scale);

  pert_order_hist = 0.;

  for(std::size_t i=0;i<measurement_period;++i){
#ifndef NDEBUG
    std::cout << " step " << step << std::endl;
#endif
    //std::cout << " step " << step << std::endl;
    step++;

    for (int i_walker=0; i_walker<num_U_scale; ++i_walker) {
      //std::cout << "   walker " << i_walker << std::endl;
      double U_scale_walker;
      int U_index_walker;
      for (int iu=0; iu<num_U_scale; ++iu) {
        if (U_scale_index[iu]==i_walker) {
          U_index_walker = iu;
          U_scale_walker = U_scale_vals[iu];
          break;
        }
      }
      //if (U_scale_walker+0.1<random()) {
        //continue;
      //}
      //if (i_walker==0) {
        //std::cout << " U_scale " << step << " " <<  U_scale_walker << std::endl;
      //}
#ifndef NDEBUG
      std::cout << " walker " << i_walker << std::endl;
#endif
      boost::timer::cpu_timer timer;

      for (int i_ins_rem=0; i_ins_rem<n_ins_rem; ++i_ins_rem) {
        update_manager.do_ins_rem_update(*walkers[i_walker], Uijkl, random, U_scale_walker);
      }

      double t_m = timer.elapsed().wall;

      for (int i_shift=0; i_shift<n_shift; ++i_shift) {
        update_manager.do_shift_update(*walkers[i_walker], Uijkl, random, !is_thermalized());
        //update_manager.do_shift_update(*walkers[i_walker], Uijkl, random, false);
      }

      for (int i_spin_flip=0; i_spin_flip<n_spin_flip; ++i_spin_flip) {
        update_manager.do_spin_flip_update(*walkers[i_walker], Uijkl, random);
      }

      t_meas[0+3*U_index_walker] += t_m;
      t_meas[1+3*U_index_walker] += (timer.elapsed().wall-t_m);

      if(step % recalc_period ==0) {
        boost::timer::cpu_timer timer2;
        submatrix_update->recompute_matrix(true);
        t_meas[2+3*U_index_walker] += timer2.elapsed().wall;
      }

      //measurement for walker with physical U
      if (U_scale_walker==1.0) {
        if(submatrix_update->pert_order()<max_order) {
          pert_hist[submatrix_update->pert_order()]++;
        }
        assert(submatrix_update->pert_order()<pert_order_hist.size());
        ++pert_order_hist[submatrix_update->pert_order()];

        for(spin_t flavor=0; flavor<n_flavors; ++flavor) {
          vertex_histograms[flavor]->count(submatrix_update->invA()[flavor].creators().size());
        }
      }
    }//i_walker

    //swaps U_SCALE
    exchange_update();
  }//for (int i_walker=0;

  //Save pertubation order for walker with physical U
  if (params.defined("PREFIX_OUTPUT_TIME_SERIES")) {
    std::valarray<double> pert_vertex(Uijkl.n_vertex_type());
    const itime_vertex_container& itime_vertices = walkers[U_scale_index[U_scale_index.size()-1]]->itime_vertices();
    assert(U_scale_vals[U_scale_index[U_scale_index.size()-1]]==1.0);
    for (std::vector<itime_vertex>::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
      assert(it->type()>=0 && it->type()<Uijkl.n_vertex_type());
      ++pert_vertex[it->type()];
    }
    for (unsigned i_vt=0; i_vt<Uijkl.n_vertex_type(); ++i_vt) {
      pert_order_dynamics.push_back(static_cast<double>(pert_vertex[i_vt]));
    }
  }

  t_meas *= 1E-6/measurement_period;
  measurements["UpdateTimeMsec"] << t_meas;

}

template<class TYPES>
void InteractionExpansion<TYPES>::exchange_update() {
  if (num_U_scale==1) {
    return;
  }
  std::valarray<double> acc(0.0, num_U_scale-1);
  for (int iu=0; iu<num_U_scale-1; ++iu) {
    const int i_walker0 = U_scale_index[iu];
    const int i_walker1 = U_scale_index[iu+1];

    const double prob = std::pow(U_scale_vals[iu+1]/U_scale_vals[iu],
                                 walkers[i_walker0]->pert_order()-walkers[i_walker1]->pert_order());
    if (prob>random()) {
      acc[iu] = 1.0;
      std::swap(U_scale_index[iu], U_scale_index[iu+1]);
    }
  }
  submatrix_update = walkers[U_scale_index[num_U_scale-1]];
  measurements["ACCEPTANCE_RATE_EXCHANGE"] << acc;

  if (num_U_scale>2 && !is_thermalized()) {
    const double magic_number = 1.01;
    std::valarray<double> diff_U_scale(num_U_scale-1);
    double sum_tmp = 0.0;
    for (int iu=0; iu<num_U_scale-1; ++iu) {
      diff_U_scale[iu] = acc[iu]>0 ? (U_scale_vals[iu+1]-U_scale_vals[iu])*magic_number : (U_scale_vals[iu+1]-U_scale_vals[iu])/magic_number;
      sum_tmp += diff_U_scale[iu];
    }
    diff_U_scale *= (1.0-min_U_scale)/sum_tmp;
    for (int iu=0; iu<num_U_scale-2; ++iu) {
      U_scale_vals[iu+1] = U_scale_vals[iu]+diff_U_scale[iu];
    }
  }
  //if (node==0) {
    //std::cout << " Pert " << step << " ";
    //for (int iu=0; iu<num_U_scale; ++iu) {
      //std::cout << walkers[iu]->pert_order() << " ";
    //}
  //}
  //std::cout << std::endl;
  //if (node==1) {
    //std::cout << " U_scale " << step << " " <<  U_scale_index[num_U_scale-1] << " " << U_scale_index[0] << std::endl;
  //}
#ifndef NDEBUG
  std::cout << " U_scale_index: ";
  for (int iu=0; iu<num_U_scale; ++iu) {
    std::cout << U_scale_index[iu] << " ";
  }
  std::cout << std::endl;
#endif
}

template<class TYPES>
void InteractionExpansion<TYPES>::measure(){
  //In the below, real physical quantities are measured.
  std::valarray<double> timings(2);
  measure_observables(timings);
  measurements["MeasurementTimeMsec"] << timings;
}

/**
 * Finalize the Monte Carlo simulation, e.g., write some data to disk.
 */
template<class TYPES>
void InteractionExpansion<TYPES>::finalize()
{

  if (pert_order_dynamics.size()>0) {
    std::ofstream ofs((params["PREFIX_OUTPUT_TIME_SERIES"].template cast<std::string>()+std::string("-pert_order-node")+boost::lexical_cast<std::string>(node)+std::string(".txt")).c_str());
    unsigned int n_data = pert_order_dynamics.size()/Uijkl.n_vertex_type();
    unsigned int i=0;
    for (unsigned int i_data=0; i_data<n_data; ++i_data) {
      ofs << i_data << " ";
      for (unsigned int i_vt = 0; i_vt < Uijkl.n_vertex_type(); ++i_vt) {
        ofs << pert_order_dynamics[i] << " ";
        ++i;
      }
      ofs << std::endl;
    }
  }

  if (Wk_dynamics.size()>0) {
    std::ofstream ofs((params["PREFIX_OUTPUT_TIME_SERIES"].template cast<std::string>()+std::string("-Wk-node")+boost::lexical_cast<std::string>(node)+std::string(".txt")).c_str());
    unsigned int n_data = Wk_dynamics.size()/(n_flavors*n_site);
    unsigned int i=0;
    for (unsigned int i_data=0; i_data<n_data; ++i_data) {
      ofs << i_data << " ";
      for (unsigned int flavor=0; flavor<n_flavors; ++flavor) {
        for (unsigned int site1 = 0; site1 < n_site; ++site1) {
          ofs << Wk_dynamics[i].real() << " " << Wk_dynamics[i].imag() << "   ";
          ++i;
        }
      }
      ofs << std::endl;
    }
  }

  if (Sl_dynamics.size()>0) {
    std::ofstream ofs((params["PREFIX_OUTPUT_TIME_SERIES"].template cast<std::string>()+std::string("-Sl-node")+boost::lexical_cast<std::string>(node)+std::string(".txt")).c_str());
    unsigned int n_data = Sl_dynamics.size()/(n_flavors*n_site);
    unsigned int i=0;
    for (unsigned int i_data=0; i_data<n_data; ++i_data) {
      ofs << i_data << " ";
      for (unsigned int flavor=0; flavor<n_flavors; ++flavor) {
        for (unsigned int site1 = 0; site1 < n_site; ++site1) {
          ofs << Sl_dynamics[i].real() << " " << Sl_dynamics[i].imag() << "   ";
          ++i;
        }
      }
      ofs << std::endl;
    }
  }

  if (params.defined("PREFIX_DUMP_CONFIG")) {
    std::ofstream os((params["PREFIX_DUMP_CONFIG"].template cast<std::string>()
                      +std::string("-node")+boost::lexical_cast<std::string>(node)+std::string(".txt")).c_str());
    dump(os, submatrix_update->itime_vertices());
  }

}


template<class TYPES>
double InteractionExpansion<TYPES>::fraction_completed() const{
  if (!is_thermalized()) {
    return 0.;
  } else {
    return ((step - therm_steps) / (double) mc_steps);
  }
}



///do all the setup that has to be done before running the simulation.
template<class TYPES>
void InteractionExpansion<TYPES>::initialize_simulation(const alps::params &parms)
{
  pert_hist.clear();
  initialize_observables();
}

template<class TYPES>
bool InteractionExpansion<TYPES>::is_quantum_number_conserved(const itime_vertex_container& vertices) {
  const int Nv = vertices.size();

  if (Nv==0)
    return true;

  std::valarray<int> qn_t(0, qn_dim), qn_max(0, qn_dim), qn_min(0, qn_dim);
  itime_vertex_container vertices_sorted(vertices);//sort vertices in decreasing order (in time)
  std::sort(vertices_sorted.begin(), vertices_sorted.end());

  for (int iv=0; iv<Nv; ++iv) {
    const vertex_definition<M_TYPE> vd = Uijkl.get_vertex(vertices_sorted[iv].type());
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
/*
template<class TYPES>
bool InteractionExpansion<TYPES>::is_quantum_number_within_range(const itime_vertex_container& vertices) {
  const int Nv = vertices.size();

  if (Nv==0)
    return true;

  std::valarray<int> qn_t(0, qn_dim), qn_max(0, qn_dim), qn_min(0, qn_dim);
  itime_vertex_container vertices_sorted(vertices);//sort vertices in decreasing order (in time)
  std::sort(vertices_sorted.begin(), vertices_sorted.end());

  for (int iv=0; iv<Nv; ++iv) {
    const vertex_definition<M_TYPE> vd = Uijkl.get_vertex(vertices_sorted[iv].type());
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
*/

template<class TYPES>
void InteractionExpansion<TYPES>::sanity_check() {
#ifndef NDEBUG
#endif
}



template<class TYPES>
void InteractionExpansion<TYPES>::prepare_for_measurement()
{
  if (num_U_scale>2) {
    std::valarray<double> U_scale_vals_reduced(0.0, num_U_scale);
    comm.barrier();
    boost::mpi::all_reduce(comm, &U_scale_vals[0], num_U_scale, &U_scale_vals_reduced[0], std::plus<double>());
    U_scale_vals_reduced /= U_scale_vals_reduced[num_U_scale-1];
    for (int iu=1; iu<num_U_scale-1; ++iu) {
      U_scale_vals[iu] = U_scale_vals_reduced[iu];
    }
  }

  /*
  this->statistics_ins.reset();
  this->statistics_rem.reset();
  this->statistics_shift.reset();
  this->simple_statistics_ins.reset();
  this->simple_statistics_rem.reset();
  this->statistics_dv_ins.reset();
  this->statistics_dv_rem.reset();
  */
}
