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

#ifndef DMFT_QMC_WEAK_COUPLING_H
#define DMFT_QMC_WEAK_COUPLING_H

#include <algorithm>
#include <fstream>
//#include <functional>
#include <cmath>

#include <boost/multi_array.hpp>
#include <boost/timer/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/multi_array.hpp>

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/alea.h>
#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

#include "submatrix.hpp"
#include "green_function.h"
#include "types.h"
#include "fouriertransform.h"
#include "U_matrix.h"
#include "operator.hpp"
#include "green_matrix.hpp"
#include "legendre.h"
#include "update_statistics.h"
#include "wang_landau.h"

/*types*/
class c_or_cdagger;
typedef class histogram simple_hist;

class histogram
{
public:
  
  histogram(unsigned int N):hist_(N, 0){}
  
  unsigned long &operator[](unsigned int n){return hist_[n];}
  const unsigned long &operator[](unsigned int n) const{return hist_[n];}
  unsigned int size() const{return hist_.size();}

  void count(size_t k) {
    if (k<size()) {
      ++hist_[k];
    }
  }
  
  unsigned int max_index() const
  { 
    unsigned int max_index=0; 
    double max=0; 
    for(unsigned int i=0;i<hist_.size();++i){
      if(max<hist_[i]){
        max=hist_[i];
        max_index=i;
      }
    }
    return max_index;
  }
  
  unsigned int top_index() const
  { 
    unsigned int top_index=0;
    for(unsigned int i=0;i<hist_.size();++i){
      if(hist_[i]!=0){
        top_index=i;
      }
    }
    return top_index;
  }
  
  double max(const unsigned int index) const
  { 
    double max=0; 
    for(unsigned int i=0;i<index;++i){
      if(max<hist_[i]){
        max=hist_[i];
      }
    }
    return max;
  }
  
  double average(const unsigned int index) const{ 
    double average=0; 
    for(unsigned int i=0;i<index;++i){
      average+=hist_[i];
    }
    return average/index;
  }
  
  bool is_flat(const unsigned int index)const{return max(index)*0.8<average(index);}
  
  void clear()
  {
    for(unsigned int i=0;i<hist_.size();++i){
      hist_[i]=0;
    }
  }

  std::valarray<double> to_valarray() {
    std::valarray<double> tmparray(hist_.size());
    for (size_t i=0; i<hist_.size(); ++i) {
      tmparray[i] = hist_[i];
    }
    return tmparray;
  }

private:
  
  std::vector<unsigned long> hist_;
};

class SymmExpDist {
public:
    SymmExpDist() : a_(0), b_(0), beta_(0), coeff_(0), coeffX_(0) {}
    SymmExpDist(double a_in, double b_in, double beta_in) : a_(a_in), b_(b_in), beta_(beta_in),
                                                            coeff_(1/((2/a_)*(1-std::exp(-a_*beta_))+b_*beta_)),
                                                            coeffX_(2*(1-std::exp(-a_*beta_))/(a_*beta_)+b_) {}

    double operator()(double dtau) const {return coeff_*bare_value(dtau);}
    double bare_value(double dtau) const {return std::exp(-a_*dtau)+std::exp(-a_*(beta_-dtau))+ b_;}
    double coeff_X(double dtau) const {return coeffX_;}

private:
    double a_, b_, beta_, coeff_, coeffX_;
};

typedef struct real_number_solver {
  typedef double M_TYPE;
  typedef double REAL_TYPE;
  typedef std::complex<double> COMPLEX_TYPE;
} real_number_solver;

typedef struct complex_number_solver {
    typedef std::complex<double> M_TYPE;
    typedef double REAL_TYPE;
    typedef std::complex<double> COMPLEX_TYPE;
} complex_number_solver;


template<typename T>
class BareGreenInterpolate {
public:
    BareGreenInterpolate(const alps::params& p);

    T operator()(const annihilator& c, const creator& cdagg) const;

private:
    const double beta_, temp_;
    const int ntau_, n_flavors_, n_sites_;
    boost::multi_array<std::pair<T,T>,4> AB_;
    double dbeta_;//beta/ntau
};

class InteractionExpansionBase: public alps::mcbase {
public:
    InteractionExpansionBase(const alps::params &p, int rank, const boost::mpi::communicator &communicator) : alps::mcbase(p,rank) {};
    virtual ~InteractionExpansionBase() {}
    virtual bool is_thermalized() const=0;
    virtual void update()=0;
    virtual void measure()=0;
    virtual void finalize()=0;
    virtual double fraction_completed() const=0;
};

template<class TYPES>
class InteractionExpansion: public InteractionExpansionBase //alps::mcbase
{
public:

  InteractionExpansion(const alps::params& p, int rank, const boost::mpi::communicator& communicator);
  ~InteractionExpansion();
  bool is_thermalized() const {return step>therm_steps;}
  void update();
  void measure();
  void finalize();
  double fraction_completed() const;

  typedef typename TYPES::M_TYPE M_TYPE;
  typedef typename TYPES::REAL_TYPE REAL_TYPE;
  typedef typename TYPES::COMPLEX_TYPE COMPLEX_TYPE;
  typedef boost::multi_array<std::complex<double>, 4> Wk_t;
  typedef green_function<typename TYPES::COMPLEX_TYPE> matsubara_green_function_t;
  typedef green_function<typename TYPES::COMPLEX_TYPE> itime_green_function_t;

protected:

    /*
  template<typename SOLVER>
  class SPLINE_G0_HELPER {
  public:
      SPLINE_G0_HELPER(SOLVER* solver) : solver_(solver) {}

      typename SOLVER::M_TYPE operator() (const annihilator& c, const creator& cdagger) {
        return solver_->green0_spline_for_M(c, cdagger);
      };
      SOLVER* solver_;
  };
     */

  /*functions*/
  /*io & initialization*/
  void initialize_simulation(const alps::params &parms); // called by constructor

  // in file io.cpp
  void print(std::ostream &os);
  
  /*green's function*/
  // in file spines.cpp
  M_TYPE green0_spline_for_M(const annihilator& c, const creator& cdagger);//with correct treatment of equal-time Green's function
  M_TYPE green0_spline_new(const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2);

  // in file observables.ipp
  void measure_observables(std::valarray<double>& timings);
  void initialize_observables(void);

  // in file selfenergy.ipp
  void compute_W_matsubara();
  void compute_Sl();
  void measure_Wk(Wk_t& Wk, const unsigned int nfreq);
  void measure_densities();

  // in file interaction_expansion.hpp
  void sanity_check();
  bool is_irreducible(const itime_vertex_container& vertices);
  bool is_quantum_number_conserved(const itime_vertex_container& vertices);
  bool is_quantum_number_within_range(const itime_vertex_container& vertices);

  void prepare_for_measurement(); //called once after thermalization is done

  /*private member variables, constant throughout the simulation*/
  const unsigned int node;
  const unsigned int max_order;                        
  const spin_t n_flavors;                                //number of flavors (called 'flavors') in InteractionExpansion
  const site_t n_site;                                //number of sites
  const frequency_t n_matsubara;        //number of matsubara freq
  const frequency_t n_matsubara_measurements;        //number of measured matsubara freq
  const itime_index_t n_tau;                        //number of imag time slices
  const itime_t n_tau_inv;                        //the inverse of n_tau
  const frequency_t n_self;                        //number of self energy (W) binning points
  const std::size_t n_legendre;
  const boost::uint64_t mc_steps;                        
  const unsigned long therm_steps;                
  const double max_time_in_seconds;
  const size_t n_multi_vertex_update;
  const int n_ins_rem;
  const int n_shift;
  const int n_spin_flip;
  const bool force_quantum_number_conservation;
  const bool single_vertex_update_non_density_type;
  const double beta;
  const double temperature;                        //only for performance reasons: avoid 1/beta computations where possible

  general_U_matrix<M_TYPE> Uijkl; //for any general two-body interaction

  /*heart of submatrix update*/
  SubmatrixUpdate<M_TYPE> *submatrix_update;

  //for measurement of Green's function
  //M is computed from A in measure_observables.
  std::vector<alps::numeric::matrix<M_TYPE> > M_flavors;

  //quantum numbers
  std::vector<std::vector<std::vector<size_t> > > groups;
  std::vector<std::vector<int> > group_map;
  std::vector<std::vector<quantum_number_t> > quantum_number_vertices;
  std::vector<int> group_dim;
  int qn_dim;

  //double vertex update
  std::vector<std::pair<int,int> > mv_update_valid_pair;
  boost::multi_array<bool,2> mv_update_valid_pair_flag;

  std::vector<bool> is_density_density_type;

  //for shift update
  std::vector<bool> shift_update_valid;

  const unsigned int recalc_period;                
  const unsigned int measurement_period;
  const unsigned int convergence_check_period;
  
  /*InteractionExpansion's roundoff threshold*/
  const double almost_zero;                        
  /*PRNG seed*/
  const int seed;
  bool is_thermalized_in_previous_step_;

  matsubara_green_function_t bare_green_matsubara;
  itime_green_function_t bare_green_itime;
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  
  //for window update
  double window_width;
  boost::random::exponential_distribution<> window_dist;
  SymmExpDist symm_exp_dist;

  simple_hist pert_hist;
  unsigned int hist_max_index;
  simple_hist **vertex_histograms;
  unsigned int vertex_histogram_size;
  
  unsigned long step;        
  time_t start_time;
  clock_t update_time;
  clock_t measurement_time;

  //Statistics about multi-vertex updates (imaginary time information)
  scalar_histogram_flavors statistics_rem, statistics_ins, statistics_shift, statistics_dv_rem, statistics_dv_ins;

  //only acceptance rate
  simple_update_statistcs simple_statistics_rem, simple_statistics_ins;
  int num_accepted_shift;

  //just for test
  update_proposer update_prop;

  LegendreTransformer legendre_transformer;

  std::valarray<double> pert_order_hist;

  const boost::mpi::communicator& comm;

  //only for test
  std::vector<typename TYPES::COMPLEX_TYPE> Wk_dynamics;
  std::vector<double> pert_order_dynamics;

  BareGreenInterpolate<M_TYPE> g0_intpl;
};

/*aux functions*/
std::ostream& operator << (std::ostream& os, const std::vector<double>& v);
std::ostream& operator << (std::ostream &os, const c_or_cdagger &c);
std::ostream& operator << (std::ostream& os, const simple_hist &h);

#include "interaction_expansion.ipp"
#include "selfenergy.ipp"
#include "io.ipp"
#include "splines.ipp"
#include "observables.ipp"

#endif
