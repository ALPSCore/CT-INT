//
// Created by H. Shinaoka on 2016/02/19.
//

#ifndef IMPSOLVER_UPDATE_MANAGER_HPP
#define IMPSOLVER_UPDATE_MANAGER_HPP

#include <vector>

#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <alps/mcbase.hpp>
#include <alps/ngs.hpp>

#include "submatrix.hpp"
#include "green_function.h"
#include "U_matrix.h"
//#include "update_statistics.h"

template<typename T, typename SPLINE_G0>
std::pair<T,T>
compute_weight(const general_U_matrix<T>& Uijkl, const SPLINE_G0& spline_G0, const itime_vertex_container& itime_vertices);

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

template<typename T>
class VertexUpdateManager {
public:
  typedef double insertion_removal_info_type;

  template<typename G0_SPLINE>
  VertexUpdateManager(const alps::params& parms, const general_U_matrix<T>& Uijkl, const G0_SPLINE& g0_spline, bool message);

  //Called once before Monte Carlo sampling
  //template<typename M>
  //void create_observables(M&);

  //Called at each Monte Carlo step once thermalization is reached
  //template<typename M>
  //void measure_observables(M&) const;

  //fix parameters
  //void finish_learning();
  template<typename R>
  T do_ins_rem_update(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random, double U_scale);

  template<typename R>
  T do_spin_flip_update(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random);

  template<typename R>
  T do_shift_update(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random, bool tune_step_size);

private:
  template<typename R>
  T insertion_step(SubmatrixUpdate<T>& submatrix, R& random, int vertex_begin, int num_vertices_ins, double U_scale);

  template<typename R>
  T removal_step(SubmatrixUpdate<T>& submatrix, R& random, double U_scale);

  template<typename R>
  std::pair<std::vector<itime_vertex>, double>
  gen_itime_vertices_insertion(const general_U_matrix<T>& Uijkl, R& random01) const;

  template<typename R>
  double acc_rate_corr_insertion(const itime_vertex_container& itime_vertices_current, const std::vector<int>& pos_vertices, R& random01) const;
  void feedback_insertion(const insertion_removal_info_type& info, bool accepted);

  template<typename R>
  double pick_up_vertices_to_be_removed(const itime_vertex_container& itime_vertices_current, R& random01, std::vector<int>& pos_vertices) const;

  template<typename R>
   T spin_flip_step(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random, int pos_vertex);

  const double beta;
  const int n_flavors;
  const int k_ins_max;
  const int max_order;
  const int num_vertex_type;
  const int n_multi_vertex_update;
  const int n_shift;

  //for single vertex update
  std::vector<int> sv_update_vertices;
  std::vector<bool> sv_update_vertices_flag;

  //for double vertex update
  std::vector<std::pair<int,int> > mv_update_valid_pair;
  boost::multi_array<bool,2> mv_update_valid_pair_flag;

  //for quantum number
  std::vector<std::vector<std::vector<size_t> > > groups;
  std::vector<std::vector<int> > group_map;
  std::vector<std::vector<quantum_number_t> > quantum_number_vertices;
  std::vector<int> group_dim;
  int qn_dim;

  //for window update
  //double window_width;
  //boost::random::exponential_distribution<> window_dist;
  SymmExpDist symm_exp_dist;

  //for shift update
  std::vector<bool> shift_update_valid;
  int num_shift_valid_vertex_types;
  double shift_step_size;

  //only acceptance rate
  simple_update_statistcs simple_statistics_rem, simple_statistics_ins;
  int num_accepted_shift;

  //Statistics about multi-vertex updates (imaginary time information)
  scalar_histogram_flavors statistics_rem, statistics_ins, statistics_shift, statistics_dv_rem, statistics_dv_ins;

  boost::random::discrete_distribution<> Nv_m1_dist;
};

template<typename T>
template<typename G0_SPLINE>
VertexUpdateManager<T>::VertexUpdateManager(const alps::params &parms, const general_U_matrix<T>& Uijkl, const G0_SPLINE &g0_spline, bool message)
  : beta(parms["BETA"]),
    n_flavors(parms["FLAVORS"]),
    k_ins_max(parms["K_INS_MAX"] | 32),
    max_order(parms["MAX_ORDER"] | 2048),
    num_vertex_type(Uijkl.get_vertices().size()),
    sv_update_vertices(),
    sv_update_vertices_flag(num_vertex_type, false),
    n_multi_vertex_update(parms["N_MULTI_VERTEX_UPDATE"] | 1),
    n_shift(parms["N_VERTEX_SHIFT"] | 1),
    shift_update_valid(num_vertex_type, false),
    num_shift_valid_vertex_types(0),
    shift_step_size(parms["VERTEX_SHIFT_STEP_SIZE"] | 0.1*beta),
    simple_statistics_ins(n_multi_vertex_update),
    simple_statistics_rem(n_multi_vertex_update),
    statistics_rem((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, n_multi_vertex_update-1),
    statistics_ins((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, n_multi_vertex_update-1),
    statistics_shift((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, 1),
    statistics_dv_rem(0, 0, 0),
    statistics_dv_ins(0, 0, 0)
{
  const double almost_zero = 1E-10;

  if (n_multi_vertex_update==1) {
    const std::vector<vertex_definition<T> >& v_defs_tmp = Uijkl.get_vertices();
    sv_update_vertices.resize(v_defs_tmp.size());
    for (int iv=0; iv<sv_update_vertices.size(); ++iv) {
      sv_update_vertices[iv] = iv;
    }
    std::fill(sv_update_vertices_flag.begin(), sv_update_vertices_flag.end(), true);
  } else {
    const std::vector<vertex_definition<T> >& v_defs_tmp = Uijkl.get_vertices();
    assert(v_defs_tmp.size()==num_vertex_type);
    sv_update_vertices.resize(0);
    for (int iv=0; iv<num_vertex_type; ++iv) {
      if (v_defs_tmp[iv].is_density_type()) {
        sv_update_vertices.push_back(iv);
        sv_update_vertices_flag[iv] = true;
      }
    }
  }

  if (n_multi_vertex_update>1) {
    quantum_number_vertices = make_quantum_numbers<T,T,G0_SPLINE>(g0_spline, Uijkl.get_vertices(), groups, group_map, almost_zero);
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

    symm_exp_dist = SymmExpDist(parms["DOUBLE_VERTEX_UPDATE_A"], parms["DOUBLE_VERTEX_UPDATE_B"], beta);

    statistics_dv_ins = scalar_histogram_flavors((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, mv_update_valid_pair.size());
    statistics_dv_rem = scalar_histogram_flavors((parms["N_TAU_UPDATE_STATISTICS"] | 100), beta, mv_update_valid_pair.size());
  }

  if (n_shift>0) {
    /*
    const std::vector<vertex_definition<T> >& v_defs_tmp = Uijkl.get_vertices();
    for (int iv=0; iv<sv_update_vertices_flag.size(); ++iv) {
      if (!sv_update_vertices_flag[iv]) {
        shift_update_valid[iv] = true;
      }
    }
    */
    std::fill(shift_update_valid.begin(), shift_update_valid.end(), true);
    num_shift_valid_vertex_types = std::accumulate(shift_update_valid.begin(),shift_update_valid.end(),0);
  }

  //set up parameters for updates
  std::vector<double> proposal_prob(n_multi_vertex_update, 1.0);
  if (parms.defined("MULTI_VERTEX_UPDATE_PROPOSAL_RATE")) {
    proposal_prob.resize(0);
    std::stringstream ss(parms["MULTI_VERTEX_UPDATE_PROPOSAL_RATE"].template cast<std::string>());
    double rtmp;
    while (ss >> rtmp) {
      proposal_prob.push_back(rtmp);
    }
    if (proposal_prob.size()!=n_multi_vertex_update)
      throw std::runtime_error("The number of elements in MULTI_VERTEX_UPDATE_PROPOSAL_RATE is different from N_MULTI_VERTEX_UPDATE");
  }
  Nv_m1_dist = boost::random::discrete_distribution<>(proposal_prob.begin(), proposal_prob.end());

  if (message) {
    std::cout << std::endl;
    std::cout << std::endl << "Single-vertex updates will be performed for the following vertices." << std::endl;
    for (int iv=0; iv<sv_update_vertices.size(); ++iv) {
      std::cout << " iv = " << sv_update_vertices[iv] << std::endl;
    }
    std::cout << std::endl;
  }

  if(message && n_multi_vertex_update>1) {
    std::cout << std::endl;
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

  //shift update
  if (num_shift_valid_vertex_types>0 && message) {
    std::cout << std::endl << "Shift updates will be performed for the following vertices." << std::endl;
    for(int iv=0; iv<shift_update_valid.size(); ++iv) {
      if (shift_update_valid[iv]) {
        std::cout << " iv = " << iv << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

template<typename T>
template<typename R>
T VertexUpdateManager<T>::do_ins_rem_update(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random, double U_scale) {

  int num_ins_try = 0;
  std::vector<bool> try_ins(2*k_ins_max, false);
  for (int i_update=0; i_update<2*k_ins_max; ++i_update) {
    if (random()<0.5) {
      try_ins[i_update] = true;
      ++num_ins_try;
    }
  }

  std::vector<int> pos_vertices_ins(num_ins_try);
  std::vector<int> num_vertices_ins(num_ins_try);

  T weight_rat = 1.0;

  std::vector<insertion_removal_info_type> ins_info(num_ins_try);

  //add non-interacting vertices
  const int Nv0 = submatrix.pert_order();
  int vertex_begin = Nv0;
  std::vector<itime_vertex> new_vertices_all;
  for (int i_ins=0; i_ins<num_ins_try; ++i_ins) {
    std::vector<itime_vertex> new_vertices;
    boost::tie(new_vertices,ins_info[i_ins]) = gen_itime_vertices_insertion(Uijkl, random);
    const int Nv = new_vertices.size();

    for (int iv=0; iv<new_vertices.size(); ++iv) {
      new_vertices[iv].set_non_interacting();//this does not modify the spin state but just hide it.
      new_vertices_all.push_back(new_vertices[iv]);
    }
    pos_vertices_ins[i_ins] = vertex_begin;
    num_vertices_ins[i_ins] = new_vertices.size();
    vertex_begin += new_vertices.size();
  }
  assert(pos_vertices_ins.size()==num_vertices_ins.size());
  assert(new_vertices_all.size()==std::accumulate(num_vertices_ins.begin(),num_vertices_ins.end(),0));

  //set starting point
  submatrix.init_update(new_vertices_all);

  //perform actual updates
  int i_ins = 0;
  for (int i_update=0; i_update<2*k_ins_max; ++i_update) {
    T rtmp;
    if (try_ins[i_update]) {
      assert(i_ins<pos_vertices_ins.size());
      if (submatrix.itime_vertices().num_interacting()+num_vertices_ins[i_ins]>max_order) {
        rtmp = 1.0;
      } else {
        rtmp = insertion_step(submatrix, random, pos_vertices_ins[i_ins], num_vertices_ins[i_ins], U_scale);
      }
      ++i_ins;
    } else {
      rtmp = removal_step(submatrix, random, U_scale);
    }
    weight_rat *= rtmp;
  }
  assert(i_ins==num_ins_try);

  //update A^{-1}
  submatrix.finalize_update();

  return weight_rat;
};

template<typename T>
template<typename R>
T VertexUpdateManager<T>::insertion_step(SubmatrixUpdate<T>& submatrix, R& random, int vertex_begin, int num_vertices_ins, double U_scale) {
  assert(vertex_begin+num_vertices_ins<=submatrix.itime_vertices().size());

  if (num_vertices_ins==0) {
    return 1.0;
  }

  T det_rat_A, f_rat, U_rat;

  const itime_vertex_container& itime_vertices = submatrix.itime_vertices();

  std::vector<int> new_spins_work(num_vertices_ins);
  std::vector<int> pos_vertices_work(num_vertices_ins);
  for (int iv=0; iv<num_vertices_ins; ++iv) {
    assert(iv+vertex_begin<submatrix.itime_vertices().size());
    assert(itime_vertices[iv+vertex_begin].is_non_interacting());
    new_spins_work[iv] = itime_vertices[iv+vertex_begin].af_state();
    pos_vertices_work[iv] = iv+vertex_begin;
  }

  boost::tie(det_rat_A,f_rat,U_rat) = submatrix.try_spin_flip(pos_vertices_work, new_spins_work);
  const double acc_corr = acc_rate_corr_insertion(itime_vertices, pos_vertices_work, random);

  T prob = det_rat_A*f_rat*U_rat*acc_corr*std::pow(U_scale, 1.0*num_vertices_ins);

  if (std::abs(prob)>random()) {
    //std::cout << "accepted " << std::endl;
    submatrix.perform_spin_flip(pos_vertices_work, new_spins_work);
    return det_rat_A*f_rat*U_rat;
  } else {
    //std::cout << "rejected " << std::endl;
    submatrix.reject_spin_flip();
    return 1.0;
  }
}

template<typename T>
template<typename R>
T VertexUpdateManager<T>::removal_step(SubmatrixUpdate<T>& submatrix, R& random, double U_scale) {
  //const int Nv = submatrix.pert_order();
  std::vector<int> pos_vertices_remove;
  const double acc_corr = pick_up_vertices_to_be_removed(submatrix.itime_vertices(), random, pos_vertices_remove);
  const int nv_rem = pos_vertices_remove.size();

  //no vertices to be removed
  if (nv_rem==0) {
    return 1.0;
  }

#ifndef NDEBUG
  for (int iv=0; iv<pos_vertices_remove.size(); ++iv) {
    assert(!submatrix.itime_vertices()[pos_vertices_remove[iv]].is_non_interacting());
  }
#endif

  std::vector<int> new_spins_remove(nv_rem, NON_INT_SPIN_STATE);

  T det_rat_A, f_rat, U_rat;
  boost::tie(det_rat_A,f_rat,U_rat) = submatrix.try_spin_flip(pos_vertices_remove, new_spins_remove);

  T prob = det_rat_A*f_rat*U_rat*acc_corr*std::pow(U_scale, -1.0*nv_rem);

  if (std::abs(prob)>random()) {
    //std::cout << "accepted " << std::endl;
    submatrix.perform_spin_flip(pos_vertices_remove, new_spins_remove);
    return det_rat_A*f_rat*U_rat;
  } else {
    //std::cout << "rejected " << std::endl;
    submatrix.reject_spin_flip();
    return 1.0;
  }
}

//template<typename T>
//template<typename M>
//void VertexUpdateManager<T>::create_observables(M& measurements) {
  //measurements << alps::ngs::SimpleRealObservable("AcceptanceRateShift");
  /*
  measurements << alps::ngs::SimpleRealObservable("AcceptanceRateShift");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexInsertion");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexRemoval");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexShift");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexInsertion_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexRemoval_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexShift_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexInsertion_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexRemoval_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsVertexShift_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexInsertion");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexInsertion_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexInsertion_sum");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexRemoval");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexRemoval_count");
  measurements << alps::ngs::SimpleRealVectorObservable("StatisticsDoubleVertexRemoval_sum");
  */
//}

//template<typename T>
//template<typename M>
//void VertexUpdateManager<T>::measure_observables(M& measurements) {
  /*
  if (n_multi_vertex_update>1) {
    measurements["StatisticsVertexInsertion"] << statistics_ins.get_mean();
    measurements["StatisticsVertexRemoval"] << statistics_rem.get_mean();
    measurements["StatisticsDoubleVertexInsertion"] << statistics_dv_ins.get_mean();
    measurements["StatisticsDoubleVertexRemoval"] << statistics_dv_rem.get_mean();

    measurements["StatisticsVertexInsertion_count"] << statistics_ins.get_counter();
    measurements["StatisticsVertexRemoval_count"] << statistics_rem.get_counter();
    measurements["StatisticsDoubleVertexInsertion_count"] << statistics_dv_ins.get_counter();
    measurements["StatisticsDoubleVertexRemoval_count"] << statistics_dv_rem.get_counter();

    measurements["StatisticsVertexInsertion_sum"] << statistics_ins.get_sumval();
    measurements["StatisticsVertexRemoval_sum"] << statistics_rem.get_sumval();
    measurements["StatisticsDoubleVertexInsertion_sum"] << statistics_dv_ins.get_sumval();
    measurements["StatisticsDoubleVertexRemoval_sum"] << statistics_dv_rem.get_sumval();
  }

  if (n_shift>0) {
    measurements["AcceptanceRateShift"] << num_accepted_shift/((double) measurement_period * (double) n_shift);
  }

  measurements["StatisticsVertexShift"] << statistics_shift.get_mean();
  measurements["StatisticsVertexShift_count"] << statistics_shift.get_counter();
  measurements["StatisticsVertexShift_sum"] << statistics_shift.get_sumval();

  statistics_ins.reset();
  statistics_rem.reset();
  statistics_shift.reset();
  statistics_dv_ins.reset();
  statistics_dv_rem.reset();

  for (int iv=0; iv<n_multi_vertex_update; ++iv){
    measurements["VertexInsertion_"+boost::lexical_cast<std::string>(iv+1)] << simple_statistics_ins.get_result(iv);
    measurements["VertexRemoval_"+boost::lexical_cast<std::string>(iv+1)] << simple_statistics_rem.get_result(iv);
  }
  simple_statistics_ins.reset();
  simple_statistics_rem.reset();
  */
//}

template<typename T>
template<typename R>
std::pair<std::vector<itime_vertex>, double>
VertexUpdateManager<T>::gen_itime_vertices_insertion(const general_U_matrix<T>& Uijkl, R &random01) const {
  const int Nv = Nv_m1_dist(random01)+1;

  std::vector<itime_vertex> vertices;
  double dist;

  if (Nv==1) {
    if (sv_update_vertices.size()==0) {
      return std::make_pair(std::vector<itime_vertex>(), 0.0);
    }
    vertices.resize(1);
    const double time = random01()*beta;
    const int iv_rnd = static_cast<int>(random01()*sv_update_vertices.size());
    const vertex_definition<T>& vdef = Uijkl.get_vertex(sv_update_vertices[iv_rnd]);
    const int af_state = static_cast<size_t>(random01()*vdef.num_af_states());
    vertices[0] = itime_vertex(vdef.id(), af_state, time, vdef.rank(), vdef.is_density_type());
  } else if (Nv==2) {
    if (mv_update_valid_pair.size()==0) {
      return std::make_pair(std::vector<itime_vertex>(), 0.0);
    }
    std::pair<int,int> v_pair;
    v_pair = mv_update_valid_pair[mv_update_valid_pair.size()*random01()];
    vertices = generate_valid_vertex_pair2(Uijkl,v_pair,random01,beta,symm_exp_dist);
  } else {
    throw std::runtime_error("Nv>2 not implemented");
  }
  return std::make_pair(vertices, dist);
}

template<typename T>
template<typename R>
double
VertexUpdateManager<T>::acc_rate_corr_insertion(const itime_vertex_container& itime_vertices_current, const std::vector<int>& pos_vertices, R& random01) const {
  const int Nv_updated = pos_vertices.size();
  if (Nv_updated==1) {
    int num_cand = 0;
    for (int iv=0; iv<itime_vertices_current.size(); ++iv) {
      if (!itime_vertices_current[iv].is_non_interacting()
          && sv_update_vertices_flag[itime_vertices_current[iv].type()]) {
        ++num_cand;
      }
    }
    return (beta*sv_update_vertices.size())/(num_cand+1.0);
  } else if (Nv_updated==2) {
    itime_vertex_container itime_vertices_new(itime_vertices_current);
    for (int iv=0; iv<Nv_updated; ++iv) {
      itime_vertices_new[pos_vertices[iv]].set_interacting();
    }

    const itime_vertex& v0 = itime_vertices_current[pos_vertices[0]];
    const itime_vertex& v1 = itime_vertices_current[pos_vertices[1]];

    const double dtau = mymod(v0.time()-v1.time(), beta);
    int n_vpair;
    double F;

    pick_up_valid_vertex_pair2(itime_vertices_new,
                               std::make_pair(v0.type(),v1.type()),
                               beta, symm_exp_dist, random01, n_vpair, F);
    if (n_vpair==0) {
      throw std::logic_error("v_pair must be larger than 0.");
    }
    return (beta*beta)*symm_exp_dist.coeff_X(dtau)/F;
  } else {
    throw std::runtime_error("Nv>2 not implemented");
  }
}

template<typename T>
template<typename R>
double VertexUpdateManager<T>::pick_up_vertices_to_be_removed(const itime_vertex_container& itime_vertices_current, R& random01,
                                                                   std::vector<int>& pos_vertices) const {
  const int Nv_updated = Nv_m1_dist(random01)+1;
  //choose vertices to be removed
  if (Nv_updated==1) {
    int num_v = 0;
    pos_vertices.resize(0);
    pos_vertices.reserve(itime_vertices_current.size());
    for (int iv=0; iv<itime_vertices_current.size(); ++iv) {
      if (!itime_vertices_current[iv].is_non_interacting()
          && sv_update_vertices_flag[itime_vertices_current[iv].type()]) {
        pos_vertices.push_back(iv);
      }
    }
    if (pos_vertices.size()>0) {
      const int num_cand = pos_vertices.size();
      std::swap(pos_vertices[0],
                pos_vertices[
                  static_cast<int>(pos_vertices.size()*random01())
                ]
      );
      pos_vertices.resize(1);
      return num_cand/(beta*sv_update_vertices.size());
    } else {
      return 0.0;
    }
  } else if (Nv_updated==2) {
    if (mv_update_valid_pair.size()==0) {
      pos_vertices.resize(0);
      return 0.0;
    }
    int n_vpair;
    double F;
    std::pair<int,int> v_pair = mv_update_valid_pair[mv_update_valid_pair.size()*random01()];
    std::pair<int,int> r = pick_up_valid_vertex_pair2(itime_vertices_current, v_pair, beta, symm_exp_dist, random01, n_vpair, F);
    if (n_vpair>0) {
      pos_vertices.resize(2);
      pos_vertices[0]=r.first;
      pos_vertices[1]=r.second;
      const double t0 = itime_vertices_current[pos_vertices[0]].time();
      const double t1 = itime_vertices_current[pos_vertices[1]].time();
      assert(itime_vertices_current[pos_vertices[0]].type()==v_pair.first);
      assert(itime_vertices_current[pos_vertices[1]].type()==v_pair.second);
      assert(!itime_vertices_current[pos_vertices[0]].is_non_interacting());
      assert(!itime_vertices_current[pos_vertices[1]].is_non_interacting());
      return F/((beta*beta)*symm_exp_dist.coeff_X(mymod(t1-t0, beta)));
    } else {
      pos_vertices.resize(0);
      return 0.0;
    }
  } else {
    throw std::runtime_error("Nv>2 not implemented");
  }
}

/***
 * Spin flip update
 */
template<typename T>
template<typename R>
T VertexUpdateManager<T>::do_spin_flip_update(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random) {

  const int Nv0 = submatrix.pert_order();
  const int nv_flip = std::min(2*k_ins_max, Nv0);

  T weight_rat = 1.0;

  std::vector<int> pos_vertices_flip = pickup_a_few_numbers(Nv0, nv_flip, random);
  assert(pos_vertices_flip.size()==nv_flip);

  //set starting point (do not add non-interacting vertices)
  submatrix.init_update(std::vector<itime_vertex>());

  //perform actual updates
  for (int i_update=0; i_update<nv_flip; ++i_update) {
    weight_rat *= spin_flip_step(submatrix, Uijkl, random, pos_vertices_flip[i_update]);
  }

  //update A^{-1}
  submatrix.finalize_update();

  return weight_rat;
};

template<typename T>
template<typename R>
T VertexUpdateManager<T>::spin_flip_step(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random, int pos_vertex) {
  T det_rat_A, f_rat, U_rat;

  std::vector<int> pos_vertices_tmp(1), new_spins_tmp(1);
  pos_vertices_tmp[0] = pos_vertex;
  const int num_af_states = Uijkl.get_vertex(submatrix.itime_vertices()[pos_vertex].type()).num_af_states();

  if (num_af_states<=1) {
    return 1.0;
  }

  new_spins_tmp[0] = (int) (random()*num_af_states);
  while (new_spins_tmp[0] == submatrix.itime_vertices()[pos_vertex].af_state()) {
    new_spins_tmp[0] = (int) (random()*num_af_states);
  }

  boost::tie(det_rat_A,f_rat,U_rat) = submatrix.try_spin_flip(pos_vertices_tmp, new_spins_tmp);
  const T prob = det_rat_A*f_rat*U_rat;

  if (std::abs(prob)>random()) {
    //std::cout << "accepted " << std::endl;
    submatrix.perform_spin_flip(pos_vertices_tmp, new_spins_tmp);
    return prob;
  } else {
    //std::cout << "rejected " << std::endl;
    submatrix.reject_spin_flip();
    return 1.0;
  }
}

template<typename T>
template<typename R>
T VertexUpdateManager<T>::do_shift_update(SubmatrixUpdate<T>& submatrix, const general_U_matrix<T>& Uijkl, R& random, bool tune_step_size) {
  const int Nv0 = submatrix.pert_order();
  const int num_shift = std::min(Nv0, k_ins_max);

  std::vector<int> num_vertices_shift(num_shift, 1);

  T weight_rat = 1.0;

  std::vector<itime_vertex> new_vertices_all;
  const std::vector<int>& pos_vertices_shift = pickup_a_few_numbers(Nv0, num_shift, random);
  for (int i_shift=0; i_shift<num_shift; ++i_shift) {
    itime_vertex new_vertex = submatrix.itime_vertices()[pos_vertices_shift[i_shift]];
    new_vertex.set_time(mymod(new_vertex.time()+(random()-0.5)*shift_step_size, beta));
    new_vertex.set_non_interacting();
    new_vertices_all.push_back(new_vertex);
  }

  //set starting point
  assert(new_vertices_all.size()==num_shift);
  submatrix.init_update(new_vertices_all);

  //perform actual updates
  const double magic_number = std::pow(1.01, 1.0/num_shift);
  std::vector<int> pos_vertices_tmp(2), new_spins_tmp(2);
  for (int i_update=0; i_update<num_shift; ++i_update) {
    T det_rat_A, f_rat, U_rat;
    pos_vertices_tmp[0] = pos_vertices_shift[i_update];
    pos_vertices_tmp[1] = i_update+Nv0;
    assert(
      submatrix.itime_vertices()[pos_vertices_tmp[0]].type()==
      submatrix.itime_vertices()[pos_vertices_tmp[1]].type()
    );
    new_spins_tmp[0] = NON_INT_SPIN_STATE;
    new_spins_tmp[1] = submatrix.itime_vertices()[pos_vertices_tmp[0]].af_state();

    boost::tie(det_rat_A,f_rat,U_rat) = submatrix.try_spin_flip(pos_vertices_tmp, new_spins_tmp);
    assert(std::abs(U_rat-1.0)<1E-8);
    const T prob = det_rat_A*f_rat*U_rat;

    if (std::abs(prob)>random()) {
      submatrix.perform_spin_flip(pos_vertices_tmp, new_spins_tmp);
      if (tune_step_size) {
        shift_step_size = std::min(shift_step_size*magic_number, beta);
      }
      weight_rat *= prob;
    } else {
      submatrix.reject_spin_flip();
      if (tune_step_size) {
        shift_step_size = std::min(shift_step_size/magic_number, beta);
      }
    }
  } //i_update

  //update A^{-1}
  submatrix.finalize_update();

  return weight_rat;
};

template<typename T, typename SPLINE_G0, typename UnaryOperator, typename R>
T global_update_impl(const general_U_matrix<T>& Uijkl, const SPLINE_G0& spline_G0, itime_vertex_container& itime_vertices,
                   const UnaryOperator& op, R& random01, const T weight_old) {
  itime_vertex_container itime_vertices_new;
  for (int iv=0; iv<itime_vertices.size(); ++iv) {
    itime_vertices_new.push_back(op(itime_vertices[iv]));
  }

  const T weight_new = compute_weight(Uijkl, spline_G0, itime_vertices_new);
  const T prob = weight_new/weight_old;
  if (std::abs(prob)>random01()) {
    std::swap(itime_vertices_new, itime_vertices);
    return weight_new;
  } else {
    return weight_old;
  }
}

/*
struct GlobalOrbitalFlip {
  GlobalOrbitalFlip(int orb1, int orb2) : orb1_(orb1), orb2_(orb2);

  itime_vertex operator()(itime_vertex v) {
  }
  int orb1_, orb2_;
};

template<typename T, typename SPLINE_G0, typename UnaryOperator, typename R>
T global_updates(const general_U_matrix<T>& Uijkl, const SPLINE_G0& spline_G0, itime_vertex_container& itime_vertices, R& random01) {

  T weight = compute_weight(Uijkl, spline_G0, itime_vertices);

  weight = global_update_impl(Uijkl, spline_G0, itime_vertices, op, random01, weight);
}
 */

#endif //IMPSOLVER_UPDATE_MANAGER_HPP
