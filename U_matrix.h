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


//Data structure for two-body interaction
#ifndef U_MATRIX_H
#define U_MATRIX_H

#include <algorithm>
#include "boost/multi_array.hpp"
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "types.h"
#include "util.h"
#include "alps/parameter.h"

typedef size_t vertex_t;
typedef size_t af_t;

class itime_vertex;
class all_type;
class density_type;
class non_density_type;
class non_density_type_in_window;

template<class T>
class vertex_definition
{
public:
   vertex_definition(size_t rank, size_t num_af_states, std::vector<spin_t>& flavors, std::vector<size_t>& sites, T Uval, boost::multi_array<T,2>& alpha_af_rank, int id)
           : rank_(rank), num_af_states_(num_af_states), flavors_(flavors), sites_(sites), Uval_(Uval), alpha_af_rank_(alpha_af_rank), id_(id) {
     assert(flavors_.size()==rank);
     assert(sites.size()==2*rank);
     assert(alpha_af_rank_.shape()[0]==num_af_states_);
     assert(alpha_af_rank_.shape()[1]==rank);
   };

   const std::vector<spin_t>& flavors() const {
     return flavors_;
   };

   const std::vector<size_t>& sites() const {
     return sites_;
   };

   double Uval() const {
     return Uval_;
   }

   size_t rank() const {
     return rank_;
   }

   size_t num_af_states() const {
     return num_af_states_;
   }

   T get_alpha(size_t af_state, size_t idx_rank) const {
     return alpha_af_rank_[af_state][idx_rank];
   }

   bool is_density_type() const {
     bool flag = true;
     for (int i_rank=0; i_rank<rank_; ++i_rank) {
       if (sites_[2*i_rank] != sites_[2*i_rank+1]) {
         flag = false;
         break;
       }
     }
     return flag;
   }

   int id() const {return id_;}

   //track changes of occupations when this vertex is applied to a vector.
   //group(flavor,site)
   void make_quantum_numbers(const std::vector<std::vector<int> >& group_map, int ndim_group_per_flavor) {
     occ_change.resize(num_af_states_);

     for (int i_af=0; i_af<num_af_states_; ++i_af) {
       occ_change[i_af].resize(0);
       int idx = 0;
       for (int i_rank=rank_-1; i_rank>=0; --i_rank) {
         int PH = 1; //PH=1 means that there is no need for PH inversion.
         const int site1 = sites_[2*i_rank];
         const int site2 = sites_[2*i_rank+1];
         const int flavor = flavors_[i_rank];
         //std::cout << "site1 " << site1 << std::endl;
         //std::cout << "site2 " << site2 << std::endl;
         //std::cout << "ndim " << ndim_group_per_flavor << std::endl;
         if (std::min(std::abs(get_alpha(i_af,i_rank)),std::abs(1-get_alpha(i_af,i_rank))) >0.5) {
           throw std::runtime_error("Please take the value of alpha sufficiently close to 0 or 1 when using quantum number conservation!");
         }
         if(site1==site2) {//if c^dagger c is not of density type.
           PH = std::abs(get_alpha(i_af, i_rank))<std::abs(1-get_alpha(i_af,i_rank)) ? 1 : -1;
         }
         if (PH==1) {
           //apply c first, then cdagger.
           occ_change[i_af].push_back(boost::make_tuple(group_map[flavor][site2]+ndim_group_per_flavor*flavor,flavor,-1));
           occ_change[i_af].push_back(boost::make_tuple(group_map[flavor][site1]+ndim_group_per_flavor*flavor,flavor,+1));
         } else {
           //apply cdagger first, then c.
           occ_change[i_af].push_back(boost::make_tuple(group_map[flavor][site1]+ndim_group_per_flavor*flavor,flavor,+1));
           occ_change[i_af].push_back(boost::make_tuple(group_map[flavor][site2]+ndim_group_per_flavor*flavor,flavor,-1));
         }
         idx += 2;
       }
       assert(occ_change[i_af].size()==2*rank_);
     }
     //for (int i_af=0; i_af<num_af_states_; ++i_af) {
       //std::cout << "i_af " << i_af << std::endl;
       //for (int i=0; i<2*rank_; ++i)
         //std::cout << " op " << i << " " << boost::get<0>(occ_change[i_af][i]) << " " << boost::get<2>(occ_change[i_af][i]) << std::endl;
     //}
   }

   void apply_occ_change(int i_af, std::valarray<int>& occ_state, std::valarray<int>& max_occ, std::valarray<int>& min_occ) const {
     assert(i_af<num_af_states_ && i_af>=0);
     for (std::vector<boost::tuple<int,int,int> >::const_iterator it=occ_change[i_af].begin(); it!=occ_change[i_af].end(); ++it) {
       const int group = boost::get<0>(*it);
       occ_state[group] += boost::get<2>(*it);
       max_occ[group] = std::max(occ_state[group],max_occ[group]);
       min_occ[group] = std::min(occ_state[group],min_occ[group]);
       assert(group<occ_state.size());
     }
   }

private:
   size_t rank_;
   std::vector<spin_t> flavors_;
   std::vector<size_t> sites_;
   std::vector<std::vector<boost::tuple<int,int,int> > >  occ_change;
   size_t num_af_states_;
   T Uval_;
   boost::multi_array<T,2> alpha_af_rank_;//first index addresses af spin state, second one addresses (cdagger c)
   int id_;
};



//Data structure for general two-body interactions for a multi-orbital cluster impurity problem
template<class T>
class general_U_matrix {
  public:
    typedef T value_type;

    general_U_matrix(const alps::Parameters &parms) :
            ns_(parms.value_or_default("SITES", 1)),
            nf_(parms.value_or_default("FLAVORS", 2))
    {
      if (!parms.defined("GENERAL_U_MATRIX_FILE")) {
        throw std::runtime_error("Error: GENERAL_U_MATRIX_FILE is not defined!");
      }
      std::string ufilename(parms["GENERAL_U_MATRIX_FILE"]);
      std::ifstream ifs(ufilename.c_str());
      size_t num_non_zero_;
      ifs >> num_nonzero_;

      //temporary
      size_t rank, num_af_states;
      T Uval_;
      std::complex<double> Uval_cmplx;
      std::vector<size_t> site_indices_;//site indices (ijkl)
      std::vector<spin_t> flavor_indices_;//flavor indices for c^dagger c
      boost::multi_array<T,2> alpha_;//the first index is auxially spin, the second denotes (ij) or (kl).
      boost::multi_array<std::complex<double>,2> alpha_cmplx;//tempolary

      for (unsigned int idx=0; idx<num_nonzero_; ++idx) {
        size_t itmp;
        ifs >> itmp >> rank >> num_af_states >> Uval_cmplx;
        Uval_ = mycast<T>(Uval_cmplx);
        assert(itmp==idx);
        assert(rank==2);

        site_indices_.resize(2*rank);
        flavor_indices_.resize(rank);
        alpha_.resize(boost::extents[num_af_states][rank]);
        alpha_cmplx.resize(boost::extents[num_af_states][rank]);

        for (size_t i_op=0; i_op<2*rank; ++i_op) {
          ifs >> site_indices_[i_op];//i, j, k, l
        }
        for (size_t i_rank=0; i_rank<rank; ++i_rank) {
          ifs >> flavor_indices_[i_rank];
        }
        for (size_t i_rank=0; i_rank<rank; ++i_rank) {
          for (size_t iaf=0; iaf<num_af_states; ++iaf) {
            ifs >> alpha_cmplx[iaf][i_rank];
          }
        }
        for (size_t i_rank=0; i_rank<rank; ++i_rank) {
          for (size_t iaf = 0; iaf < num_af_states; ++iaf) {
              alpha_[iaf][i_rank] = mycast<T>(static_cast<std::complex<double> >(alpha_cmplx[iaf][i_rank]));
          }
        }

        vertex_list.push_back(vertex_definition<T>(rank, num_af_states, flavor_indices_, site_indices_, Uval_, alpha_, idx));
      }

      find_non_density_vertices();
    }

    size_t n_vertex_type() const{return vertex_list.size();}
    spin_t nf()const {return nf_;}
    spin_t ns()const {return ns_;}

    vertex_definition<T>& get_vertex(size_t vertex_idx) {
      assert(vertex_idx<n_vertex_type());
      return vertex_list[vertex_idx];
    }

    const vertex_definition<T>& get_vertex(size_t vertex_idx) const {
      assert(vertex_idx<n_vertex_type());
      return vertex_list[vertex_idx];
    }

    const std::vector<vertex_definition<T> >& get_vertices() const {
      return vertex_list;
    }

    const std::vector<vertex_definition<T> >& get_vertices(all_type& pred) const {
      return vertex_list;
    }

    const std::vector<vertex_definition<T> >& get_vertices(non_density_type& pred) const {
      return non_density_vertices;
    }

    const std::vector<vertex_definition<T> >& get_vertices(density_type& pred) const {
      return density_vertices;
    }

    const std::vector<vertex_definition<T> >& get_vertices(non_density_type_in_window& pred) const {
        return non_density_vertices;
    }

    const std::vector<vertex_definition<T> >& get_non_density_vertex_defs() const {
      return non_density_vertices;
    }

    int num_vertex_type(all_type& pred) const {
      return n_vertex_type();
    }

    int num_vertex_type(non_density_type& pred) const {
        return non_density_vertices.size();
    }

    int num_vertex_type(density_type& pred) const {
      return density_vertices.size();
    }

    int num_vertex_type(non_density_type_in_window& pred) const {
        return non_density_vertices.size();
    }

    int num_density_vertex_type() const {
      return density_vertices.size();
    }


  private:
    unsigned int ns_, nf_, num_nonzero_;
    std::vector<vertex_definition<T> > vertex_list, non_density_vertices, density_vertices;
    //std::vector<int> non_density_vertices;
    std::vector<bool> is_density_type;

    void find_non_density_vertices() {
      is_density_type.resize(vertex_list.size());
      non_density_vertices.clear();
      for (int iv=0; iv<vertex_list.size(); ++iv) {
        is_density_type[iv] = vertex_list[iv].is_density_type();
        if (!is_density_type[iv]) {
          non_density_vertices.push_back(vertex_list[iv]);
        } else {
          density_vertices.push_back(vertex_list[iv]);
        }
      }
    }
 };

//to remember what vertices are on the imaginary time axis..
typedef struct itime_vertex {
public:
  itime_vertex()
             : vertex_type_(-1),
               af_state_(-1),
               time_(-1),
               rank_(-1),
               is_density_type_(false)
  {}

  itime_vertex(int vertex_type, int af_state, double time, int rank, bool is_density_type)
           : vertex_type_(vertex_type),
             af_state_(af_state),
             time_(time),
             rank_(rank),
             is_density_type_(is_density_type)
  {}

  int af_state() const { return af_state_; }
  int vertex_type() const {return vertex_type_;}
  int type() const {return vertex_type_;}
  int rank() const {return rank_;}
  double time() const {return time_;}
  void set_time(double new_time) {time_ = new_time;}
  bool is_density_type() const {return is_density_type_;}

private:
  int vertex_type_, af_state_, rank_;
  double time_;
  bool is_density_type_;
} itime_vertex;

inline bool operator<(const itime_vertex& v1, const itime_vertex& v2) {
  return (v1.time()<v2.time());
}
//template<class T>
//itime_vertex generate_itime_vertex(vertex_definition<T> vertex_df
//template<class V>
//int num_non_denisty_vertices(const V& itime_vertices) {
  //int r = 0;
  //const std::vector<bool>& tmp = Uijkl.get_is_denisty_type();
  //for (typename V::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    //if (tmp[it->type()]) ++r;
  //}
  //return r;
//}

class density_type : public std::unary_function<itime_vertex,bool> {
public:
    bool operator()(itime_vertex v) {return v.is_density_type();}

    double random_time(double random01, double beta) const {
      assert(random01>=0 && random01<=1);
      return random01*beta;
    }
};

class non_density_type : public std::unary_function<itime_vertex,bool> {
public:
    bool operator()(itime_vertex v) {return !v.is_density_type();}

    double random_time(double random01, double beta) const {
        assert(random01>=0 && random01<=1);
        return random01*beta;
    }
};

class non_density_type_in_window : public std::unary_function<itime_vertex,bool> {
public:
    non_density_type_in_window() : ts_(0), w_(0), t_small1_(0), t_small2_(0), t_large1_(1E+100), t_large2_(1E+100) {}

    non_density_type_in_window(double ts, double w, double beta) : ts_(ts), w_(w), beta_(beta) {
        assert(w<=beta);
        if (ts+w<=beta) {
            t_small1_ = t_small2_= ts;
            t_large1_ = t_large2_ = ts+w;
        } else {
            t_small1_ = ts;
            t_large1_ = beta;
            t_small2_ = 0.0;
            t_large2_ = w+ts-beta;
        }
        //std::cout << "t_s1,t_l1" << t_small1_ << " " << t_large1_ << std::endl;
        //std::cout << "t_s2,t_l2" << t_small2_ << " " << t_large2_ << std::endl;
    }

    bool operator()(itime_vertex v) const {
        const double t = v.time();
        return !v.is_density_type() && ((t_small1_<=t && t<=t_large1_) || (t_small2_<=t && t<=t_large2_));
    }

    double random_time(double random01, double beta)  const{
        assert(random01>=0 && random01<=1);
        assert(beta_==beta);
        assert(w_>0);
        double t = random01*w_+ts_;
        if (t>beta_) t -= beta_;
        return t;
    }

    double width() {
        assert(w_>0);
        return w_;
    }

private:
    double ts_, w_, beta_;
    double t_small1_, t_large1_;
    double t_small2_, t_large2_;
};

class all_type : public std::unary_function<itime_vertex,bool> {
public:
    bool operator()(itime_vertex v) {return true;}

    double random_time(double random01, double beta)  const{
        assert(random01>=0 && random01<=1);
        return random01*beta;
    }
};

template<class T, class R, class UnaryPredicate>
std::vector<itime_vertex> generate_valid_vertex_pair(const general_U_matrix<T>& Uijkl, const std::vector<boost::tuple<int,int> >& valid_pair,
  R& random01, double beta, const UnaryPredicate& pred) {
  const int n_vertices_add = 2;

  std::vector<itime_vertex> itime_vertices;
  itime_vertices.reserve(n_vertices_add);

  if (valid_pair.size()==0)
    return std::vector<itime_vertex>();

  int v1,v2;
  boost::tie(v1,v2) = valid_pair[static_cast<int>(random01()*valid_pair.size())];
  int vtypes[] = {v1, v2};

  for (int iv=0; iv<n_vertices_add; ++iv) {
    const double time = pred.random_time(random01(), beta);
    const int v_type = vtypes[iv];
    const int rank = Uijkl.get_vertices()[v_type].rank();
    const int af_state = random01()*Uijkl.get_vertices()[v_type].num_af_states();
    itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank, false));
    if(Uijkl.get_vertices()[v_type].is_density_type()) {
      throw std::logic_error("Error found density type vertex");
    }
  }
  return itime_vertices;
}

template<class T, class R, class P>
std::vector<itime_vertex> generate_valid_vertex_pair2(const general_U_matrix<T>& Uijkl, const std::pair<int,int> v_pair,
  R& random01, double beta, const P& normalized_prob_dist) {
  const int n_vertices_add = 2;

  std::vector<itime_vertex> itime_vertices;
  itime_vertices.reserve(n_vertices_add);

  int vtypes[] = {v_pair.first, v_pair.second};
  double times[2];
  times[0] = mymod(beta*random01(), beta);
  times[1] = mymod(times[0]+gen_rand_rejection_method(normalized_prob_dist, normalized_prob_dist(0.0), random01, beta), beta);
  assert(vtypes[0]<vtypes[1]);

  for (int iv=0; iv<n_vertices_add; ++iv) {
    const double time = times[iv];
    const int v_type = vtypes[iv];
    const int rank = Uijkl.get_vertices()[v_type].rank();
    const int af_state = random01()*Uijkl.get_vertices()[v_type].num_af_states();
    itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank, false));
    if(Uijkl.get_vertices()[v_type].is_density_type()) {
      throw std::logic_error("Error found density type vertex");
    }
  }
  return itime_vertices;
}


template<class T, class R, class UnaryPredicate>
std::vector<itime_vertex> generate_itime_vertices(const general_U_matrix<T>& Uijkl, R& random01, double beta, int n_vertices_add, UnaryPredicate pred) {
  std::vector<itime_vertex> itime_vertices;
  itime_vertices.reserve(n_vertices_add);

  const std::vector<vertex_definition<T> >& valid_vs = Uijkl.get_vertices(pred);
  const int n_valid_vs = valid_vs.size();
  if (n_valid_vs==0)
    return std::vector<itime_vertex>();

  for (int iv=0; iv<n_vertices_add; ++iv) {
    const double time = pred.random_time(random01(), beta);
    const int iv_rnd = static_cast<int>(random01()*n_valid_vs);
    const int v_type = valid_vs[iv_rnd].id();
    const int rank = valid_vs[iv_rnd].rank();
    const int af_state = static_cast<size_t>(random01()*valid_vs[iv_rnd].num_af_states());
    itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank, valid_vs[iv_rnd].is_density_type()));
  }
  return itime_vertices;
}

template<class UnaryPredicate>
int count_valid_vertex_pair(const std::vector<itime_vertex>& itime_vertices, const boost::multi_array<bool,2>& valid_pair_flag, const UnaryPredicate& window) {
  std::vector<int> vs_win; vs_win.reserve(itime_vertices.size());
  std::vector<int> pos; pos.reserve(itime_vertices.size());

  int idx = 0;
  for (std::vector<itime_vertex>::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    if (window(*it) && !it->is_density_type()) {
      vs_win.push_back(it->type());
      pos[vs_win.size()-1] = idx;
    }
    ++idx;
  }

  const int N = vs_win.size();
  int N_valid_pair = 0;
  for (int iv=0; iv<N; ++iv) {
    for (int iv2=0; iv2<N; ++iv2) {
      if (iv<=iv2)
        continue;

      if (valid_pair_flag[vs_win[iv]][vs_win[iv2]])
        ++N_valid_pair;
    }
  }

  return N_valid_pair;
}

template<class UnaryPredicate, class R>
std::vector<int>
pick_up_valid_vertex_pair(const std::vector<itime_vertex>& itime_vertices, const boost::multi_array<bool,2>& valid_pair_flag, const UnaryPredicate& window, R& random01, int& N_valid_pair) {
  std::vector<int> vs_win; vs_win.reserve(itime_vertices.size());
  std::vector<int> pos; pos.reserve(itime_vertices.size());

  int idx = 0;
  for (std::vector<itime_vertex>::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
    if (window(*it) && !it->is_density_type()) {
      vs_win.push_back(it->type());
      pos.push_back(idx);
    }
    ++idx;
  }

  const int N = vs_win.size();
  std::vector<std::pair<int,int> > pos_pair;
  N_valid_pair = 0;
  for (int iv=0; iv<N; ++iv) {
    for (int iv2=0; iv2<N; ++iv2) {
      if (iv<=iv2)
        continue;

      if (valid_pair_flag[vs_win[iv]][vs_win[iv2]]) {
        ++N_valid_pair;
        pos_pair.push_back(std::make_pair(pos[iv],pos[iv2]));
      }
    }
  }

  if (N_valid_pair==0)
    return std::vector<int>();

  assert(pos_pair.size()==N_valid_pair);
  //const int selected = static_cast<int>(random01()*N_valid_pair);
  //std::cout << " debug " << selected << " " << N_valid_pair << " " << pos_pair.size() << std::endl;
  //assert(selected<N_valid_pair);
  std::pair<int,int> tmp = pos_pair[static_cast<int>(random01()*N_valid_pair)];
  assert(!itime_vertices[tmp.first].is_density_type());
  assert(!itime_vertices[tmp.second].is_density_type());
  std::vector<int> r(2);
  r[0] = tmp.first;
  r[1] = tmp.second;
  return r;
}

template<class P, class R>
std::pair<int,int>
pick_up_valid_vertex_pair2(const std::vector<itime_vertex>& itime_vertices, std::pair<int,int> v_pair,
                           double beta, P& p, R& random01, int& N_valid_pair, double& F) {

  std::vector<itime_vertex> v1, v2;
  std::vector<int> pos_v1, pos_v2;
  std::vector<double> prob;

  const int Nv = itime_vertices.size();
  for (int iv=0; iv<Nv; ++iv) {
    if (itime_vertices[iv].type()==v_pair.first) {
      v1.push_back(itime_vertices[iv]);
      pos_v1.push_back(iv);
    } else if (itime_vertices[iv].type()==v_pair.second) {
      v2.push_back(itime_vertices[iv]);
      pos_v2.push_back(iv);
    }
  }

  N_valid_pair = v1.size()*v2.size();
  if (N_valid_pair==0) {
    return std::make_pair(-1,-1);
  }

  assert(v1.size()==pos_v1.size());
  assert(v2.size()==pos_v2.size());
  assert(v1.size()==v2.size());
  if (v1.size()!=v2.size())
    throw std::logic_error("v1.size() != v2.size()");
  prob.resize(0); prob.reserve(v1.size()*v2.size());
  for (int iv1=0; iv1<v1.size(); ++iv1) {
    for (int iv2=0; iv2<v2.size(); ++iv2) {
      assert(itime_vertices[pos_v2[iv2]].type()==v_pair.second);
      assert(itime_vertices[pos_v1[iv1]].type()==v_pair.first);
      prob.push_back(p.bare_value(
          mymod(itime_vertices[pos_v2[iv2]].time()-itime_vertices[pos_v1[iv1]].time(), beta)
      ));
    }
  }

  F = std::accumulate(prob.begin(), prob.end(), 0.0);
  const int idx = boost::random::discrete_distribution<>(prob.begin(), prob.end())(random01.engine());
  assert(idx/v2.size()<v1.size());
  assert(idx%v2.size()<v2.size());
  return std::make_pair(pos_v1[idx/v2.size()], pos_v2[idx%v2.size()]);
}


template<class R, class UnaryPredicate>
std::vector<int> pick_up_itime_vertices(const std::vector<itime_vertex>& itime_vertices,
                                              R& random01,
                                              int n_vertices_rem,
                                              UnaryPredicate pred) {
  const int n_active_vertices = std::count_if(itime_vertices.begin(), itime_vertices.end(), pred);
  if (n_active_vertices<n_vertices_rem)
    return std::vector<int>();

  std::vector<int> pos(n_active_vertices);
  int idx=0;
  for (int iv=0; iv<itime_vertices.size(); ++iv) {
    if (pred(itime_vertices[iv])) {
      assert(idx<pos.size());
      pos[idx] = iv;
      ++idx;
    }
  }
  assert(idx==n_active_vertices);

  const std::vector<int>& indices = pickup_a_few_numbers(n_active_vertices, n_vertices_rem, random01);
  std::vector<int> indices2(n_vertices_rem);
  for (int i=0; i<n_vertices_rem; ++i) {
    assert(indices[i]<pos.size());
    indices2[i] = pos[indices[i]];
  }

#ifndef NDEBUG
  for (int i=0; i<n_vertices_rem; ++i) {
    assert(indices2[i]<itime_vertices.size());
    assert(pred(itime_vertices[indices2[i]]));
  }
#endif

  return indices2;
};

//void
//find_valid_pair_multi_vertex_update(const std::vector<std::vector<std::valarray<int> > >& quantum_numbers, std::vector<boost::tuple<int,int,int,int> >& v_pair, boost::multi_array<bool,4>& v_pair_flag);
template<class T>
void
find_valid_pair_multi_vertex_update(const std::vector<vertex_definition<T> >& vertex_defs, const std::vector<std::vector<std::valarray<int> > >& quantum_numbers, std::vector<std::pair<int,int> >& v_pair, boost::multi_array<bool,2>& v_pair_flag) {
  v_pair.resize(0);
  v_pair_flag.resize(boost::extents[quantum_numbers.size()][quantum_numbers.size()]);

  std::valarray<int> qn(quantum_numbers[0][0].size());
  const int i_af1=0, i_af2=0;
  for (int v1=0; v1<quantum_numbers.size(); ++v1) {
    for (int v2=0; v2<quantum_numbers.size(); ++v2) {
      qn = quantum_numbers[v1][i_af1] + quantum_numbers[v2][i_af2];
      bool flag = is_all_zero<int>(qn) && (!vertex_defs[v1].is_density_type()) && (!vertex_defs[v2].is_density_type());

      v_pair_flag[v1][v2] = flag;
      if (v1 < v2 && flag)
        v_pair.push_back(std::make_pair(v1, v2));
    }
  }
};

std::ostream &operator<<(std::ostream &os, const itime_vertex &v);
void print_vertices(std::ostream &os, const std::vector<itime_vertex> &v);

//U_MATRIX_H
#endif 
