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

#include "types.h"
#include "util.h"
#include "alps/parameter.h"

typedef size_t vertex_t;
typedef size_t af_t;

class itime_vertex;
class all_type;
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

private:
   size_t rank_;
   std::vector<spin_t> flavors_;
   std::vector<size_t> sites_;
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
        assert(num_af_states==2);

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

    const vertex_definition<T>& get_vertex(size_t vertex_idx) const {
      assert(vertex_idx<n_vertex_type());
      return vertex_list[vertex_idx];
    }

    const std::vector<vertex_definition<T> >& get_vertices() const {
      return vertex_list;
    }

    //const std::vector<int>& get_non_density_vertices() const {
      //return non_density_vertices;
    //}

    const std::vector<vertex_definition<T> >& get_vertices(all_type& pred) const {
      return vertex_list;
    }

    const std::vector<vertex_definition<T> >& get_vertices(non_density_type& pred) const {
      return non_density_vertices;
    }

    const std::vector<vertex_definition<T> >& get_vertices(non_density_type_in_window& pred) const {
        return non_density_vertices;
    }

    int num_vertex_type(all_type& pred) const {
      return n_vertex_type();
    }

    int num_vertex_type(non_density_type& pred) const {
        return non_density_vertices.size();
    }

    int num_vertex_type(non_density_type_in_window& pred) const {
        return non_density_vertices.size();
    }

    //int num_vertex_type(non_density_type& pred) {
      //return non_density_vertices.size();
    //}
    //const std::vector<bool>& get_is_denisty_type() const {
      //return is_density_type;
    //}

  private:
    unsigned int ns_, nf_, num_nonzero_;
    std::vector<vertex_definition<T> > vertex_list, non_density_vertices;
    //std::vector<int> non_density_vertices;
    std::vector<bool> is_density_type;

    void find_non_density_vertices() {
      is_density_type.resize(vertex_list.size());
      non_density_vertices.clear();
      for (int iv=0; iv<vertex_list.size(); ++iv) {
        is_density_type[iv] = vertex_list[iv].is_density_type();
        if (!is_density_type[iv])
          non_density_vertices.push_back(vertex_list[iv]);
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
  bool is_density_type() const {return is_density_type_;}

private:
  int vertex_type_, af_state_, rank_;
  double time_;
  bool is_density_type_;
} itime_vertex;

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
};

class non_density_type : public std::unary_function<itime_vertex,bool> {
public:
    bool operator()(itime_vertex v) {return !v.is_density_type();}

    double random_time(double random01, double beta) {
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

    bool operator()(itime_vertex v) {
        const double t = v.time();
        return !v.is_density_type() && ((t_small1_<=t && t<=t_large1_) || (t_small2_<=t && t<=t_large2_));
    }

    double random_time(double random01, double beta) {
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

    double random_time(double random01, double beta) {
        assert(random01>=0 && random01<=1);
        return random01*beta;
    }
};

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


std::ostream &operator<<(std::ostream &os, const itime_vertex &v);

//U_MATRIX_H
#endif 
