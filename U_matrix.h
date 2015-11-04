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

#include "boost/multi_array.hpp"

#include "types.h"
#include "util.h"
#include "alps/parameter.h"

//Data structure for repulsion for density density bands.
class U_matrix{
public:
  U_matrix(const alps::Parameters &parms) :
    ns_(parms.value_or_default("SITES", 1)),
    nf_(parms.value_or_default("FLAVORS", 2)),
      n_nonzero_(0), mu_shift_(0)
  {
    val_ = new double[nf_*nf_];
    for(unsigned i=0; i<nf_*nf_; ++i)
      val_[i]=0; //default: non-interacting.
    if(parms.defined("U_MATRIX")){
      std::string ufilename(parms["U_MATRIX"]);
      std::ifstream u_file(ufilename.c_str());
      assert(u_file.is_open());
      int i;
      int j;
      double U_ij;
      while(u_file>>i>>j>>U_ij){
        operator()(i,j)=U_ij;
      }
    } else if (nf_==1) {
      //special case: only 1 orbital
      operator()(0,0)=(double)(parms["U"]);
    }else if (nf_==2) {
      //your ordinary two site problem
      assert(parms.defined("U"));
      double U=(double)(parms["U"]);
      operator()(0,0)=0; operator()(1,1)=0;
      operator()(0,1)=U; operator()(1,0)=U;
    } else {
      assert(parms.defined("U") && parms.defined("J"));
      double U=(double)(parms["U"]);
      double J=(double)(parms["J"]);
      double Uprime = parms.value_or_default("U'", U-2*J);
      assemble(U, Uprime, J);
    }
    for (unsigned i=0; i<nf_*nf_; ++i)
      if (val_[i]!=0)
        n_nonzero_++;
    for (unsigned i=0; i<nf_; ++i)
      mu_shift_ += operator()(i,0);
    mu_shift_ /= 2;
  }


  void assemble(const double U, const double Uprime, const double J){
    //this implements the U matrix for the special case of n_flavor/2 degenerate bands
    assert(ns_==1);
    assert(nf_%2==0);
    for(spin_t i=0;i<nf_;i+=2){
      operator()(i  , i  ) = 0; //Pauli
      operator()(i+1, i+1) = 0; //Pauli
      operator()(i  , i+1) = U; //Hubbard repulsion same band
      operator()(i+1, i  ) = U; //Hubbard repulsion same band
      for(spin_t j=0; j<nf_; j+=2){
        if(j==i)
          continue;
        operator()(i  ,j  ) = Uprime-J; //Hubbard repulsion interband same spin
        operator()(i+1,j+1) = Uprime-J; //Hubbard repulsion interband same spin
        operator()(i  ,j+1) = Uprime; //Hubbard repulsion interband opposite spin (this used to be '+J', the rest of the world uses '-J' -> changed to be consistent).
        operator()(i+1,j  ) = Uprime; //Hubbard repulsion interband opposite spin
      }
    }
  }

  ~U_matrix(){
    delete[] val_;
  }

  double &operator()(spin_t flavor_i, spin_t flavor_j){
    return val_[flavor_i*nf_+flavor_j];
    }

  const double &operator() (spin_t flavor_i, spin_t flavor_j)const {
    return val_[flavor_i*nf_+flavor_j];
  }

  spin_t nf()const {return nf_;}
  spin_t ns()const {return ns_;}
  double mu_shift() const { return mu_shift_; }

  inline int n_nonzero() const{return n_nonzero_;}

private:
  double *val_;
  size_t ns_;
  size_t nf_;
  int n_nonzero_;
  double mu_shift_;
};

typedef size_t vertex_t;
typedef size_t af_t;

template<class T>
class vertex_definition
 {
 public:
    vertex_definition(size_t rank, size_t num_af_states, std::vector<spin_t>& flavors, std::vector<size_t>& sites, T Uval, boost::multi_array<T,2>& alpha_af_rank)
            : rank_(rank), num_af_states_(num_af_states), flavors_(flavors), sites_(sites), Uval_(Uval), alpha_af_rank_(alpha_af_rank) {
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

 private:
    size_t rank_;
    std::vector<spin_t> flavors_;
    std::vector<size_t> sites_;
    size_t num_af_states_;
    T Uval_;
    boost::multi_array<T,2> alpha_af_rank_;//first index addresses af spin state, second one addresses (cdagger c)
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
            //ifs >> alpha_[iaf][i_rank];
            ifs >> alpha_cmplx[iaf][i_rank];
            //std::cout << " i_rank, i_af " << alpha_cmplx[iaf][i_rank] << std::endl;
          }
        }
        for (size_t i_rank=0; i_rank<rank; ++i_rank) {
          for (size_t iaf = 0; iaf < num_af_states; ++iaf) {
              alpha_[iaf][i_rank] = mycast<T>(static_cast<std::complex<double> >(alpha_cmplx[iaf][i_rank]));
          }
        }

        vertex_list.push_back(vertex_definition<T>(rank, num_af_states, flavor_indices_, site_indices_, Uval_, alpha_));
      }
    }

    size_t n_vertex_type() const{return vertex_list.size();}
    spin_t nf()const {return nf_;}
    spin_t ns()const {return ns_;}

    const vertex_definition<T>& get_vertex(size_t vertex_idx) const {
      assert(vertex_idx<n_vertex_type());
      return vertex_list[vertex_idx];
    }

    std::vector<vertex_definition<T> > get_vertices() const {
      return vertex_list;
    }

  private:
    unsigned int ns_, nf_, num_nonzero_;
    std::vector<vertex_definition<T> > vertex_list;
 };

 //to remember what vertices are on the imaginary time axis..
 typedef struct itime_vertex {
 public:
   itime_vertex()
             : vertex_type_(-1),
               af_state_(-1),
               time_(-1),
               rank_(-1) {}

   itime_vertex(size_t vertex_type, size_t af_state, double time, size_t rank)
           : vertex_type_(vertex_type),
             af_state_(af_state),
             time_(time),
             rank_(rank) {}

   size_t af_state() const { return af_state_; }
   size_t vertex_type() const {return vertex_type_;}
   size_t type() const {return vertex_type_;}
   size_t rank() const {return rank_;}
   double time() const {return time_;}

 private:
   size_t vertex_type_, af_state_, rank_;
   double time_;
 } itime_vertex;

template<class T, class R>
std::vector<itime_vertex> generate_vertices(const general_U_matrix<T>& Uijkl, R& random01, double beta, int n_vertices_add) {
  std::vector<itime_vertex> itime_vertices;
  itime_vertices.reserve(n_vertices_add);
  for (int iv=0; iv<n_vertices_add; ++iv) {
    const double time = beta * random01();
    const size_t v_type = static_cast<size_t>(random01() * Uijkl.n_vertex_type());
    const vertex_definition<T> new_vertex_type = Uijkl.get_vertex(v_type);
    const size_t rank = new_vertex_type.rank();
    const size_t af_state = static_cast<size_t>(random01() * new_vertex_type.num_af_states());
    itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank));
  }
  return itime_vertices;
}

std::ostream &operator<<(std::ostream &os, const itime_vertex &v);

std::ostream &operator<<(std::ostream &os, const U_matrix &U);
//U_MATRIX_H
#endif 
