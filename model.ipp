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

template<class TYPES>
std::pair<typename TYPES::M_TYPE, typename TYPES::M_TYPE>
InteractionExpansion<TYPES>::try_add(size_t n_vertices_add, std::vector<itime_vertex>& new_vertices)
{
  assert(n_vertices_add==new_vertices.size());

  //work array
  add_helper.clear();

  //add vertices one by one
  M_TYPE prod_Uval_wm = 1.0;
  for (size_t iv=0; iv<n_vertices_add; ++iv) {
    const size_t v_type = new_vertices[iv].type();
    const vertex_definition<M_TYPE> new_vertex_type = Uijkl.get_vertex(v_type);
    const double time = new_vertices[iv].time();
    const size_t rank = new_vertex_type.rank();
    const size_t af_state = new_vertices[iv].af_state();
    prod_Uval_wm *= -new_vertex_type.Uval();

    for (size_t i_rank=0; i_rank<rank; ++i_rank) {
      const size_t flavor_rank = new_vertex_type.flavors()[i_rank];

      operator_time op_t(time, -i_rank);
      M[flavor_rank].creators().push_back(creator(flavor_rank, new_vertex_type.sites()[2*i_rank], op_t, n_matsubara));
      M[flavor_rank].annihilators().push_back(annihilator(flavor_rank, new_vertex_type.sites()[2*i_rank+1], op_t, n_matsubara));
      M[flavor_rank].alpha_push_back(new_vertex_type.get_alpha(af_state, i_rank));
      M[flavor_rank].vertex_info().push_back(std::pair<vertex_t,size_t>(v_type,i_rank));

      ++(add_helper.num_new_rows[flavor_rank]);
    }
    itime_vertices.push_back(new_vertices[iv]);
  }

  //fast update for each flavor
  M_TYPE lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (add_helper.num_new_rows[flavor]>0) {
      lambda_prod *= fastupdate_up(flavor, true, add_helper.num_new_rows[flavor]); // true means compute_only_weight
    }
  }
  add_helper.det_rat_ = lambda_prod;
  return std::pair<M_TYPE,M_TYPE>(prod_Uval_wm*lambda_prod, lambda_prod);
}


template<class TYPES>
void InteractionExpansion<TYPES>::perform_add(size_t n_vertices_add)
{
  M_TYPE det_rat = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (add_helper.num_new_rows[flavor]>0) {
      det_rat *= fastupdate_up(flavor,false,add_helper.num_new_rows[flavor]);
    }
  }
  assert(std::abs(det_rat/add_helper.det_rat_-1.0)<1E-8);
  M.sanity_check(itime_vertices, Uijkl);
}


template<class TYPES>
void InteractionExpansion<TYPES>::reject_add(size_t n_vertices_add)
{
  //get rid of the operators
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    for (int i=0; i<add_helper.num_new_rows[flavor]; ++i) {
      M[flavor].pop_back_op();
    }
  }
  for(int i=0; i<n_vertices_add; ++i) {
    itime_vertices.pop_back();
  }
  M.sanity_check(itime_vertices, Uijkl);
}

template<class TYPES>
std::pair<typename TYPES::M_TYPE,typename TYPES::M_TYPE>
InteractionExpansion<TYPES>::try_remove(const std::vector<int>& vertices_nr)
{
  //boost::timer::cpu_timer timer;
  //get weight
  //figure out and remember which rows (columns) are to be removed
  // lists of rows and cols will be sorted in ascending order
  const size_t nv_old = itime_vertices.size();
  M.sanity_check(itime_vertices, Uijkl);
  remove_helper.clear();
  typename TYPES::M_TYPE prod_Uval = 1.0;
  for (size_t iv=0; iv<vertices_nr.size(); ++iv) {
    vertex_definition<M_TYPE> vertex_def = Uijkl.get_vertex(itime_vertices[vertices_nr[iv]].type());
    prod_Uval *= -vertex_def.Uval();
    for (size_t i_rank=0; i_rank<vertex_def.rank(); ++i_rank) {
      const size_t flavor = vertex_def.flavors()[i_rank];
      int r = M[flavor].find_row_col(itime_vertices[vertices_nr[iv]].time(), itime_vertices[vertices_nr[iv]].type(), i_rank);
      assert(r<num_rows(M[flavor].matrix()));
      remove_helper.rows_cols_removed[flavor].push_back(r);
    }
  }
  remove_helper.sort_rows_cols();

  M_TYPE lambda_prod = 1.0;
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (remove_helper.rows_cols_removed[flavor].size()>0) {
      lambda_prod *= fastupdate_down(remove_helper.rows_cols_removed[flavor], flavor, true);  // true means compute_only_weight
    }
  }
  remove_helper.det_rat_ = lambda_prod;
  assert(itime_vertices.size()==nv_old);
  return std::pair<typename TYPES::M_TYPE,typename TYPES::M_TYPE>(lambda_prod/prod_Uval,lambda_prod);
}


template<class TYPES>
void InteractionExpansion<TYPES>::perform_remove(const std::vector<int>& vertices_nr)
{
  //boost::timer::cpu_timer timer;
  //perform fastupdate down
  M_TYPE det_rat = 1.0;
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (remove_helper.rows_cols_removed[flavor].size()>0) {
      //remove rows and columns
      det_rat *= fastupdate_down(remove_helper.rows_cols_removed[flavor], flavor, false);  // false means really perform, not only compute weight
      //get rid of operators
      for (int iop=0; iop<remove_helper.rows_cols_removed[flavor].size(); ++iop) {
        M[flavor].pop_back_op();
      }
    }
  }
  assert(std::abs(det_rat/remove_helper.det_rat_-1.0)<1E-8);
  //get rid of vertex list entries.
  remove_elements_from_vector(itime_vertices, vertices_nr);
  M.sanity_check(itime_vertices, Uijkl);
}


template<class TYPES>
void InteractionExpansion<TYPES>::reject_remove()
{
  //do nothing
  return;
}

//This is specialized for pair hopping and spin flip.
template<class TYPES>
typename TYPES::M_TYPE InteractionExpansion<TYPES>::try_shift(int idx_vertex, double new_time) {
  assert(idx_vertex<itime_vertices.size());

  itime_vertex& v = itime_vertices[idx_vertex];
  v.set_time(new_time);
  shift_helper.find_rows_cols_set_time(v.rank(), v.type(), Uijkl.get_vertex(v.type()).flavors(), shift_helper.get_old_time(), v.time(), M);

  typename TYPES::M_TYPE lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (shift_helper.num_rows_cols_updated[flavor]==0)
      continue;

    const int Nv = M[flavor].matrix().num_cols();
    for (int i=0; i<shift_helper.rows_cols_updated[flavor].size(); ++i) {
      const int idx = shift_helper.rows_cols_updated[flavor][i];
      const int n_rows_cols_updated = shift_helper.rows_cols_updated[flavor].size();
      assert(idx<M[flavor].creators().size());

      operator_time cr_op_time = M[flavor].creators()[idx].t(); cr_op_time.set_time(v.time());
      M[flavor].creators()[idx].set_time(cr_op_time);

      operator_time ann_op_time = M[flavor].annihilators()[idx].t(); ann_op_time.set_time(v.time());
      M[flavor].annihilators()[idx].set_time(ann_op_time);
    }

    //actual fast update
    fastupdate_shift_init(flavor, shift_helper.rows_cols_updated[flavor]);
    M_TYPE tmp = fastupdate_shift(flavor, shift_helper.rows_cols_updated[flavor], true);
    lambda_prod *= tmp;
  }
  return lambda_prod;
}

template<class TYPES>
void InteractionExpansion<TYPES>::perform_shift(int idx_vertex) {
  assert(idx_vertex<itime_vertices.size());

  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (shift_helper.num_rows_cols_updated[flavor]==0)
      continue;

    fastupdate_shift(flavor, shift_helper.rows_cols_updated[flavor], false);
    fastupdate_shift_finalize(flavor, shift_helper.rows_cols_updated[flavor]);
  }
}

template<class TYPES>
void InteractionExpansion<TYPES>::reject_shift(int idx_vertex) {
  for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
    if (shift_helper.num_rows_cols_updated[flavor]==0)
      continue;

    fastupdate_shift_finalize(flavor, shift_helper.rows_cols_updated[flavor]);

    itime_vertices[idx_vertex].set_time(shift_helper.get_old_time());
    for (int i=0; i<shift_helper.rows_cols_updated[flavor].size(); ++i) {
      const int idx = shift_helper.rows_cols_updated[flavor][i];

      operator_time cr_op_time = M[flavor].creators()[idx].t();
      assert(cr_op_time.time()==shift_helper.get_new_time());
      cr_op_time.set_time(shift_helper.get_old_time());
      M[flavor].creators()[idx].set_time(cr_op_time);

      operator_time ann_op_time = M[flavor].annihilators()[idx].t();
      assert(ann_op_time.time()==shift_helper.get_new_time());
      ann_op_time.set_time(shift_helper.get_old_time());
      M[flavor].annihilators()[idx].set_time(ann_op_time);
    }
  }
}

/*
template<class TYPES>
void InteractionExpansion<TYPES>::try_spin_flip(int iv, int new_af_state) {
  itime_vertex& v = itime_vertices[iv];
  const vertex_definition<M_TYPE> v_def = Uijkl.get_vertex(v.type());

  //change status
  v.set_af_state(new_af_state);
  std::vector<bool> updated(n_flavors, false);
  for (int i_rank=0; i_rank<v_def.rank(); ++i_rank) {
    const int flavor_rank = v_def.flavors()[i_rank];
    int pos = M[flavor_rank].find_row_col(v.time(), v.type(), i_rank);
    if (M[flavor_rank].bare_alpha_at(pos)!=v_def.get_alpha(new_af_state,i_rank)) {
      M[flavor_rank].set_alpha(pos, v_def.get_alpha(new_af_state,i_rank));
      updated[flavor_rank] = true;
    }
  }

  //update M
  double det_rat = 1.0;
  for (int flavor=0; flavor<n_flavors; ++flavor) {

  }
}
 */
