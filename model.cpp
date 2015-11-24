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

std::pair<double,double> HubbardInteractionExpansion::try_add(fastupdate_add_helper& helper, size_t n_vertices_add, std::vector<itime_vertex>& new_vertices)
{
  assert(n_vertices_add==new_vertices.size());

  //work array
  helper.clear();

  //add vertices one by one
  double prod_Uval_wm = 1.0;
  for (size_t iv=0; iv<n_vertices_add; ++iv) {
    const size_t v_type = new_vertices[iv].type();
    const vertex_definition<GTYPE> new_vertex_type = Uijkl.get_vertex(v_type);
    const double time = new_vertices[iv].time();
    const size_t rank = new_vertex_type.rank();
    const size_t af_state = new_vertices[iv].af_state();
    prod_Uval_wm *= -new_vertex_type.Uval();

    for (size_t i_rank=0; i_rank<rank; ++i_rank) {
      const size_t flavor_rank = new_vertex_type.flavors()[i_rank];

      M[flavor_rank].creators().push_back(creator(flavor_rank, new_vertex_type.sites()[2*i_rank], time, n_matsubara));
      M[flavor_rank].annihilators().push_back(annihilator(flavor_rank, new_vertex_type.sites()[2*i_rank+1], time, n_matsubara));
      M[flavor_rank].alpha_push_back(new_vertex_type.get_alpha(af_state, i_rank));
      M[flavor_rank].vertex_info().push_back(std::pair<vertex_t,size_t>(v_type,i_rank));

      ++(helper.num_new_rows[flavor_rank]);
    }
    itime_vertices.push_back(new_vertices[iv]);
  }

  //fast update for each flavor
  double lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]>0) {
      lambda_prod *= fastupdate_up(flavor, true, helper.num_new_rows[flavor]); // true means compute_only_weight
    }
  }
  helper.det_rat_ = lambda_prod;
  return std::pair<double,double>(prod_Uval_wm*lambda_prod, lambda_prod);
}


void HubbardInteractionExpansion::perform_add(fastupdate_add_helper& helper, size_t n_vertices_add)
{
  double det_rat = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]>0) {
      det_rat *= fastupdate_up(flavor,false,helper.num_new_rows[flavor]);
    }
  }
  assert(std::abs(det_rat/helper.det_rat_-1)<1E-8);
  M.sanity_check(itime_vertices);
}


void HubbardInteractionExpansion::reject_add(fastupdate_add_helper& helper, size_t n_vertices_add)
{
  //get rid of the operators
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    for (int i=0; i<helper.num_new_rows[flavor]; ++i) {
      M[flavor].pop_back_op();
    }
  }
  for(int i=0; i<n_vertices_add; ++i) {
    itime_vertices.pop_back();
  }
  M.sanity_check(itime_vertices);
}


std::pair<double,double> HubbardInteractionExpansion::try_remove(const std::vector<int>& vertices_nr, fastupdate_remove_helper& helper)
{
  //get weight
  //figure out and remember which rows (columns) are to be removed
  // lists of rows and cols will be sorted in ascending order
  const size_t nv_old = itime_vertices.size();
  M.sanity_check(itime_vertices);
  helper.clear();
  double prod_Uval = 1.0;
  for (size_t iv=0; iv<vertices_nr.size(); ++iv) {
    vertex_definition<GTYPE> vertex_def = Uijkl.get_vertex(itime_vertices[vertices_nr[iv]].type());
    prod_Uval *= -vertex_def.Uval();
    for (size_t i_rank=0; i_rank<vertex_def.rank(); ++i_rank) {
      const size_t flavor = vertex_def.flavors()[i_rank];
      int r = M[flavor].find_row_col(itime_vertices[vertices_nr[iv]].time(), itime_vertices[vertices_nr[iv]].type(), i_rank);
      assert(r<num_rows(M[flavor].matrix()));
      helper.rows_cols_removed[flavor].push_back(r);
    }
  }
  helper.sort_rows_cols();

  GTYPE lambda_prod = 1.0;
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.rows_cols_removed[flavor].size()>0) {
      lambda_prod *= fastupdate_down(helper.rows_cols_removed[flavor], flavor, true);  // true means compute_only_weight
    }
  }
  helper.det_rat_ = lambda_prod;
  assert(itime_vertices.size()==nv_old);
  //double r1 = vertices_nr.size()== 1 ?
      //permutation(itime_vertices.size(),vertices_nr.size()) :
      //permutation(std::count_if(itime_vertices.begin(), itime_vertices.end(), helper.op), vertices_nr.size());
  //double r2 = vertices_nr.size()== 1 ?
      //pow(-beta*Uijkl.n_vertex_type(),(double)vertices_nr.size()) :
      //pow(-helper.op.width()*Uijkl.num_vertex_type(helper.op),(double)vertices_nr.size());
  return std::pair<double,double>(lambda_prod/prod_Uval,lambda_prod);
}


void HubbardInteractionExpansion::perform_remove(const std::vector<int>& vertices_nr, fastupdate_remove_helper& helper)
{
  //perform fastupdate down
  double det_rat = 1.0;
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.rows_cols_removed[flavor].size()>0) {
      //remove rows and columns
      det_rat *= fastupdate_down(helper.rows_cols_removed[flavor], flavor, false);  // false means really perform, not only compute weight
      //get rid of operators
      for (int iop=0; iop<helper.rows_cols_removed[flavor].size(); ++iop) {
        M[flavor].pop_back_op();
      }
    }
  }
  assert(std::abs(det_rat/helper.det_rat_-1)<1E-8);
  //get rid of vertex list entries.
  remove_elements_from_vector(itime_vertices, vertices_nr);
  M.sanity_check(itime_vertices);
}


void HubbardInteractionExpansion::reject_remove(fastupdate_remove_helper& helper)
{
  //do nothing
  return;
}

double HubbardInteractionExpansion::try_shift(int idx_vertex, double new_time) {
  assert(idx_vertex<itime_vertices.size());

  itime_vertex& v = itime_vertices[idx_vertex];
  shift_helper.old_time = v.time();
  v.set_time(new_time);
  shift_helper.find_rows_cols_set_time(v.rank(), v.type(), Uijkl.get_vertex(v.type()).flavors(), shift_helper.old_time, v.time(), M);

  double lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (shift_helper.num_rows_cols_updated[flavor]==0)
      continue;

    //std::cout << "Flavor " << flavor << std::endl;
    for (int i=0; i<shift_helper.rows_cols_updated[flavor].size(); ++i) {
      const int idx = shift_helper.rows_cols_updated[flavor][i];
      assert(idx<M[flavor].creators().size());
      M[flavor].creators()[idx].set_time(v.time());
      M[flavor].annihilators()[idx].set_time(v.time());
    }

    //actual fast update
    shift_helper.M_old[flavor] = M[flavor].matrix();
    lambda_prod *= fastupdate_shift(flavor, shift_helper.rows_cols_updated[flavor]);
  }
  return lambda_prod;
}

void HubbardInteractionExpansion::perform_shift(int idx_vertex) {
  assert(idx_vertex<itime_vertices.size());

  //itime_vertex& v = itime_vertices[idx_vertex];
  //std::vector<int> rows_cols_updated(v.rank());
  //for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    //if (shift_helper.num_rows_cols_updated[flavor]==0)
      //continue;
    //fastupdate_shift(flavor, shift_helper.rows_cols_updated[flavor], false);
  //}
}

void HubbardInteractionExpansion::reject_shift(int idx_vertex) {
  for (spin_t flavor = 0; flavor < n_flavors; ++flavor) {
    if (shift_helper.num_rows_cols_updated[flavor] == 0)
      continue;

    //swap_cols_rows(M[flavor].matrix(), shift_helper.swap_list[flavor].rbegin(), shift_helper.swap_list[flavor].rend());
    //M[flavor].swap_ops2(shift_helper.swap_list[flavor].rbegin(), shift_helper.swap_list[flavor].rend());
    M[flavor].matrix() = shift_helper.M_old[flavor];

    itime_vertices[idx_vertex].set_time(shift_helper.old_time);
    for (int i=0; i<shift_helper.rows_cols_updated[flavor].size(); ++i) {
      const int idx = shift_helper.rows_cols_updated[flavor][i];
      M[flavor].creators()[idx].set_time(shift_helper.old_time);
      M[flavor].annihilators()[idx].set_time(shift_helper.old_time);
    }
  }
  //std::cout << "debug rejected (2) " << std::endl;
}
