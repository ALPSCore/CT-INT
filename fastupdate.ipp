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
#include <boost/numeric/bindings/blas/level2/ger.hpp>

#include <alps/numeric/matrix/algorithms.hpp>

#include "fastupdate_formula.h"

#include "util.h"

/// @brief This function performs the InteractionExpansion update of adding one vertex to the set of vertices. If the 
///        Green's function for the measurement has to be tracked it calls the function that does that.
/// @param vertex_nr specify which of the U n_up n_down you want to look at.
/// @param compute_only_weight Do not perform the fastupdate, but only return the Mnn entry (1/lambda), which is used for the acceptance weight.
/// @param track_green_matsubara Track Green function in Matsubara frequencies or imaginary time if required.

template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::fastupdate_up(const int flavor, bool compute_only_weight, size_t n_vertices_add)
{
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  unsigned int noperators = num_rows(M[flavor].matrix());
  //current size of M: number of vertices - n_vertices_add. We need to add the last vertex.
  //A pointer to the creator and annihilator is already stored in creators_ and annihilators_
  //at positions M[flavor].size() ... M[flavor].size()+n_vertices_add-1;
  alps::numeric::matrix<M_TYPE> Green0_n_n(n_vertices_add, n_vertices_add);
  alps::numeric::matrix<M_TYPE> Green0_n_j(n_vertices_add, noperators);
  alps::numeric::matrix<M_TYPE> Green0_j_n(noperators, n_vertices_add);
  for(unsigned int i=0;i<noperators;++i) {
    for (size_t iv=0; iv<n_vertices_add; ++iv) {
      Green0_n_j(iv,i) = mycast<M_TYPE>(green0_spline_for_M(flavor, noperators+iv, i));
    }
  }
  for(unsigned int i=0;i<noperators;++i){
    for (size_t iv=0; iv<n_vertices_add; ++iv) {
      Green0_j_n(i,iv) = mycast<M_TYPE>(green0_spline_for_M(flavor, i, noperators+iv));
    }
  }
  for (size_t iv2=0; iv2<n_vertices_add; ++iv2) {
    for (size_t iv=0; iv<n_vertices_add; ++iv) {
      Green0_n_n(iv, iv2) = mycast<M_TYPE>(green0_spline_for_M(flavor, noperators+iv, noperators+iv2));
    }
  }
  for (size_t iv=0; iv<n_vertices_add; ++iv) {
    Green0_n_n(iv, iv) -= mycast<M_TYPE>(M[flavor].alpha_at(noperators+iv));
  }

  //B: Green0_j_n
  //C: Green0_n_j
  //D: Green0_n_n
  //invA: M[flavor].matrix()
  if(compute_only_weight){
    return compute_det_ratio_up(Green0_j_n, Green0_n_j, Green0_n_n, M[flavor].matrix());
  } else {
    //boost::timer::cpu_timer timer;
    const M_TYPE r = compute_inverse_matrix_up2(Green0_j_n, Green0_n_j, Green0_n_n, M[flavor].matrix(), M[flavor].matrix());
    //std::cout << "matrix_up " << timer.elapsed().wall*1E-6 << std::endl;
    return r;
  }
}


///Fastupdate formulas, remove order by one (remove a vertex). If necessary
///also take track of the Green's function.
template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::fastupdate_down(const std::vector<size_t>& rows_cols_removed, const int flavor, bool compute_only_weight)
{
  using std::swap;
  using alps::numeric::column_view;
  using alps::numeric::vector;
  using alps::numeric::matrix;
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  const size_t n_vertices_remove = rows_cols_removed.size();

  if(compute_only_weight) {
    return compute_det_ratio_down(n_vertices_remove, rows_cols_removed, M[flavor].matrix());
  } else {
    std::vector<std::pair<size_t,size_t> > rows_cols_swap_list;
    M_TYPE det_rat = compute_inverse_matrix_down(n_vertices_remove,rows_cols_removed,M[flavor].matrix(),rows_cols_swap_list);

    assert(rows_cols_swap_list.size()==n_vertices_remove);
    for (size_t s=0; s<n_vertices_remove; ++s) {
      const size_t idx1 = rows_cols_swap_list[s].first;
      const size_t idx2 = rows_cols_swap_list[s].second;
      M[flavor].swap_ops(idx1, idx2);
    }
    return det_rat;
  }
}

template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::fastupdate_shift_init(const int flavor, const std::vector<int>& rows_cols_updated) {
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  const int num_rows_cols_updated = rows_cols_updated.size();
  const int noperators = num_rows(M[flavor].matrix());

  for (int i=0; i<num_rows_cols_updated; ++i) {
    M[flavor].swap_rows_cols(rows_cols_updated[num_rows_cols_updated-1-i], noperators-1-i);
  }
}

//VERY UGLY IMPLEMENTATION
template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::fastupdate_shift(const int flavor, const std::vector<int>& rows_cols_updated, bool compute_only_weight) {
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  const int num_rows_cols_updated = rows_cols_updated.size();
  const int noperators = num_rows(M[flavor].matrix());
  const int noperators_rest = noperators-num_rows_cols_updated;

  alps::numeric::matrix<M_TYPE> Green0_n_n(num_rows_cols_updated, num_rows_cols_updated);//S
  alps::numeric::matrix<M_TYPE> Green0_n_j(num_rows_cols_updated, noperators_rest);//R
  alps::numeric::matrix<M_TYPE> Green0_j_n(noperators_rest, num_rows_cols_updated);//Q

  for(int i=0;i<noperators_rest;++i) {
    for (int iv=0; iv<num_rows_cols_updated; ++iv) {
      Green0_n_j(iv,i) = mycast<M_TYPE>(green0_spline_for_M(flavor, noperators_rest+iv, i));
    }
  }
  for(int i=0;i<noperators_rest;++i){
    for (int iv=0; iv<num_rows_cols_updated; ++iv) {
      Green0_j_n(i,iv) = mycast<M_TYPE>(green0_spline_for_M(flavor, i, noperators_rest+iv));
    }
  }
  for (int iv2=0; iv2<num_rows_cols_updated; ++iv2) {
    for (int iv=0; iv<num_rows_cols_updated; ++iv) {
      Green0_n_n(iv, iv2) = mycast<M_TYPE>(green0_spline_for_M(flavor, noperators_rest+iv, noperators_rest+iv2));
    }
  }
  for (int iv=0; iv<num_rows_cols_updated; ++iv) {
    Green0_n_n(iv, iv) -= mycast<M_TYPE>(M[flavor].alpha_at(iv+noperators_rest));
  }

  if (compute_only_weight) {
    return shift_helper.det_rat = compute_det_ratio_replace_rows_cols(
        M[flavor].matrix(), Green0_j_n, Green0_n_j, Green0_n_n, shift_helper.Mmat[flavor], shift_helper.inv_tSp[flavor]);
  } else {
    return shift_helper.det_rat = compute_inverse_matrix_replace_rows_cols(
        M[flavor].matrix(), Green0_j_n, Green0_n_j, Green0_n_n, shift_helper.Mmat[flavor], shift_helper.inv_tSp[flavor]);
  }
}

template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::fastupdate_shift_finalize(const int flavor, const std::vector<int>& rows_cols_updated) {
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  const int num_rows_cols_updated = rows_cols_updated.size();
  const int noperators = num_rows(M[flavor].matrix());

  for (int i=num_rows_cols_updated-1; i>=0; --i) {
    M[flavor].swap_rows_cols(rows_cols_updated[num_rows_cols_updated-1-i], noperators-1-i);
  }
}

template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::fastupdate_spin_flip(const int flavor, const std::vector<int>& rows_cols_updated, const std::vector<M_TYPE>& old_alpha,
                                                  bool compute_only_weight) {
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  assert(old_alpha.size()==rows_cols_updated.size());

  const int num_rows_cols_updated = rows_cols_updated.size();

  std::vector<M_TYPE> diff(num_rows_cols_updated);

  for (int iv=0; iv<num_rows_cols_updated; ++iv) {
    diff[iv] = -M[flavor].alpha_at(rows_cols_updated[iv])+old_alpha[iv];
    //std::cout << " debug pos " << rows_cols_updated[iv] << " diff " << diff[iv] << std::endl;
  }

  if (compute_only_weight) {
    return compute_det_ratio_replace_diaognal_elements(M[flavor].matrix(), num_rows_cols_updated, rows_cols_updated, diff, true);
  } else {
    compute_det_ratio_replace_diaognal_elements(M[flavor].matrix(), num_rows_cols_updated, rows_cols_updated, diff, false);
  }
}
