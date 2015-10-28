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

double InteractionExpansion::fastupdate_up(const int flavor, bool compute_only_weight, size_t n_vertices_add=1)
{
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  unsigned int noperators = num_rows(M[flavor].matrix());
  //current size of M: number of vertices - n_vertices_add. We need to add the last vertex.
  //A pointer to the creator and annihilator is already stored in creators_ and annihilators_
  //at positions M[flavor].size() ... M[flavor].size()+n_vertices_add-1;
  alps::numeric::matrix<GTYPE> Green0_n_n(n_vertices_add, n_vertices_add);
  alps::numeric::matrix<GTYPE> Green0_n_j(n_vertices_add, noperators);
  alps::numeric::matrix<GTYPE> Green0_j_n(noperators, n_vertices_add);
  for(unsigned int i=0;i<noperators;++i) {
    for (size_t iv=0; iv<n_vertices_add; ++iv) {
      Green0_n_j(iv,i) = green0_spline_new(M[flavor].annihilators()[noperators+iv], M[flavor].creators()[i]);
    }
  }
  for(unsigned int i=0;i<noperators;++i){
    for (size_t iv=0; iv<n_vertices_add; ++iv) {
      Green0_j_n(i,iv) = green0_spline_new(M[flavor].annihilators()[i], M[flavor].creators()[noperators+iv]);
    }
  }
  for (size_t iv2=0; iv2<n_vertices_add; ++iv2) {
    for (size_t iv=0; iv<n_vertices_add; ++iv) {
      Green0_n_n(iv, iv2) = green0_spline_new(M[flavor].annihilators()[noperators+iv], M[flavor].creators()[noperators+iv2]);
    }
  }

  //B: Green0_j_n
  //C: Green0_n_j
  //D: Green0_n_n
  //invA: M[flavor].matrix()
  if(compute_only_weight){
    return compute_det_ratio_up(Green0_j_n, Green0_n_j, Green0_n_n, M[flavor].matrix());
  } else {
    return compute_inverse_matrix_up2(Green0_j_n, Green0_n_j, Green0_n_n, M[flavor].matrix(), M[flavor].matrix());
  }
}



///Fastupdate formulas, remove order by one (remove a vertex). If necessary
///also take track of the Green's function.
double InteractionExpansion::fastupdate_down(const int operator_nr, const int flavor, bool compute_only_weight)
{
  using std::swap;
  using alps::numeric::column_view;
  using alps::numeric::vector;
  using alps::numeric::matrix;
  assert(num_rows(M[flavor].matrix()) == num_cols(M[flavor].matrix()));
  //perform updates according to formula 21.1, 21.2
  unsigned int noperators=num_rows(M[flavor].matrix());  //how many operators do we have in total?
  if(compute_only_weight){
    return M[flavor].matrix()(operator_nr,operator_nr);
  }
  //swap rows and colums of M <-> move selected vertex to the end.
  for(unsigned int i=0;i<noperators;++i){
    swap(M[flavor].matrix()(i,noperators-1), M[flavor].matrix()(i,operator_nr));
  }
  for(unsigned int i=0;i<noperators;++i){
    swap(M[flavor].matrix()(noperators-1,i),M[flavor].matrix()(operator_nr,i));
  }
  //swap creator and annihilator
  swap(M[flavor].creators()[operator_nr],     M[flavor].creators()[noperators-1]);
  swap(M[flavor].annihilators()[operator_nr], M[flavor].annihilators()[noperators-1]);
  swap(M[flavor].alpha()[operator_nr],        M[flavor].alpha()[noperators-1]);
  double Mnn=M[flavor].matrix()(noperators-1,noperators-1);
  //now perform fastupdate of M
  vector<double> lastrow(noperators-1);
  vector<double> lastcolumn(column_view<matrix<double> >(M[flavor].matrix(),noperators-1));
  for(unsigned int j=0;j<noperators-1;++j){
    lastrow[j]=M[flavor].matrix()(noperators-1, j);
  }
  if(noperators>1)
    lastcolumn *= -1./Mnn;
  resize(M[flavor].matrix(),noperators-1,noperators-1);  //lose the last row and last column, reduce size by one, but keep contents.
  if(noperators>1)
    boost::numeric::bindings::blas::ger(1.0,lastcolumn,lastrow,M[flavor].matrix());
  return Mnn;  //the determinant ratio det D_{k-1}/det D_{k}
}




