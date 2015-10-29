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
#include <boost/math/special_functions/binomial.hpp>

double HubbardInteractionExpansion::try_add(fastupdate_add_helper& helper, size_t n_vertices_add=1)
{
  //work array
  helper.clear();

  //add vertices one by one
  double prod_Uval = 1.0;
  for (size_t iv=0; iv<n_vertices_add; ++iv) {
    const double time = beta*random();
    const size_t v_type = static_cast<size_t>(random()*Uijkl.n_vertex_type());
    const vertex_definition<GTYPE> new_vertex_type = Uijkl.get_vertex(v_type);
    const size_t rank = new_vertex_type.rank();
    const size_t af_state = static_cast<size_t>(random()* new_vertex_type.num_af_states());
    prod_Uval *= new_vertex_type.Uval();
    assert(rank==2);

    for (size_t i_rank=0; i_rank<rank; ++i_rank) {
      const size_t flavor_rank = new_vertex_type.flavors()[i_rank];
      assert(new_vertex_type.sites()[2*i_rank]== new_vertex_type.sites()[2*i_rank+1]);//assume onsite U

      M[flavor_rank].creators().push_back(creator(flavor_rank, new_vertex_type.sites()[2*i_rank], time, n_matsubara));
      M[flavor_rank].annihilators().push_back(annihilator(flavor_rank, new_vertex_type.sites()[2*i_rank+1], time, n_matsubara));
      M[flavor_rank].alpha().push_back(new_vertex_type.get_alpha(af_state, i_rank));
      M[flavor_rank].vertex_info().push_back(std::pair<vertex_t,size_t>(v_type,i_rank));

      ++(helper.num_new_rows[flavor_rank]);
    }
    itime_vertices.push_back(itime_vertex(v_type, af_state, time, rank));
  }

  //fast update for each flavor
  double lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]>0) {
      lambda_prod *= fastupdate_up(flavor, true, helper.num_new_rows[flavor]); // true means compute_only_weight
    }
  }
  const double rtmp = pow(-beta*Uijkl.n_vertex_type(),(double)n_vertices_add)
                      /boost::math::binomial_coefficient<double>(itime_vertices.size(),n_vertices_add);
  return rtmp*prod_Uval*lambda_prod;
}


void HubbardInteractionExpansion::perform_add(fastupdate_add_helper& helper, size_t n_vertices_add=1)
{
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]>0) {
      fastupdate_up(flavor,false,helper.num_new_rows[flavor]);
    }
  }
  M.sanity_check(itime_vertices);
  //ssert(num_rows(M[0].matrix())== itime_vertices.size());
}


void HubbardInteractionExpansion::reject_add(fastupdate_add_helper& helper, size_t n_vertices_add=1)
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


double HubbardInteractionExpansion::try_remove(const std::vector<size_t>& vertices_nr, fastupdate_remove_helper& helper)
{
  //get weight
  //figure out and remember which rows (columns) are to be removed
  // lists of rows and cols will be sorted in ascending order
  M.sanity_check(itime_vertices);
  helper.clear();
  double prod_Uval = 1.0;
  for (size_t iv=0; iv<vertices_nr.size(); ++iv) {
    vertex_definition<GTYPE> vertex_def = Uijkl.get_vertex(itime_vertices[vertices_nr[iv]].type());
    prod_Uval *= vertex_def.Uval();
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
  double r1 = boost::math::binomial_coefficient<double>(itime_vertices.size(),vertices_nr.size());
  double r2 = pow(-beta*Uijkl.n_vertex_type()*prod_Uval,(double)vertices_nr.size());
  return (r1/r2)*lambda_prod;
}


void HubbardInteractionExpansion::perform_remove(const std::vector<size_t>& vertices_nr, fastupdate_remove_helper& helper)
{
  //perform fastupdate down
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.rows_cols_removed[flavor].size()>0) {
      //remove rows and columns
      fastupdate_down(helper.rows_cols_removed[flavor], flavor, false);  // false means really perform, not only compute weight
      //get rid of operators
      for (int iop=0; iop<helper.rows_cols_removed[flavor].size(); ++iop) {
        M[flavor].pop_back_op();
      }
    }
  }
  //get rid of vertex list entries. I know this is crapy.
  remove_elements_from_vector(itime_vertices, vertices_nr);
  M.sanity_check(itime_vertices);
}


void HubbardInteractionExpansion::reject_remove(fastupdate_remove_helper& helper)
{
  //do nothing
  return;
}
