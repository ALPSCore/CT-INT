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

    std::vector<size_t> pos_c_in_smallM(rank);//remember where we insert creation operators in M matrices.
    for (size_t i_rank=0; i_rank<rank; ++i_rank) {
      const size_t flavor_rank = new_vertex_type.flavors()[i_rank];
      pos_c_in_smallM[i_rank] = M[flavor_rank].creators().size();
      assert(new_vertex_type.sites()[2*i_rank]== new_vertex_type.sites()[2*i_rank+1]);//assume onsite U
      M[flavor_rank].creators().push_back(creator(flavor_rank, new_vertex_type.sites()[2*i_rank], time, n_matsubara));
      M[flavor_rank].annihilators().push_back(annihilator(flavor_rank, new_vertex_type.sites()[2*i_rank+1], time, n_matsubara));
      M[flavor_rank].alpha().push_back(new_vertex_type.get_alpha(af_state, i_rank));
      ++(helper.num_new_rows[flavor_rank]);
    }
    vertices_new.push_back(itime_vertex(v_type, af_state, time));
  }

  //fast update for each flavor
  double lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]==0) {
      continue;
    }
    lambda_prod *= fastupdate_up(flavor, true, helper.num_new_rows[flavor]); // true means compute_only_weight
  }
  //bit afraid of overlow. better to use log.
  const double rtmp = pow(-beta*Uijkl.n_vertex_type(),(double)n_vertices_add)/boost::math::binomial_coefficient<double>(vertices_new.size(),n_vertices_add);
  return rtmp*prod_Uval*lambda_prod;
}


void HubbardInteractionExpansion::perform_add(fastupdate_add_helper& helper, size_t n_vertices_add=1)
{
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]>0) {
      assert(helper.num_new_rows[flavor]==1);
      fastupdate_up(flavor,false,1);
    }
  }
  assert(num_rows(M[0].matrix())==vertices_new.size());
}


void HubbardInteractionExpansion::reject_add(fastupdate_add_helper& helper, size_t n_vertices_add=1)
{
  //get rid of the operators
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    for (int i=0; i<helper.num_new_rows[flavor]; ++i) {
      M[flavor].creators().pop_back();
      M[flavor].annihilators().pop_back();
      M[flavor].alpha().pop_back();
    }
  }
  //get rid of the vertex from vertex list
  //vertices.pop_back();
  vertices_new.pop_back();
  assert(num_rows(M[0].matrix())==vertices_new.size());
}


double HubbardInteractionExpansion::try_remove(unsigned int vertex_nr, fastupdate_remove_helper& helper)
{
  //get weight
  helper.clear();
  vertex_definition<GTYPE> vertex_def = Uijkl.get_vertex(vertices_new[vertex_nr].vertex_type());
  for (size_t i_rank=0; i_rank<vertex_def.rank(); ++i_rank) {
    const size_t flavor = vertex_def.flavors()[i_rank];
    //figure out and remember which rows (columns) are to be removed
    helper.rows_cols_removed[flavor].push_back(M[flavor].find_row_col(vertices_new[vertex_nr].time()));
  }
  GTYPE lambda_prod = 1.0;
  //sort lists of rows and cols to be removed in ascending order
  helper.sort_rows_cols();
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.rows_cols_removed[flavor].size()>0) {
      lambda_prod *= fastupdate_down(helper.rows_cols_removed[flavor], flavor, true);  // true means compute_only_weight
    }
  }
  double pert_order=num_rows(M[0].matrix());
  assert(num_rows(M[0].matrix())==vertices_new.size());
  assert(Uijkl.n_vertex_type()==n_site);
  //std::swap(vertices_new[vertex_nr], vertices_new[vertices_new.size()-1]);
  return  -pert_order/(beta*onsite_U*Uijkl.n_vertex_type())*lambda_prod;
}


void HubbardInteractionExpansion::perform_remove(unsigned int vertex_nr, fastupdate_remove_helper& helper)
{
  //perform fastupdate down
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.rows_cols_removed[flavor].size()>0) {
      //remove row and column
      fastupdate_down(helper.rows_cols_removed[flavor], flavor, false);  // false means really perform, not only compute weight
      //get rid of operators
      M[flavor].creators().pop_back();
      M[flavor].annihilators().pop_back();
      M[flavor].alpha().pop_back();
    }
  }
  //get rid of vertex list entries
  {
    std::vector<itime_vertex>::iterator it = vertices_new.begin();
    std::advance(it, vertex_nr);
    vertices_new.erase(it);//could be costly..
  }
}


void HubbardInteractionExpansion::reject_remove(fastupdate_remove_helper& helper)
{
  //do nothing
  return;
}
