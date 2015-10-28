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

double HubbardInteractionExpansion::try_add(fastupdate_add_helper& helper)
{
  //work array
  helper.clear();

  //select a new vertex
  double t = beta*random();
  const size_t v_type = static_cast<size_t>(random()*Uijkl.n_vertex_type());
  const vertex_definition<GTYPE> new_vertex = Uijkl.get_vertex(v_type);
  const size_t rank = new_vertex.rank();
  const size_t af_state = static_cast<size_t>(random()*new_vertex.num_af_states());
  assert(rank==2);

  const double Uval = new_vertex.Uval();

  std::vector<size_t> pos_c_in_smallM(rank);//remember where we insert creation operators in M matrices.
  //"i_rank" stands for a pair of creation and annihilation operators
  for (size_t i_rank=0; i_rank<rank; ++i_rank) {
    const size_t flavor_rank = new_vertex.flavors()[i_rank];
    pos_c_in_smallM[i_rank] = M[flavor_rank].creators().size();
    assert(new_vertex.sites()[2*i_rank]==new_vertex.sites()[2*i_rank+1]);
    M[flavor_rank].creators().push_back(creator(flavor_rank, new_vertex.sites()[2*i_rank], t, n_matsubara));
    M[flavor_rank].annihilators().push_back(annihilator(flavor_rank, new_vertex.sites()[2*i_rank+1], t, n_matsubara));
    M[flavor_rank].alpha().push_back(new_vertex.get_alpha(af_state, i_rank));
    ++(helper.num_new_rows[flavor_rank]);
  }

  //keep track of vertex list
  vertices_new.push_back(itime_vertex(v_type, af_state, pos_c_in_smallM));
  {
    spin_t flavor0 = 0;
    spin_t flavor1 = 1;
    size_t site = new_vertex.sites()[0];
    vertices.push_back(vertex(flavor0, site, M[0].creators().size()-1, M[0].annihilators().size()-1,
                              flavor1, site, M[1].creators().size()-1, M[1].annihilators().size()-1, Uval));
  }

  //fast update for each flavor
  double lambda_prod = 1.0;
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    if (helper.num_new_rows[flavor]==0) {
      continue;
    }
    assert(helper.num_new_rows[flavor]==1);
    lambda_prod *= fastupdate_up(flavor, true); // true means compute_only_weight
  }
  //std::cout << "debug " <<  -beta*Uval*Uijkl.n_vertex_type()/(vertices_new.size())*lambda_prod << std::endl;
  return -beta*Uval*Uijkl.n_vertex_type()/(vertices_new.size())*lambda_prod;
}


void HubbardInteractionExpansion::perform_add(fastupdate_add_helper& helper)
{
  for (spin_t flavor=0; flavor<n_flavors; ++flavor) {
    //std::cout << " debug_add " << flavor << " " << helper.num_new_rows[flavor] << std::endl;
    if (helper.num_new_rows[flavor]>0) {
      assert(helper.num_new_rows[flavor]==1);
      fastupdate_up(flavor,false);
    }
  }
}


void HubbardInteractionExpansion::reject_add(fastupdate_add_helper& helper)
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
  vertices.pop_back();
  vertices_new.pop_back();
}


double HubbardInteractionExpansion::try_remove(unsigned int vertex_nr, fastupdate_remove_helper& helper)
{
  //get weight
  helper.clear();
  vertex_definition<GTYPE> vertex_def = Uijkl.get_vertex(vertices_new[vertex_nr].vertex_type());
  for (size_t i_rank=0; i_rank<vertex_def.rank(); ++i_rank) {
    const size_t flavor = vertex_def.flavors()[i_rank];
    ++helper.num_removed_rows[flavor];
  }
  GTYPE lambda_prod = 1.0;
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    assert(helper.num_removed_rows[flavor]<=1);
    if (helper.num_removed_rows[flavor]==1) {
      lambda_prod *= fastupdate_down(vertex_nr, flavor, true);  // true means compute_only_weight
    }
  }
  double pert_order=num_rows(M[0].matrix());
  assert(num_rows(M[0].matrix())==vertices_new.size());
  assert(Uijkl.n_vertex_type()==n_site);
  //swap vertices
  std::swap(vertices[vertex_nr], vertices[vertices.size()-1]);
  std::swap(vertices_new[vertex_nr], vertices_new[vertices_new.size()-1]);
  return  -pert_order/(beta*onsite_U*Uijkl.n_vertex_type())*lambda_prod;
}


void HubbardInteractionExpansion::perform_remove(unsigned int vertex_nr, fastupdate_remove_helper& helper)
{
  //perform fastupdate down
  for (size_t flavor=0; flavor<n_flavors; ++flavor) {
    assert(helper.num_removed_rows[flavor]<=1);
    if (helper.num_removed_rows[flavor]>0) {
      assert(helper.num_removed_rows[flavor]==1);
      //remove row and column
      fastupdate_down(vertex_nr, flavor, false);  // false means really perform, not only compute weight
      //get rid of operators
      M[flavor].creators().pop_back();
      M[flavor].annihilators().pop_back();
      M[flavor].alpha().pop_back();
    }
  }
  //get rid of vertex list entries
  vertices.pop_back();
  vertices_new.pop_back();
}


void HubbardInteractionExpansion::reject_remove(fastupdate_remove_helper& helper)
{
  //do nothing
  return;
}
