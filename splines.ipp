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

///as the name says: linearly interpolate between two points. this is VERY
//EXPENSIVE because of the division.

template<class X, class Y> inline Y linear_interpolate(const X x0, const X x1, const Y y0, const Y y1, const X x)
{
  return y0 + (x-x0)/(x1-x0)*(y1-y0);
}

//with correct treatment of equal-time Green's function
//take care of the ambiguity of the equal-time Green's function.
//We appreciate the ordering of the creation and annihilation operators in the vertex.
template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::green0_spline_for_M(const annihilator& c, const creator& cdagger)
{
  assert(c.flavor()==cdagger.flavor());
  const int flavor = c.flavor();
  typename TYPES::M_TYPE result;
  if (c.t().time()!=cdagger.t().time()) {
    itime_t delta_t=c.t().time()-cdagger.t().time();
    site_t site1 = c.s();
    site_t site2 = cdagger.s();
    result = green0_spline_new(delta_t, flavor, site1, site2);
  } else {
    itime_t delta_t=c.t().time()-cdagger.t().time();
    site_t site1 = c.s();
    site_t site2 = cdagger.s();
    itime_t time_shift = c.t().small_index() > cdagger.t().small_index() ? beta*1E-10 : -beta*1E-10;
    result = green0_spline_new(delta_t+time_shift, flavor, site1, site2);
  }
  if(std::abs(result-g0_intpl(c,cdagger))>1E-8) {
    std::cout << " time " << c.t().time() - cdagger.t().time() << std::endl;
    std::cout << " small_inex " << c.t().small_index() << " " <<  cdagger.t().small_index() << std::endl;
    std::cout << " debug " << result << " " << g0_intpl(c, cdagger) << std::endl;
  }
  assert(std::abs(result-g0_intpl(c,cdagger))<1E-8);
  return result;
}

///Compute the bare green's function for a given flavor, site, and imaginary time.
template<class TYPES>
typename TYPES::M_TYPE
InteractionExpansion<TYPES>::green0_spline_new(const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2)
{
  assert(delta_t<= beta);
  assert(delta_t>=-beta);
  if(delta_t*delta_t < almost_zero && delta_t>0){
    return mycast<M_TYPE>(bare_green_itime(0, site1, site2, flavor));
  }
  else if(delta_t*delta_t < almost_zero && delta_t<0){
    return mycast<M_TYPE>(-bare_green_itime(n_tau, site1, site2, flavor));
  }
  else if(delta_t>0){
    int time_index_1 = (int)(delta_t*n_tau*temperature);
    int time_index_2 = time_index_1+1;
    return linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,
                              mycast<M_TYPE>(bare_green_itime(time_index_1, site1, site2, flavor)),
                              mycast<M_TYPE>(bare_green_itime(time_index_2, site1, site2, flavor)),
                              delta_t);
  }
  else{
    int time_index_1 = (int)(delta_t*n_tau*temperature+n_tau);
    int time_index_2 = time_index_1+1;
    return -linear_interpolate((double)time_index_1*beta*n_tau_inv, (double)time_index_2*beta*n_tau_inv,
                               mycast<M_TYPE>(bare_green_itime(time_index_1,site1,site2,flavor)),
                               mycast<M_TYPE>(bare_green_itime(time_index_2,site1,site2,flavor)),
                               delta_t+beta);
  }
}
