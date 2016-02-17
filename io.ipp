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
void InteractionExpansion<TYPES>::print(std::ostream &os){
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"*** ALPS InteractionExpansion solver for multi-orbital cluster impurity model                           ***"<<std::endl;
  os<<"*** Hiroshi Shinaoka                                                                                    ***"<<std::endl;
  os<<"*** This code is based on a previous version written by                                                 ***"<<std::endl;
  os<<"*** Emanuel Gull, Philipp Werner, Sebastian Fuchs, Brigitte Surer, Thomas Pruschke, and Matthias Troyer ***"<<std::endl;
  os<<"*** Published under the ALPS cite-me license in:                                                        ***"<<std::endl;
  os<<"*** ***************** Computer Physics Communications 182, 1078 (2011). ******************************* ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"***                                                                                                     ***"<<std::endl;
  os<<"*** implementing the interaction expansion algorithm by Rubtsov et al., JETP Letters 80, 61.            ***"<<std::endl;
  os<<"***                                                                                                     ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"max order\t"<<max_order<<"\tn_flavors: "
    <<n_flavors<<"\tn_site: "<<n_site
    <<"\tn_matsubara: "<<n_matsubara<<std::endl;
  os<<"n_tau: "<<n_tau<<"\tmc steps: "<<mc_steps
    <<"\ttherm steps: "<<therm_steps<<std::endl;
  
  os<<"beta: "<<beta<<std::endl;

  os<<"recalc period: "<<recalc_period<<"\tmeasurement period: "<< measurement_period
    <<"\tconvergence period: "<< convergence_check_period<<std::endl;
  os<<"almost zero: "<<almost_zero<<std::endl;
}
