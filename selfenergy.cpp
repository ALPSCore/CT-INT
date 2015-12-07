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
#include "operator.hpp"
#include <complex>
#include <alps/alea.h>
#include <alps/alea/simpleobseval.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>


typedef alps::SignedObservable<alps::RealVectorObservable> signed_vec_obs_t;
typedef alps::RealVectorObservable vec_obs_t;
typedef alps::SimpleObservable<double,alps::DetailedBinning<double> > simple_obs_t;
typedef const alps::SimpleObservable<double,alps::DetailedBinning<double> > const_simple_obs_t;


#ifdef SSE
#include<emmintrin.h>
class twocomplex{
public:
  inline twocomplex(){};
  inline twocomplex(const std::complex<double>&p, const std::complex<double> &q){
    r=_mm_loadl_pd(r, &(p.real())); //load two complex numbers.
    r=_mm_loadh_pd(r, &(q.real()));
    i=_mm_loadl_pd(i, &(p.imag()));
    i=_mm_loadh_pd(i, &(q.imag()));
  }
  inline void store(std::complex<double> &p, std::complex<double> &q){
    _mm_store_sd(&(p.real()), r);
    _mm_store_sd(&(p.imag()), i);
    _mm_storeh_pd(&(q.real()), r);
    _mm_storeh_pd(&(q.imag()), i);
  }
  __m128d r;
  __m128d i;
};

inline twocomplex fastcmult(const twocomplex &a, const twocomplex &b)
{
  twocomplex c;
  c.r = _mm_sub_pd(_mm_mul_pd(a.r,b.r), _mm_mul_pd(a.i,b.i));
  c.i = _mm_add_pd(_mm_mul_pd(a.r,b.i), _mm_mul_pd(a.i,b.r));
  return c;
}
#endif



void InteractionExpansion::compute_W_matsubara()
{
  Wk_t Wk(boost::extents[n_flavors][n_site][n_site][n_matsubara]);
  measure_Wk(Wk, n_matsubara_measurements);
  measure_densities();
}

void InteractionExpansion::measure_Wk(Wk_t& Wk, const unsigned int nfreq)
{
  //clear contents of Wk
  std::fill(Wk.origin(), Wk.origin()+Wk.num_elements(), 0);

  //allocate memory for work
  size_t max_mat_size = 0;
  for (unsigned int z=0; z<n_flavors; ++z) {
    assert( num_rows(M[z].matrix()) == num_cols(M[z].matrix()) );
    max_mat_size = std::max(max_mat_size, num_rows(M[z].matrix()));
  }
  alps::numeric::matrix<std::complex<double> > GR(max_mat_size, n_site),
          GL(n_site, max_mat_size), MGR(max_mat_size, n_site), GLMGR(max_mat_size, max_mat_size);

  //Note: vectorized for the loops over operators (should be optimal for large beta)
  for (size_t z=0; z<n_flavors; ++z) {
    const size_t Nv = num_rows(M[z].matrix());

    if (Nv==0) {
      continue;
    }

    //We do not need any more?
    //for(size_t p=0;p<Nv;++p) {
      //M[z].creators()[p].compute_exp(n_matsubara, -1);
    //}
    //for(size_t q=0;q<Nv;++q) {
      //M[z].annihilators()[q].compute_exp(n_matsubara, +1);
    //}

    for(unsigned int i_freq=0; i_freq <nfreq; ++i_freq) {
      GR.resize(Nv, n_site);
      MGR.resize(Nv, n_site);
      GL.resize(n_site, Nv);
      GLMGR.resize(n_site, n_site);

      //GR
      for(unsigned int p=0;p<Nv;++p) {
        const size_t site_p = M[z].annihilators()[p].s();
        const double phase = M[z].annihilators()[p].t().time()*(2*i_freq+1)*M_PI/beta;
        const std::complex<double> exp = std::complex<double>(std::cos(phase), -std::sin(phase));
        for (size_t site=0; site<n_site; ++site) {
          GR(p, site) = bare_green_matsubara(i_freq, site_p, site, z)*exp;
        }
      }

      //GL
      for(unsigned int q=0;q<Nv;++q) {
        const size_t site_q = M[z].creators()[q].s();
        const double phase = M[z].creators()[q].t().time()*(2*i_freq+1)*M_PI/beta;
        const std::complex<double> exp = std::complex<double>(std::cos(phase), std::sin(phase));
        //const std::complex<double> exp = M[z].annihilators()[q].exp_iomegat()[i_freq];
        for (size_t site=0; site<n_site; ++site) {
          GL(site, q) = bare_green_matsubara(i_freq, site, site_q, z)*exp;
        }
      }

      //clear MGR, GLMGR
      std::fill(MGR.get_values().begin(), MGR.get_values().end(), 0);
      std::fill(GLMGR.get_values().begin(), GLMGR.get_values().end(), 0);

      gemm(M[z].matrix(), GR, MGR);
      gemm(GL, MGR, GLMGR);

      for (unsigned int site1=0; site1<n_site; ++site1) {
        for (unsigned int site2=0; site2<n_site; ++site2) {
          Wk[z][site1][site2][i_freq] += GLMGR(site1, site2);
        }
      }

    }//i_freq
  }//z

  int num_nd = 0;
  for (int iv=0; iv<itime_vertices.size(); ++iv) {
    if (!itime_vertices[iv].is_density_type())
      ++num_nd;
  }

  if (params.defined("OUTPUT_Wk") ? params["OUTPUT_Wk"] : false) {
    std::valarray<double> pert_vertex(Uijkl.n_vertex_type());
    for (std::vector<itime_vertex>::const_iterator it=itime_vertices.begin(); it!=itime_vertices.end(); ++it) {
      assert(it->type()>=0 && it->type()<Uijkl.n_vertex_type());
      ++pert_vertex[it->type()];
    }
    std::cout << " Wk_node= " << node << " "
      << Wk[0][0][0][0].real() << " " << Wk[0][0][0][0].imag() << " "
      << Wk[0][1][1][0].real() << " " << Wk[0][1][1][0].imag() << " "
      << Wk[0][2][2][0].real() << " " << Wk[0][2][2][0].imag() << " "
      << Wk[0][3][3][0].real() << " " << Wk[0][3][3][0].imag() << " "
      << Wk[1][0][0][0].real() << " " << Wk[1][0][0][0].imag() << " "
      << Wk[1][1][1][0].real() << " " << Wk[1][1][1][0].imag() << " "
      << Wk[1][2][2][0].real() << " " << Wk[1][2][2][0].imag() << " "
      << Wk[1][3][3][0].real() << " " << Wk[1][3][3][0].imag() << " "
      << " step " << step;
    std::cout << " pert_PH " << pert_vertex[6] << " " << pert_vertex[4] << " pert_SF " << pert_vertex[8] << " " << pert_vertex[1] << std::endl;
  }

  //if (!is_quantum_number_conserved(itime_vertices)) {
    //std::cout << "w/o check range " << is_quantum_number_conserved(itime_vertices,false) << std::endl;
    //std::vector<itime_vertex> itime_vertices_sorted(itime_vertices);
    //std::sort(itime_vertices_sorted.begin(), itime_vertices_sorted.end());
    //print_vertices(std::cout, itime_vertices_sorted);
  //}

  for(unsigned int flavor=0;flavor<n_flavors;++flavor) {
    for (unsigned int site1 = 0; site1 < n_site; ++site1) {
      for (unsigned int site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Wk_real_name, Wk_imag_name;
        Wk_real_name << "Wk_real_" << flavor << "_" << site1 << "_" << site2;
        Wk_imag_name << "Wk_imag_" << flavor << "_" << site1 << "_" << site2;
        std::valarray<double> Wk_real(nfreq);
        std::valarray<double> Wk_imag(nfreq);
        for (unsigned int w = 0; w < nfreq; ++w) {
          Wk_real[w] = Wk[flavor][site1][site2][w].real();
          Wk_imag[w] = Wk[flavor][site1][site2][w].imag();
        }
        measurements[Wk_real_name.str().c_str()] << static_cast<std::valarray<double> > (Wk_real * sign);
        measurements[Wk_imag_name.str().c_str()] << static_cast<std::valarray<double> > (Wk_imag * sign);
      }
    }
  }
}

void InteractionExpansion::compute_Sl() {
  static boost::multi_array<std::complex<double>,3> Sl(boost::extents[n_site][n_site][n_legendre]);
  const size_t num_random_walk = 100;

  //Work arrays
  size_t max_mat_size = 0;
  for (unsigned int z=0; z<n_flavors; ++z) {
    assert( num_rows(M[z].matrix()) == num_cols(M[z].matrix()) );
    max_mat_size = std::max(max_mat_size, num_rows(M[z].matrix()));
  }
  std::vector<double> legendre_vals(n_legendre), sqrt_vals(n_legendre);
  for(unsigned int i_legendre=0; i_legendre<n_legendre; ++i_legendre) {
    sqrt_vals[i_legendre] = std::sqrt(2.0*i_legendre+1.0);
  }
  alps::numeric::matrix<GTYPE> gR(max_mat_size, n_site), M_gR(max_mat_size, n_site);

  for (unsigned int z=0; z<n_flavors; ++z) {
    const size_t Nv = num_rows(M[z].matrix());
    gR.resize(Nv, n_site);
    M_gR.resize(Nv, n_site);
    std::fill(Sl.origin(),Sl.origin()+Sl.num_elements(),0.0);//clear the content for safety

    //shift times of operators by time_shift
    for (std::size_t random_walk=0; random_walk<num_random_walk; ++random_walk) {
      const double time_shift = beta * random();

      for (unsigned int p = 0; p < Nv; ++p) {//annihilation operators
        const double tmp = M[z].annihilators()[p].t().time() + time_shift;
        const double time_a_shifted = tmp < beta ? tmp : tmp - beta;
        const double coeff = tmp < beta ? 1 : -1;

        //interpolate G0
        for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
          gR(p, site_B) = coeff*green0_spline_new(time_a_shifted, z, M[z].annihilators()[p].s(), site_B);
        }
      }

      gemm(M[z].matrix(), gR, M_gR);

      for (unsigned int q = 0; q < Nv; ++q) {//creation operators
        const unsigned int site_c = M[z].creators()[q].s();
        const double tmp = M[z].creators()[q].t().time() + time_shift;
        const double time_c_shifted = tmp < beta ? tmp : tmp - beta;
        const double coeff = tmp < beta ? 1 : -1;

        legendre_transformer.compute_legendre(2 * time_c_shifted/beta - 1.0, legendre_vals);//P_l[x(tau_q)]
        for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
          for (unsigned int i_legendre = 0; i_legendre < n_legendre; ++i_legendre) {
            Sl[site_c][site_B][i_legendre] -= coeff*sqrt_vals[i_legendre] * legendre_vals[i_legendre] * M_gR(q, site_B);
          }
        }
      }
    }//random_walk

    //pass data to ALPS library
    for (unsigned int site1 = 0; site1 < n_site; ++site1) {
      for (unsigned int site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Sl_real_name, Sl_imag_name;
        Sl_real_name << "Sl_real_" << z << "_" << site1 << "_" << site2;
        Sl_imag_name << "Sl_imag_" << z << "_" << site1 << "_" << site2;
        std::valarray<double> Sl_real(n_legendre);
        std::valarray<double> Sl_imag(n_legendre);
        for (unsigned int i_legendre = 0; i_legendre < n_legendre; ++i_legendre) {
          Sl_real[i_legendre] = Sl[site1][site2][i_legendre].real()/num_random_walk;
          Sl_imag[i_legendre] = Sl[site1][site2][i_legendre].imag()/num_random_walk;
        }
        measurements[Sl_real_name.str().c_str()] << static_cast<std::valarray<double> > (Sl_real * sign);
        measurements[Sl_imag_name.str().c_str()] << static_cast<std::valarray<double> > (Sl_imag * sign);
      }//site2
    }//site1

  }//z
}


void InteractionExpansion::measure_densities()
{
  std::vector< std::vector<double> > dens(n_flavors);
  for(unsigned int z=0;z<n_flavors;++z){
    dens[z].resize(n_site);
    memset(&(dens[z][0]), 0., sizeof(double)*(n_site));
  }
  double tau = beta*random();
  for (unsigned int z=0; z<n_flavors; ++z) {
    const size_t Nv = num_rows(M[z].matrix());
    alps::numeric::vector<double> g0_tauj(Nv);
    alps::numeric::vector<double> M_g0_tauj(Nv);
    alps::numeric::vector<double> g0_taui(Nv);
    for (unsigned int s=0;s<n_site;++s) {             
      for (unsigned int j=0;j<Nv;++j)
        g0_tauj[j] = green0_spline_new(M[z].annihilators()[j].t().time()-tau, z, M[z].annihilators()[j].s(), s);//CHECK THE TREATMENT OF EQUAL-TIME Green's function
      for (unsigned int i=0;i<Nv;++i)
        g0_taui[i] = green0_spline_new(tau-M[z].creators()[i].t().time(),z, s, M[z].creators()[i].s());
      if (num_rows(M[z].matrix())>0)
          gemv(M[z].matrix(),g0_tauj,M_g0_tauj);
      dens[z][s] += green0_spline_new(-beta*1E-10,z,s,s);//tau=-0
      for (unsigned int j=0;j<Nv;++j)
        dens[z][s] -= g0_taui[j]*M_g0_tauj[j]; 
    }
  }
  std::valarray<double> densities(0., n_flavors);
  for (unsigned int z=0; z<n_flavors; ++z) {                  
    std::valarray<double> densmeas(n_site);
    for (unsigned int i=0; i<n_site; ++i) {
      densities[z] += dens[z][i];
      densmeas[i] = dens[z][i];
    }
    measurements["densities_"+boost::lexical_cast<std::string>(z)] << static_cast<std::valarray<double> > (densmeas*sign);
    densities[z] /= n_site;
    densities[z] = densities[z];
  }
  measurements["densities"] << static_cast<std::valarray<double> > (densities*sign);
  double density_correlation = 0.;
  for (unsigned int i=0; i<n_site; ++i) {
    density_correlation += (dens[0][i])*(dens[1][i]);
  }
  density_correlation /= n_site;
  measurements["density_correlation"] << (density_correlation*sign);
  std::valarray<double> ninj(n_site*n_site*4);
  for (unsigned int i=0; i<n_site; ++i) {
    for (unsigned int j=0; j<n_site; ++j) {
      ninj[i*n_site+j] = (dens[0][i])*(dens[0][j]);
      ninj[i*n_site+j+1] = (dens[0][i])*(dens[1][j]);
      ninj[i*n_site+j+2] = (dens[1][i])*(dens[0][j]);
      ninj[i*n_site+j+3] = (dens[1][i])*(dens[1][j]);
    }
  }
  measurements["n_i n_j"] << static_cast<std::valarray<double> > (ninj*sign);
}

void evaluate_selfenergy_measurement_matsubara(const alps::results_type<HubbardInteractionExpansion>::type &results,
                                                                        matsubara_green_function_t &green_matsubara_measured,
                                                                        const matsubara_green_function_t &bare_green_matsubara,
                                                                        std::vector<double>& densities,
                                                                        const double &beta, std::size_t n_site,
                                                                        std::size_t n_flavors, std::size_t n_matsubara, std::size_t n_matsubara_meas) {
  green_matsubara_measured.clear();
  std::cout << "evaluating self energy measurement: matsubara, reciprocal space" << std::endl;
  double sign = results["Sign"].mean<double>();
  for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
    for (std::size_t site1 = 0; site1 < n_site; ++site1) {
      for (std::size_t site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Wk_real_name, Wk_imag_name;
        Wk_real_name << "Wk_real_" << flavor1 << "_" << site1 << "_" << site2;
        Wk_imag_name << "Wk_imag_" << flavor1 << "_" << site1 << "_" << site2;
        std::vector<double> mean_real = results[Wk_real_name.str().c_str()].mean<std::vector<double> >();
        std::vector<double> mean_imag = results[Wk_imag_name.str().c_str()].mean<std::vector<double> >();
        assert(mean_real.size()==n_matsubara_meas && mean_imag.size()==n_matsubara_meas);
        if (!(mean_real.size()==n_matsubara_meas && mean_imag.size()==n_matsubara_meas)) {
          throw std::logic_error("Logic error in evaluate_selfenergy_measurement_matsubara.");
        }
        for (unsigned int w = 0; w < n_matsubara_meas; ++w) {
          green_matsubara_measured(w, site1, site2, flavor1) = -std::complex<double>(mean_real[w], mean_imag[w]) / (sign*beta);
        }
      }
    }
  }
  for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
    for (std::size_t site1 = 0; site1 < n_site; ++site1) {
      for (std::size_t site2 = 0; site2 < n_site; ++site2) {
        for (std::size_t w = 0; w < n_matsubara; ++w) {
          green_matsubara_measured(w, site1, site2, flavor1) += bare_green_matsubara(w, site1, site2, flavor1);
        }
      }
    }
  }
  std::vector<double> dens = results["densities"].mean<std::vector<double> >();
  for (std::size_t z=0; z<n_flavors; ++z)
    densities[z] = dens[z];

}

void evaluate_selfenergy_measurement_legendre(const alps::results_type<HubbardInteractionExpansion>::type &results,
                                               matsubara_green_function_t &green_matsubara_measured,
                                               matsubara_green_function_t &sigma_green_matsubara_measured,
                                               itime_green_function_t &green_itime_measured,
                                               const matsubara_green_function_t &bare_green_matsubara,
                                               std::vector<double>& densities,
                                               const double &beta, std::size_t n_site,
                                               std::size_t n_flavors, std::size_t n_matsubara, std::size_t n_tau, std::size_t n_legendre) {
  green_matsubara_measured.clear();
  green_itime_measured.clear();
  std::cout << "evaluating self energy measurement: lengendre, real space" << std::endl;
  double sign = results["Sign"].mean<double>();

  //Legendre expansion utils
  LegendreTransformer legendre_transformer(n_matsubara,n_legendre);
  const alps::numeric::matrix<std::complex<double> > & Tnl(legendre_transformer.Tnl());

  boost::multi_array<std::complex<double>,4> Sl(boost::extents[n_legendre][n_site][n_site][n_flavors]);
  boost::multi_array<std::complex<double>,4> Sw(boost::extents[n_matsubara][n_site][n_site][n_flavors]);

  //load S_l
  for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
    for (std::size_t site1 = 0; site1 < n_site; ++site1) {
      for (std::size_t site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Sl_real_name, Sl_imag_name;
        Sl_real_name << "Sl_real_" << flavor1 << "_" << site1 << "_" << site2;
        Sl_imag_name << "Sl_imag_" << flavor1 << "_" << site1 << "_" << site2;
        std::vector<double> mean_real = results[Sl_real_name.str().c_str()].mean<std::vector<double> >();
        std::vector<double> mean_imag = results[Sl_imag_name.str().c_str()].mean<std::vector<double> >();
        for (unsigned int i_l = 0; i_l < n_legendre; ++i_l) {
          Sl[i_l][site1][site2][flavor1] = std::complex<double>(mean_real[i_l], mean_imag[i_l])/sign;
        }
      }
    }
  }

  //compute S(iomega_n)
  {
//    std::cout << "n_matsubara " << n_matsubara << std::endl;
    alps::numeric::matrix<std::complex<double> > Sl_vec(n_legendre,1), Sw_vec(n_matsubara,1);//tmp arrays for gemm
    for (std::size_t flavor1 = 0; flavor1 < n_flavors; ++flavor1) {
      for (std::size_t site1 = 0; site1 < n_site; ++site1) {
        for (std::size_t site2 = 0; site2 < n_site; ++site2) {
          for(std::size_t i_l = 0; i_l < n_legendre;  ++ i_l) {
            Sl_vec(i_l,0) = Sl[i_l][site1][site2][flavor1];
          }
          std::fill(Sw_vec.get_values().begin(), Sw_vec.get_values().end(), 0.0);
          alps::numeric::gemm(Tnl,Sl_vec,Sw_vec);
          for (std::size_t w = 0; w < n_matsubara; ++w) {
            Sw[w][site1][site2][flavor1] = Sw_vec(w,0);
          }
        }
      }
    }
  }
//  for(std::size_t w = 0; w < n_matsubara;  ++w) {
//    std::cout << " w " << w << Sw[w][0][0][0] << std::endl;
//  }

  //compute G(iomega_n) by Dyson eq.
  alps::numeric::matrix<std::complex<double> > g0_mat(n_site, n_site, 0.0), s_mat(n_site, n_site, 0.0);//tmp arrays for gemm
  alps::numeric::matrix<std::complex<double> > g_mat(n_site, n_site, 0.0);
  for (std::size_t w = 0; w < n_matsubara; ++w) {
    for (std::size_t z = 0; z < n_flavors; ++z) {

      //prepare a matrix for S and bare G
      for (std::size_t site1 = 0; site1 < n_site; ++site1) {
        for (std::size_t site2 = 0; site2 < n_site; ++site2) {
          s_mat(site1, site2) = Sw[w][site1][site2][z];
          g0_mat(site1, site2) = bare_green_matsubara(w, site1, site2, z);
          g_mat(site1, site2) = 0.0;//just for safety
        }
      }

      //solve Dyson equation
      alps::numeric::gemm(g0_mat, s_mat, g_mat);
      g_mat += g0_mat;

      //write back full G
      for (std::size_t site1 = 0; site1 < n_site; ++site1) {
        for (std::size_t site2 = 0; site2 < n_site; ++site2) {
          sigma_green_matsubara_measured(w, site1, site2, z) = s_mat(site1, site2);
          green_matsubara_measured(w, site1, site2, z) = g_mat(site1, site2);
        }
      }
    }
  }
}


double green0_spline(const itime_green_function_t &green0, const itime_t delta_t,
                                              const int s1, const int s2, const spin_t z, int n_tau, double beta);


