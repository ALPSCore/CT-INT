#ifndef IMPSOLVER_SELFENERGY_HPP
#define IMPSOLVER_SELFENERGY_HPP

#include "interaction_expansion.hpp"

template<class TYPES>
void InteractionExpansion<TYPES>::compute_W_matsubara()
{
  Wk_t Wk(boost::extents[n_flavors][n_site][n_site][n_matsubara]);
  measure_Wk(Wk, n_matsubara_measurements);
  measure_densities();
}

template<class TYPES>
void InteractionExpansion<TYPES>::measure_Wk(Wk_t& Wk, const unsigned int nfreq)
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

  //Note: vectorized for the loops over operators (should be optimal at large beta and U)
  // if computing sin and cos is too expensive, consider using linear interpolation.
  for (size_t z=0; z<n_flavors; ++z) {
    const size_t Nv = num_rows(M[z].matrix());

    if (Nv==0) {
      continue;
    }

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
    for (unsigned int flavor=0; flavor<n_flavors; ++flavor) {
      for (unsigned int site1 = 0; site1 < n_site; ++site1) {
        Wk_dynamics.push_back(Wk[flavor][site1][site1][0]*sign);
      }
    }
  }

  for(unsigned int flavor=0;flavor<n_flavors;++flavor) {
    for (unsigned int site1 = 0; site1 < n_site; ++site1) {
      for (unsigned int site2 = 0; site2 < n_site; ++site2) {
        std::stringstream Wk_real_name, Wk_imag_name;
        Wk_real_name << "Wk_real_" << flavor << "_" << site1 << "_" << site2;
        Wk_imag_name << "Wk_imag_" << flavor << "_" << site1 << "_" << site2;
        std::valarray<double> Wk_real(nfreq);
        std::valarray<double> Wk_imag(nfreq);
        for (unsigned int w = 0; w < nfreq; ++w) {
          std::complex<double> ztmp = Wk[flavor][site1][site2][w]*sign;
          Wk_real[w] = ztmp.real();
          Wk_imag[w] = ztmp.imag();
          //Wk_real[w] = Wk[flavor][site1][site2][w].real();
          //Wk_imag[w] = Wk[flavor][site1][site2][w].imag();
        }
        measurements[Wk_real_name.str().c_str()] << Wk_real;
        measurements[Wk_imag_name.str().c_str()] << Wk_imag;
        //measurements[Wk_real_name.str().c_str()] << static_cast<std::valarray<double> > (Wk_real);
        //measurements[Wk_imag_name.str().c_str()] << static_cast<std::valarray<double> > (Wk_imag);
      }
    }
  }
}

template<class TYPES>
void InteractionExpansion<TYPES>::compute_Sl() {
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
  alps::numeric::matrix<M_TYPE> gR(max_mat_size, n_site), M_gR(max_mat_size, n_site);

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
          gR(p, site_B) = mycast<M_TYPE>(coeff*green0_spline_new(time_a_shifted, z, M[z].annihilators()[p].s(), site_B));
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
          const std::complex<double> ztmp = (Sl[site1][site2][i_legendre]*sign)/static_cast<double>(num_random_walk);
          Sl_real[i_legendre] = ztmp.real();
          Sl_imag[i_legendre] = ztmp.imag();
          //Sl_real[i_legendre] = Sl[site1][site2][i_legendre].real()/num_random_walk;
          //Sl_imag[i_legendre] = Sl[site1][site2][i_legendre].imag()/num_random_walk;
        }
        measurements[Sl_real_name.str().c_str()] << Sl_real;
        measurements[Sl_imag_name.str().c_str()] << Sl_imag;
      }//site2
    }//site1

  }//z
}


template<class TYPES>
void InteractionExpansion<TYPES>::measure_densities()
{
  std::vector< std::vector<double> > dens(n_flavors);
  for(unsigned int z=0;z<n_flavors;++z){
    dens[z].resize(n_site);
    memset(&(dens[z][0]), 0., sizeof(double)*(n_site));
  }
  double tau = beta*random();
  double sign_real = mycast<double>(sign);
  for (unsigned int z=0; z<n_flavors; ++z) {
    const size_t Nv = num_rows(M[z].matrix());
    alps::numeric::vector<M_TYPE> g0_tauj(Nv);
    alps::numeric::vector<M_TYPE> M_g0_tauj(Nv);
    alps::numeric::vector<M_TYPE> g0_taui(Nv);
    for (unsigned int s=0;s<n_site;++s) {
      for (unsigned int j=0;j<Nv;++j)
        g0_tauj[j] = mycast<M_TYPE>(green0_spline_new(M[z].annihilators()[j].t().time()-tau, z, M[z].annihilators()[j].s(), s));//CHECK THE TREATMENT OF EQUAL-TIME Green's function
      for (unsigned int i=0;i<Nv;++i)
        g0_taui[i] = mycast<M_TYPE>(green0_spline_new(tau-M[z].creators()[i].t().time(),z, s, M[z].creators()[i].s()));
      if (num_rows(M[z].matrix())>0)
        gemv(M[z].matrix(),g0_tauj,M_g0_tauj);
      dens[z][s] += mycast<double>(green0_spline_new(-beta*1E-10,z,s,s));//tau=-0
      for (unsigned int j=0;j<Nv;++j)
        dens[z][s] -= mycast<double>(g0_taui[j]*M_g0_tauj[j]);
    }
  }
  std::valarray<double> densities(0., n_flavors);
  for (unsigned int z=0; z<n_flavors; ++z) {
    std::valarray<double> densmeas(n_site);
    for (unsigned int i=0; i<n_site; ++i) {
      densities[z] += dens[z][i];
      densmeas[i] = dens[z][i];
    }
    //measurements["densities_"+boost::lexical_cast<std::string>(z)] << static_cast<std::valarray<double> > (densmeas*sign);
    measurements["densities_"+boost::lexical_cast<std::string>(z)] << static_cast<std::valarray<double> >(densmeas*sign_real);
    densities[z] /= n_site;
    densities[z] = densities[z];
  }
  measurements["densities"] << static_cast<std::valarray<double> > (densities*sign_real);
  double density_correlation = 0.;
  for (unsigned int i=0; i<n_site; ++i) {
    density_correlation += (dens[0][i])*(dens[1][i]);
  }
  density_correlation /= n_site;
  measurements["density_correlation"] << (density_correlation*sign_real);
  std::valarray<double> ninj(n_site*n_site*4);
  for (unsigned int i=0; i<n_site; ++i) {
    for (unsigned int j=0; j<n_site; ++j) {
      ninj[i*n_site+j] = (dens[0][i])*(dens[0][j]);
      ninj[i*n_site+j+1] = (dens[0][i])*(dens[1][j]);
      ninj[i*n_site+j+2] = (dens[1][i])*(dens[0][j]);
      ninj[i*n_site+j+3] = (dens[1][i])*(dens[1][j]);
    }
  }
  measurements["n_i n_j"] << static_cast<std::valarray<double> > (ninj*sign_real);
}

#endif //IMPSOLVER_SELFENERGY_HPP
