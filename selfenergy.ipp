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
  const M_TYPE sign = submatrix_update->sign();

  //clear contents of Wk
  std::fill(Wk.origin(), Wk.origin()+Wk.num_elements(), 0);

  //allocate memory for work
  size_t max_mat_size = 0;
  for (unsigned int z=0; z<n_flavors; ++z) {
    max_mat_size = std::max(max_mat_size, M_flavors[z].num_cols());
  }
  alps::numeric::matrix<std::complex<double> > GR(max_mat_size, n_site),
      GL(n_site, max_mat_size), MGR(max_mat_size, n_site), GLMGR(max_mat_size, max_mat_size);

  //Note: vectorized for the loops over operators (should be optimal at large beta and U)
  // if computing sin and cos is too expensive, consider using linear interpolation.
  for (size_t z=0; z<n_flavors; ++z) {
    const size_t Nv = num_rows(M_flavors[z]);

    if (Nv==0) {
      continue;
    }

    const std::vector<annihilator>& annihilators = submatrix_update->invA()[z].annihilators();
    const std::vector<creator>& creators = submatrix_update->invA()[z].creators();

    for(unsigned int i_freq=0; i_freq <nfreq; ++i_freq) {
      GR.resize(Nv, n_site);
      MGR.resize(Nv, n_site);
      GL.resize(n_site, Nv);
      GLMGR.resize(n_site, n_site);

      //GR
      for(unsigned int p=0;p<Nv;++p) {
        const size_t site_p = annihilators[p].s();
        const double phase = annihilators[p].t().time()*(2*i_freq+1)*M_PI/beta;
        const std::complex<double> exp = std::complex<double>(std::cos(phase), -std::sin(phase));
        for (size_t site=0; site<n_site; ++site) {
          GR(p, site) = bare_green_matsubara(i_freq, site_p, site, z)*exp;
        }
      }

      //GL
      for(unsigned int q=0;q<Nv;++q) {
        const size_t site_q = creators[q].s();
        const double phase = creators[q].t().time()*(2*i_freq+1)*M_PI/beta;
        const std::complex<double> exp = std::complex<double>(std::cos(phase), std::sin(phase));
        for (size_t site=0; site<n_site; ++site) {
          GL(site, q) = bare_green_matsubara(i_freq, site, site_q, z)*exp;
        }
      }

      //clear MGR, GLMGR
      std::fill(MGR.get_values().begin(), MGR.get_values().end(), 0);
      std::fill(GLMGR.get_values().begin(), GLMGR.get_values().end(), 0);

      gemm(M_flavors[z], GR, MGR);
      gemm(GL, MGR, GLMGR);

      for (unsigned int site1=0; site1<n_site; ++site1) {
        for (unsigned int site2=0; site2<n_site; ++site2) {
          Wk[z][site1][site2][i_freq] += GLMGR(site1, site2);
        }
      }

    }//i_freq
  }//z

  //int num_nd = 0;
  //for (int iv=0; iv<itime_vertices.size(); ++iv) {
    //if (!itime_vertices[iv].is_density_type())
      //++num_nd;
  //}

  if (params.defined("PREFIX_OUTPUT_TIME_SERIES")) {
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
        }
        measurements[Wk_real_name.str().c_str()] << Wk_real;
        measurements[Wk_imag_name.str().c_str()] << Wk_imag;
      }
    }
  }
}

template<class TYPES>
void InteractionExpansion<TYPES>::compute_Sl() {
  static boost::multi_array<std::complex<double>,3> Sl(boost::extents[n_site][n_site][n_legendre]);
  const size_t num_random_walk = 100;

  const M_TYPE sign = submatrix_update->sign();
  const double temperature = 1.0/beta;

  //Work arrays
  size_t max_mat_size = 0;
  for (unsigned int z=0; z<n_flavors; ++z) {
    max_mat_size = std::max(max_mat_size, M_flavors[z].num_cols());
  }
  //std::vector<double> legendre_vals(n_legendre), sqrt_vals(n_legendre);
  //std::vector<double> sqrt_vals(n_legendre);
  //for(unsigned int i_legendre=0; i_legendre<n_legendre; ++i_legendre) {
    //sqrt_vals[i_legendre] = std::sqrt(2.0*i_legendre+1.0);
  //}
  const std::vector<double>& sqrt_vals = legendre_transformer.get_sqrt_2l_1();

  std::vector<double> x_vals;
  boost::multi_array<double,2> legendre_vals_all; //, legendre_vals_trans_all;

  alps::numeric::matrix<M_TYPE> gR(max_mat_size, n_site), M_gR(max_mat_size, n_site);

  for (unsigned int z=0; z<n_flavors; ++z) {
    std::fill(Sl.origin(),Sl.origin()+Sl.num_elements(),0.0);//clear the content for safety
    const size_t Nv = M_flavors[z].num_cols();

    if (Nv==0) {
      continue;
    }
    gR.resize(Nv, n_site);
    M_gR.resize(Nv, n_site);

    x_vals.resize(Nv);
    legendre_vals_all.resize(boost::extents[n_legendre][Nv]);

    const std::vector<annihilator>& annihilators = submatrix_update->invA()[z].annihilators();
    const std::vector<creator>& creators = submatrix_update->invA()[z].creators();

    //shift times of operators by time_shift
    for (std::size_t random_walk=0; random_walk<num_random_walk; ++random_walk) {
      const double time_shift = beta * random();

      for (unsigned int p = 0; p < Nv; ++p) {//annihilation operators
        const double tmp = annihilators[p].t().time() + time_shift;
        const double time_a_shifted = tmp < beta ? tmp : tmp - beta;
        const double coeff = tmp < beta ? 1 : -1;

        //interpolate G0
        for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
          gR(p, site_B) = mycast<M_TYPE>(coeff*g0_intpl(time_a_shifted, z, annihilators[p].s(), site_B));
        }
      }

      gemm(M_flavors[z], gR, M_gR);

      //compute legendre coefficients
      for (unsigned int q = 0; q < Nv; ++q) {//creation operators
        const double tmp = creators[q].t().time() + time_shift;
        const double time_c_shifted = tmp < beta ? tmp : tmp - beta;
        x_vals[q] = 2 * time_c_shifted*temperature - 1.0;
      }
      legendre_transformer.compute_legendre(x_vals, legendre_vals_all);//P_l[x(tau_q)]

      for (unsigned int q = 0; q < Nv; ++q) {//creation operators
        const unsigned int site_c = creators[q].s();
        const double tmp = creators[q].t().time() + time_shift;
        const double coeff = tmp < beta ? 1 : -1;

        const int n_legendre_tmp = n_legendre;
        for (unsigned int site_B = 0; site_B < n_site; ++site_B) {
          for (unsigned int i_legendre = 0; i_legendre < n_legendre_tmp; ++i_legendre) {
            Sl[site_c][site_B][i_legendre] -= coeff*sqrt_vals[i_legendre] * legendre_vals_all[i_legendre][q] * M_gR(q, site_B);
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
  const M_TYPE sign = submatrix_update->sign();

  std::vector< std::vector<double> > dens(n_flavors);
  for(unsigned int z=0;z<n_flavors;++z){
    dens[z].resize(n_site);
    memset(&(dens[z][0]), 0., sizeof(double)*(n_site));
  }
  double tau = beta*random();
  double sign_real = mycast<double>(sign);
  for (unsigned int z=0; z<n_flavors; ++z) {
    const size_t Nv = M_flavors[z].num_cols();
    const std::vector<annihilator>& annihilators = submatrix_update->invA()[z].annihilators();
    const std::vector<creator>& creators = submatrix_update->invA()[z].creators();

    alps::numeric::vector<M_TYPE> g0_tauj(Nv);
    alps::numeric::vector<M_TYPE> M_g0_tauj(Nv);
    alps::numeric::vector<M_TYPE> g0_taui(Nv);
    for (unsigned int s=0;s<n_site;++s) {
      for (unsigned int j=0;j<Nv;++j)
        g0_tauj[j] = mycast<M_TYPE>(green0_spline_new(annihilators[j].t().time()-tau, z, annihilators[j].s(), s));//CHECK THE TREATMENT OF EQUAL-TIME Green's function
      for (unsigned int i=0;i<Nv;++i)
        g0_taui[i] = mycast<M_TYPE>(green0_spline_new(tau-creators[i].t().time(),z, s, creators[i].s()));
      if (M_flavors[z].num_cols()>0) {
        gemv(M_flavors[z],g0_tauj,M_g0_tauj);
      }
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
    measurements["densities_"+boost::lexical_cast<std::string>(z)] << static_cast<std::valarray<double> >(densmeas*sign_real);
    densities[z] /= n_site;
    densities[z] = densities[z];
  }
  measurements["densities"] << static_cast<std::valarray<double> > (densities*sign_real);

  if (n_flavors==2) {
    double density_correlation = 0.;
    for (unsigned int i=0; i<n_site; ++i) {
      density_correlation += (dens[0][i])*(dens[1][i]);
    }
    density_correlation /= n_site;
    measurements["density_correlation"] << (density_correlation*sign_real);
  }

  {
    std::valarray<double> ninj(n_site*n_site*n_flavors*n_flavors);
    int pos = 0;
    for (unsigned int i=0; i<n_site; ++i) {
      for (unsigned int j=0; j<n_site; ++j) {
        for (unsigned int flavor1=0; flavor1<n_flavors; ++flavor1) {
          for (unsigned int flavor2=0; flavor2<n_flavors; ++flavor2) {
            ninj[pos] = (dens[flavor1][i])*(dens[flavor2][j]);
            ++pos;
          }
        }
      }
    }
    measurements["n_i n_j"] << static_cast<std::valarray<double> > (ninj*sign_real);
  }
}

#endif //IMPSOLVER_SELFENERGY_HPP
