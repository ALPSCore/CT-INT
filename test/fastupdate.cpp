#include <algorithm>

#include "gtest.h"
#include "common.hpp"

#include <alps/params.hpp>

#include "../src/submatrix.hpp"
#include "../src/operator.hpp"
#include "../src/update_manager.hpp"
#include "../src/program_options.hpp"

using namespace alps::ctint;

class simple_vertex {
public:
    simple_vertex(double time) : time_(time) {}
    double time() const {return time_;}

private:
    double time_;
};


//Interpolation of G0
template<typename T>
struct DiagonalG0 {
    DiagonalG0 (double beta) : beta_(beta) {}

    int nsite() const {return 1;}
    int nflavor() const {return 2;}
    bool is_zero(int site1, int site2, int flavor, double eps) const {return site1!=site2;}

    T operator() (double dt) const {
      //const double dt = c_op.t().time()-cdagg_op.t().time();
      if (dt>=0 && dt<beta_) {
        return -(1.0 / (beta_ * beta_)) * (dt - 0.5 * beta_) * (dt - 0.5 * beta_) - 0.25;
      } else if (dt>=beta_ && dt<2*beta_) {
        double dt2 = dt-beta_;
        return -(-(1.0 / (beta_ * beta_)) * (dt2 - 0.5 * beta_) * (dt2 - 0.5 * beta_) - 0.25);
      } else if (dt>=-beta_ && dt<0) {
        double dt2 = dt+beta_;
        return -(-(1.0 / (beta_ * beta_)) * (dt2 - 0.5 * beta_) * (dt2 - 0.5 * beta_) - 0.25);
      }  else {
        throw std::runtime_error("Digonaal G0: Not implemented!");
      }
    }

    T operator() (const annihilator& c_op, const creator& cdagg_op) {
      if (c_op.s()!=cdagg_op.s())
        return 0.0;
      if (c_op.flavor()!=cdagg_op.flavor())
        return 0.0;

      //std::cout << "--- " <<  " ---- " << std::endl;
      double dt = c_op.t().time()-cdagg_op.t().time();
      if (std::abs(dt)<1E-10) {
        dt += c_op.t().small_index() > cdagg_op.t().small_index() ? 1E-8 : -1E-8;
      }
      return operator()(dt);
    }

    double beta_;
};

//Interpolation of G0
template<typename T>
struct OffDiagonalG0 {
    OffDiagonalG0 (double beta, int n_site, const std::vector<double>& E, const boost::multi_array<T,2>& phase) : beta_(beta), n_site_(n_site), E_(E), phase_(phase) {}

    int nsite() const {return n_site_;}
    int nflavor() const {return 2;}
    bool is_zero(int site1, int site2, int flavor, double eps) const {return false;}

    T operator() (const annihilator& c_op, const creator& cdagg_op) {
      const double dt = c_op.t().time()-cdagg_op.t().time();
      double dt_tmp = dt;
      if (dt_tmp > beta_) dt_tmp -= beta_;
      if (dt_tmp < 0) dt_tmp += beta_;

      if (c_op.s()==cdagg_op.s()) {
        //return -(1.0 / (beta_ * beta_)) * (dt_tmp - 0.5 * beta_) * (dt_tmp - 0.5 * beta_) - 0.25 + (-dt+0.5*beta_)/(10*beta_);
        const double E_tmp = E_[c_op.s()];
        return  -std::exp((beta_-dt_tmp)*E_tmp)/(1.0+std::exp(beta_*E_tmp));
      } else {
        //std::cout << myconj(phase_[c_op.s()]) << " " << phase_[cdagg_op.s()] << myconj(phase_[c_op.s()])*phase_[cdagg_op.s()] << std::endl;
        return (-dt+0.5*beta_)/(2*beta_)*phase_[c_op.s()][cdagg_op.s()];
      }
    }

    double beta_;
    int n_site_;
    std::vector<double> E_;
    boost::multi_array<T,2> phase_;
};

TEST(FastUpdate, BlockMatrixReplaceRowsColsSingular) {
    typedef alps::numeric::matrix<double> matrix_t;

    for (int i=0; i<20; ++i) {
        double alpha = pow(0.1, 20*i);
        double a=4., b=8.;

        matrix_t G2(2,2), G2p(2,2), M2p(2,2);
        std::vector<std::pair<int, int> > swap_list;

        G2(0,0) = alpha;
        G2(1,1) = alpha;
        G2(1,0) = a;
        G2(0,1) = a;

        M2p(0,0) = alpha;
        M2p(1,0) = -b;
        M2p(0,1) = -b;
        M2p(1,1) = alpha;
        M2p.block() /= (alpha+b)*(alpha-b);

        matrix_t Q(1,1), R(1,1), S(1,1);
        Q(0,0) = b;
        R(0,0) = b;
        S(0,0) = alpha;

        matrix_t M2p_fast = inverse(G2);
        matrix_t Mmat, inv_tSp;
        matrix_t tPp, tQp, tRp, tSp;

    }
}

TEST(FastUpdate, ReplaceDiagonalElements) {
    typedef double T;
    typedef alps::numeric::matrix<T> matrix_t;

    const int N=10, m=2, offset=2;
    //const int N=2, m=1, offset=0;
    assert(m+offset<=N);

    matrix_t A_old(N,N), A_new(N,N), new_elems(m,1);
    std::vector<T> elems_diff(m);
    std::vector<int> pos(m);

    randomize_matrix(A_old, 100);
    randomize_matrix(new_elems, 200);
    A_new = A_old;
    matrix_t invA_old = inverse(A_old);
    for (int i=0; i<m; ++i) {
        pos[i] = i+offset;
    }
    for (int i=0; i<m; ++i) {
        elems_diff[i] = new_elems(i,0)-A_old(pos[i],pos[i]);
        A_new(pos[i],pos[i]) = new_elems(i,0);
    }

    const T det_rat = determinant(A_new)/determinant(A_old);
    const T det_rat_fast = compute_det_ratio_replace_diaognal_elements(invA_old, m, pos, elems_diff, true);
    ASSERT_TRUE(std::abs((det_rat-det_rat_fast)/det_rat)<1E-8);

    /* inverse matrix update */
    matrix_t invA_new = inverse(A_new);
    matrix_t invA_new_fast = invA_old;
    compute_det_ratio_replace_diaognal_elements(invA_new_fast, m, pos, elems_diff, false);

    ASSERT_TRUE(std::abs(alps::fastupdate::norm_square(invA_new-invA_new_fast))<1E-5);
}

TEST(SubmatrixUpdate, single_vertex_insertion_spin_flip)
{
  typedef std::complex<double> T;
  const int n_sites = 3;
  const double U = 2.0;
  const double alpha = 1E-2;
  const double beta = 200.0;
  const int Nv_max = 2;
  const int n_spins = 2;
  const int k_ins_max = 32;
  const int n_update = 5;
  const int seed = 100;

  std::vector<double> E(n_sites);
  boost::multi_array<T,2> phase(boost::extents[n_sites][n_sites]);

  for (int i=0; i<n_sites; ++i) {
    E[i] = (double) i;
    //std::cout << phase[i] << std::endl;
  }
  for (int i=0; i<n_sites; ++i) {
    for (int j=i; j<n_sites; ++j) {
      phase[i][j] = std::exp(std::complex<double>(0.0, 1.*i*(2*j+1.0)));
      phase[j][i] = myconj(phase[i][j]);
    }
  }

  general_U_matrix<T> Uijkl(n_sites, U, alpha);

  itime_vertex_container itime_vertices_init;
  itime_vertices_init.push_back(itime_vertex(0, 0, 0.5*beta, 2, true));

  /* initialize submatrix_update */
  //SubmatrixUpdate<T> submatrix_update(k_ins_max, n_spins, DiagonalG0<T>(beta), &Uijkl, beta, itime_vertices_init);
  SubmatrixUpdate<T> submatrix_update(k_ins_max, n_spins, OffDiagonalG0<T>(beta, n_sites, E, phase), &Uijkl, beta, itime_vertices_init);

  submatrix_update.sanity_check();

  /* init udpate_manager */
  alps::params params;
  define_ctint_options(params);
  params["model.beta"] = beta;
  params["model.spins"] = n_spins;
  params["update.n_multi_vertex_update"] = Nv_max;
  params["update.double_vertex_update_A"] = 1.0/beta;
  params["update.double_vertex_update_B"] = 1.0e-2;
  params["update.vertex_shift_step_size"] = 0.1*beta;
  VertexUpdateManager<T> manager(params, Uijkl, OffDiagonalG0<T>(beta, n_sites, E, phase), false);

  /* initialize RND generator */
  boost::random::uniform_smallint<> dist(1,Nv_max);
  boost::random::uniform_01<> dist01;
  boost::random::mt19937 gen(seed);
  boost::random::variate_generator<boost::random::mt19937&, boost::random::uniform_smallint<> > Nv_prob(gen, dist);
  boost::random::variate_generator<boost::random::mt19937&, boost::random::uniform_01<> > random01(gen, dist01);

  std::vector<alps::numeric::matrix<T> > M(n_spins), M_scratch(n_spins);

  for (int i_update=0; i_update<n_update; ++i_update) {
    T sign_from_M0, weight_from_M0;
    boost::tie(sign_from_M0,weight_from_M0) = submatrix_update.compute_M_from_scratch(M_scratch);

    const T weight_rat = manager.do_ins_rem_update(submatrix_update, Uijkl, random01, 1.0);
    const T sign_bak = submatrix_update.sign();

    ASSERT_TRUE(submatrix_update.sanity_check());
    submatrix_update.recompute_matrix(true);
    submatrix_update.compute_M(M);
    T sign_from_M, weight_from_M;
    boost::tie(sign_from_M,weight_from_M) = submatrix_update.compute_M_from_scratch(M_scratch);

    ASSERT_TRUE(my_equal(weight_from_M/weight_from_M0, weight_rat, 1E-5));

    ASSERT_TRUE(std::abs(sign_bak-submatrix_update.sign())<1.0e-5);
    ASSERT_TRUE(std::abs(sign_from_M-submatrix_update.sign())<1.0e-5);
    for (int flavor=0; flavor<n_spins; ++flavor) {
      if (M[flavor].size2()>0) {
        ASSERT_TRUE(alps::fastupdate::norm_square(M[flavor]-M_scratch[flavor])/alps::fastupdate::norm_square(M[flavor])<1E-8);
      }
    }

    const T weight_rat2 = manager.do_spin_flip_update(submatrix_update, Uijkl, random01);

    T sign_from_M2, weight_from_M2;
    boost::tie(sign_from_M2,weight_from_M2) = submatrix_update.compute_M_from_scratch(M_scratch);
    ASSERT_TRUE(my_equal(weight_from_M2/weight_from_M, weight_rat2, 1E-5));

  }
}
