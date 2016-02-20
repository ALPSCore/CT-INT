#include <complex>
#include <limits>

#include <boost/random.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <boost/range/irange.hpp>

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

#include"legendre.h"
#include"util.h"
#include"fastupdate_formula.h"
#include"green_function.h"
#include"update_statistics.h"
#include"operator.hpp"
#include"submatrix.hpp"
#include "update_manager.hpp"

template<class T>
void randomize_matrix(alps::numeric::matrix<T>& mat, size_t seed=100) {
    boost::random::mt19937 gen;
    boost::random::uniform_01<double> dist;
    gen.seed(seed);

    //std::cout << " debug " << num_cols(mat) << " " << num_rows(mat) << std::endl;
    for (int j=0; j<num_cols(mat); ++j) {
        for (int i=0; i<num_rows(mat); ++i) {
            mat(i,j) = dist(gen);
        }
    }
}


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
