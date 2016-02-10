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

  T operator() (annihilator& c_op, creator& cdagg_op) {
    if (c_op.s()!=cdagg_op.s())
      return 0.0;
    if (c_op.flavor()!=cdagg_op.flavor())
      return 0.0;

    //std::cout << "--- " <<  " ---- " << std::endl;
    return operator()(c_op.t().time()-cdagg_op.t().time());
  }

  double beta_;
};

