//
// Created by Hiroshi Shinaoka on 2015/10/15.
//

#pragma once

#include <math.h>
#include <vector>
#include <valarray>
#include <complex>

#include <Eigen/LU>

#include <alps/numeric/real.hpp>
#include "matrix.hpp"

namespace alps {
    namespace ctint {

        template<typename T> T mycast(std::complex<double> val);

        template<>
        double mycast(std::complex<double> val) {
            return val.real();
        }

        template<>
        std::complex<double> mycast(std::complex<double> val) {
            return val;
        }

        template<typename T> T myconj(T val);

        template<>
        double myconj(double val) {
            return val;
        }

        template<>
        std::complex<double> myconj(std::complex<double> val) {
            return std::conj(val);
        }

        template<typename T>
        T mysign(T x) {
          return x/std::abs(x);
        }


        inline double permutation(size_t N, size_t k) {
          assert(k>0);
          double r=1.0;
          for(size_t i=N-k+1; i<=N; ++i) {
            r *= static_cast<double>(i);
          }
          return r;
        }


        template<class R>
        std::vector<int> pickup_a_few_numbers(int N, int n, R& random01) {
            std::vector<int> flag(N,0), list(n);

            for (int i=0; i<n; ++i) {
                int itmp = 0;
                while(true) {
                    itmp = static_cast<int>(random01()*N);
                    assert(itmp<flag.size());
                    if (flag[itmp]==0) {
                        break;
                    }
                }
                assert(i<list.size());
                assert(itmp<flag.size());
                list[i] = itmp;
                flag[itmp] = 1;
            }
            return list;
        }


//very crapy way to remove elements at given positions from a std::vector
        template<class V>
        void remove_elements_from_vector(V& vec, std::vector<int> elements_to_be_removed) {
            size_t ipos = 0;
            for(typename V::iterator it = vec.begin(); it != vec.end();) {
                if(std::find(elements_to_be_removed.begin(), elements_to_be_removed.end(), ipos)!=elements_to_be_removed.end()) {
                    it = vec.erase(it);
                    ++ipos;
                } else {
                    ++it;
                    ++ipos;
                }
            }
        }

    }
}


/**
 * Just for compatibility. To be removed in future version
 */

namespace alps {
    namespace numeric {
        template <typename MatrixA, typename MatrixB>
        void my_copy_block(MatrixA const& A, typename MatrixA::size_type ai, typename MatrixA::size_type aj,
                           MatrixB& B, typename MatrixB::size_type bi, typename MatrixB::size_type bj,
                           typename MatrixA::difference_type m, typename MatrixA::difference_type n);


        template<typename T>
        size_t num_cols(const alps::numeric::submatrix_view<T>& matrix) {
          return matrix.num_cols();
        }

        template<typename T>
        size_t num_rows(const alps::numeric::submatrix_view<T>& matrix) {
          return matrix.num_rows();
        }


        //it's modified from alps::numeric::copy_block so that it accepcts a submatrix view
        template <typename MatrixA, typename MatrixB>
        void my_copy_block(MatrixA const& A, typename MatrixA::size_type ai, typename MatrixA::size_type aj,
                           MatrixB& B, typename MatrixB::size_type bi, typename MatrixB::size_type bj,
                           typename MatrixA::difference_type m, typename MatrixA::difference_type n)
        {
          assert(num_cols(B) >= bj+n);
          assert(num_rows(B) >= bi+m);
          for(typename MatrixA::difference_type j=0; j<n; ++j)
            std::copy(A.col(aj+j).first + ai, A.col(aj+j).first + ai + m, B.col(bj+j).first + bi);
        }

        template <typename T>
        typename alps::numeric::real_type<T>::type norm_square(const alps::numeric::submatrix_view<T>& M){
          using alps::numeric::real;
          typename alps::numeric::real_type<T>::type ret(0);
          for (std::size_t c = 0; c < num_cols(M); ++c)
            for (std::size_t r = 0; r < num_rows(M); ++r)
              ret += real(conj(M(r,c)) * M(r,c));
          return ret;
        }
    }
}

//
// An adaptor for the matrix to the boost::numeric::bindings
//

template<class T>
void blas_swap_cols(alps::numeric::matrix<T>& mat, int i, int j) {
  mat.swap_col(i, j);
}

template<class T>
void blas_swap_rows(alps::numeric::matrix<T>& mat, int i, int j) {
  mat.swap_row(i, j);
}

template<class T, class InputIterator>
void swap_cols_rows(alps::numeric::matrix<T>& mat, InputIterator first, InputIterator end) {
  for (InputIterator it=first; it!=end; ++it) {
    mat.swap_row_col(it->first, it->second);
  }
}

//double mymod(double x, double beta);

template<class T> alps::numeric::matrix<T>
mygemm(const alps::numeric::matrix<T>& A, const alps::numeric::matrix<T>& B) {
  alps::numeric::matrix<T> AB(A.size1(), B.size2(), 0.0);
  alps::fastupdate::gemm(A, B, AB);
  return AB;
}

template<typename MatrixA, typename MatrixB, typename MatrixC>
inline void mygemm(const typename MatrixA::type alpha,
                   const MatrixA& a, const MatrixB& b,
                   const typename MatrixA::type beta,
                   MatrixC& c ) {
  assert(a.size1()==c.size1());
  assert(a.size2()==b.size1());
  assert(b.size2()==c.size2());
  assert(a.size1()>0);
  assert(a.size2()>0);
  assert(b.size1()>0);
  assert(b.size2()>0);

  c = alpha * a * b + beta * c;
}

template<class T>
double
is_identity(const alps::numeric::matrix<T>& M) {
  if (num_cols(M)!=num_rows(M))
    throw std::runtime_error("num_rows != num_cols");

  const int Nv = num_cols(M);
  double max_diff = 0;
  for (size_t q=0; q<Nv; ++q) {
    for (size_t p=0; p<Nv; ++p) {
      if (p == q) {
        max_diff = std::max(max_diff, std::abs(M(p, q) - 1.));
      } else {
        max_diff = std::max(max_diff, std::abs(M(p, q)));
      }
    }
  }
  return max_diff;
}

template<class T>
bool is_all_zero(const std::valarray<T>& array) {
  bool flag = true;
  for (int i=0; i<array.size(); ++i) {
    if (array[i]!=0) {
      flag = false;
      break;
    }
  }
  return flag;
}

//norm_dist: distribution function normalized for [0,beta]
template<class P, class R>
double gen_rand_rejection_method(const P& norm_dist, double max_norm_dist, R& random01, double beta) {
  double inv_max_norm_dist = 1/max_norm_dist;
  double x;
  while(true) {
    x = beta*random01();
    //std::cout << " debug " << x << " " << random01() << " " <<  max_norm_dist << std::endl;
    const double Px = norm_dist(x);
    assert(Px<=max_norm_dist);
    if (random01() < inv_max_norm_dist*Px) {
      break;
    }
  }
  return x;
};


template<typename T>
bool my_equal(T x, T y, double eps=1E-8) {
  return std::abs(x-y)/std::max(std::abs(x),std::abs(y))<eps;
}

template<typename T>
bool my_rdiff(T x, T y) {
  return std::abs(x-y)/std::max(std::abs(x),std::abs(y));
}

inline
double mymod(double x, double beta) {
  if (x>=0) {
    return x-beta*static_cast<int>(x/beta);
  } else {
    return x+beta*(static_cast<int>(-x/beta)+1);
  }
}

