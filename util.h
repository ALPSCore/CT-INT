//
// Created by Hiroshi Shinaoka on 2015/10/15.
//

#ifndef IMPSOLVER_UTIL_H
#define IMPSOLVER_UTIL_H

#include <math.h>
#include <vector>
#include <complex>

#include <alps/numeric/matrix.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

template<typename T> T mycast(std::complex<double> val);
template<typename T> T myconj(T val);

//typedef alps::numeric::matrix<double> dense_matrix;
//typedef alps::numeric::matrix<std::complex<double> > complex_dense_matrix;
//
//template<class MAT, class SCALAR>
//inline void invert(MAT & A, SCALAR & det);

//return permutation P(N, k)
double permutation(size_t N, size_t k);

template<class R>
std::vector<int> pickup_a_few_numbers(int N, int n, R& random01) {
    std::vector<int> flag(N,0), list(n);

    for (int i=0; i<n; ++i) {
        int itmp = 0;
        while(true) {
            itmp = static_cast<int>(random01()*N);
            if (flag[itmp]==0) {
                break;
            }
        }
        list[i] = itmp;
        flag[itmp] = 1;
    }
    return list;
}


//very crapy way to remove elements at given positions from a std::vector
template<class T>
void remove_elements_from_vector(std::vector<T>& vec, std::vector<int> elements_to_be_removed) {
    size_t ipos = 0;
    for(typename std::vector<T>::iterator it = vec.begin(); it != vec.end();) {
        if(std::find(elements_to_be_removed.begin(), elements_to_be_removed.end(), ipos)!=elements_to_be_removed.end()) {
            it = vec.erase(it);
            ++ipos;
        } else {
            ++it;
            ++ipos;
        }
    }
}


namespace alps {
    namespace numeric {

        template<class Matrix>
        typename Matrix::value_type determinant(Matrix M) {
            std::vector<int> ipiv(num_rows(M));
            const size_t N = num_rows(M);

            if (N==0) {
                return 1.0;
            }

            int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            double det = 1.0;
            for (size_t i=0; i<N; ++i) {
                det *= M(i,i);
            }
            int p = 0;
            for (size_t i=0; i<N-1; ++i) {
                if (ipiv[i] != i+1) {
                    ++p;
                }
            }
            return (p%2==0 ? det : -det);
        }

        template<class Matrix>
        typename Matrix::value_type safe_determinant(Matrix M) {
            std::vector<int> ipiv(num_rows(M));
            const int N = num_rows(M);

            double norm = std::sqrt(norm_square(M)/(N*N));
            M /= norm;

            if (N==0) {
                return 1.0;
            }

            int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            double det = 1.0;
            for (size_t i=0; i<N; ++i) {
                det *= M(i,i);
            }
            int p = 0;
            for (size_t i=0; i<N-1; ++i) {
                if (ipiv[i] != i+1) {
                    ++p;
                }
            }
            double norm_pow = std::pow(norm, static_cast<double>(N));
            return (p%2==0 ? det*norm_pow : -det*norm_pow);
        }


        template<class T> alps::numeric::matrix<T>
        safe_inverse(const alps::numeric::matrix<T>& A) {
            using namespace alps::numeric;
            assert(num_rows(A)==num_cols(A));
            const int N = num_rows(A);
            double norm = std::sqrt(alps::numeric::norm_square(A)/(N*N));
            alps::numeric::matrix<T> A_norm(A);
            A_norm /= norm;
            A_norm = alps::numeric::inverse(A_norm);
            A_norm /= norm;
            return A_norm;
        }
    }
}

//template<class T>
//void swap_cols_rows(alps::numeric::matrix<T>& mat, const std::vector<std::pair<int,int> >& swap_list) {
    //for (int i=0; i<swap_list.size(); ++i) {
        //mat.swap_cols(swap_list[i].first, swap_list[i].second);
        ////mat.swap_rows(swap_list[i].first, swap_list[i].second);
    //}
//}

template<class T, class InputIterator>
void swap_cols_rows(alps::numeric::matrix<T>& mat, InputIterator first, InputIterator end) {
    for (InputIterator it=first; it!=end; ++it) {
        //mat.swap_cols(swap_list[i].first, swap_list[i].second);
        //mat.swap_rows(swap_list[i].first, swap_list[i].second);
        mat.swap_cols(it->first, it->second);
        mat.swap_rows(it->first, it->second);
    }
}

double mymod(double x, double beta);

template<class T> alps::numeric::matrix<T>
mygemm(const alps::numeric::matrix<T>& A, const alps::numeric::matrix<T>& B) {
    alps::numeric::matrix<T> AB(alps::numeric::num_rows(A), alps::numeric::num_cols(B), 0.0);
    alps::numeric::gemm(A, B, AB);
    return AB;
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


#endif //IMPSOLVER_UTIL_H
