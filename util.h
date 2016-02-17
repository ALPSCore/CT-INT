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
#include <boost/numeric/bindings/stride.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

template<typename T> T mycast(std::complex<double> val);
template<typename T> T myconj(T val);
template<typename T> T mysign(T x);

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


namespace alps {
    namespace numeric {
        template <typename MatrixA, typename MatrixB>
        void my_copy_block(MatrixA const& A, typename MatrixA::size_type ai, typename MatrixA::size_type aj,
                           MatrixB& B, typename MatrixB::size_type bi, typename MatrixB::size_type bj,
                           typename MatrixA::difference_type m, typename MatrixA::difference_type n);

        template<class Matrix>
        typename Matrix::value_type determinant(const Matrix& M) {
            assert(num_rows(M)==num_cols(M));
            std::vector<int> ipiv(num_rows(M));
            const size_t N = num_rows(M);

            if (N==0) {
                return 1.0;
            } else if (N==1) {
                return M(0,0);
            } else if (N==2) {
                return M(0,0)*M(1,1)-M(0,1)*M(1,0);
            }

            Matrix M_copy(M);

            int info = boost::numeric::bindings::lapack::getrf(M_copy, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            typename Matrix::value_type det = 1.0;
            for (size_t i=0; i<N; ++i) {
                det *= M_copy(i,i);
            }
            int p = 0;
            for (size_t i=0; i<N-1; ++i) {
                if (ipiv[i] != i+1) {
                    ++p;
                }
            }
            return (p%2==0 ? det : -det);
        }

        template<class Matrix, class Matrix2>
        typename Matrix::value_type determinant_no_copy(const Matrix& M, Matrix2& workspace) {
            assert(num_rows(M)==num_cols(M));
            std::vector<int> ipiv(num_rows(M));
            const size_t N = num_rows(M);

            if (N==0) {
                return 1.0;
            } else if (N==1) {
                return M(0,0);
            } else if (N==2) {
                return M(0,0)*M(1,1)-M(0,1)*M(1,0);
            }

            assert(num_rows(workspace)==N);
            assert(num_cols(workspace)==N);
            my_copy_block(M,0,0,workspace,0,0,N,N);

            int info = boost::numeric::bindings::lapack::getrf(workspace, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            typename Matrix::value_type det = 1.0;
            for (size_t i=0; i<N; ++i) {
                det *= workspace(i,i);
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

            if (N==0) {
                return 1.0;
            } else if (N==1) {
                return M(0,0);
            } else if (N==2) {
                return M(0,0)*M(1,1)-M(0,1)*M(1,0);
            }

            double norm = std::sqrt(norm_square(M)/(N*N));
            M /= norm;

            int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            typename Matrix::value_type det = 1.0;
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

        template<class Matrix>
        typename Matrix::value_type sign_determinant(Matrix M) {
            std::vector<int> ipiv(num_rows(M));
            const int N = num_rows(M);

            if (N==0) {
                return 1.0;
            } else if (N==1) {
                return mysign<typename Matrix::value_type>(M(0,0));
            } else if (N==2) {
                return mysign<typename Matrix::value_type>(M(0,0)*M(1,1)-M(0,1)*M(1,0));
            }

            double norm = std::sqrt(norm_square(M)/(N*N));
            M /= norm;

            int info = boost::numeric::bindings::lapack::getrf(M, ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            typename Matrix::value_type sign_det = 1.0;
            for (size_t i=0; i<N; ++i) {
                sign_det *= mysign<typename Matrix::value_type>(M(i,i));
            }
            int p = 0;
            for (size_t i=0; i<N-1; ++i) {
                if (ipiv[i] != i+1) {
                    ++p;
                }
            }
            return (p%2==0 ? sign_det : -sign_det);
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

        template<class Matrix>
        void inverse_in_place(Matrix& M)
        {
            std::vector<int> ipiv(num_rows(M));

            int info = boost::numeric::bindings::lapack::getrf(M,ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRF !");

            info = boost::numeric::bindings::lapack::getri(M,ipiv);
            if (info != 0)
                throw std::runtime_error("Error in GETRI !");
        }

        //a view for a submatrix of alps::matrix<T> for using gemm in blas
        //Assumed column major order
        template<typename T>
        class submatrix_view {
        public:
            typedef T value_type;
            typedef typename alps::numeric::matrix<T>::size_type size_type;
            typedef typename alps::numeric::matrix<T>::difference_type difference_type;

            //types for iterators
            typedef value_type*
                col_element_iterator;                         ///< Iterator to iterate through the elements of a columns of the matrix
            typedef value_type const*
                const_col_element_iterator;                   ///< Const version of col_element_iterator

            submatrix_view(const alps::numeric::matrix<T>& matrix, std::size_t start_row, std::size_t start_col, std::size_t num_rows, std::size_t num_cols)
                : value(const_cast<T*>(&matrix(start_row,start_col))),
                  num_rows_(num_rows),
                  num_cols_(num_cols),
                  stride1_(matrix.stride1()),
                  stride2_(matrix.stride2()) {
                assert(matrix.stride1()==1);
                assert(matrix.num_rows()>=start_row+num_rows);
                assert(matrix.num_cols()>=start_col+num_cols);
                assert(start_row>=0);
                assert(start_col>=0);
            }

            T* get_value() const {
                return value;
            }

            size_t stride1() const {
                return stride1_;
            }

            size_t stride2() const {
                return stride2_;
            }

            size_t num_cols() const {
                return num_cols_;
            }

            size_t num_rows() const {
                return num_rows_;
            }

            value_type& operator()(const std::size_t i, const std::size_t j) {
                assert(i<num_rows());
                assert(j<num_cols());
                return *(value+i+stride2_*j);
            }

            value_type operator()(const std::size_t i, const std::size_t j) const {
                assert(i<num_rows());
                assert(j<num_cols());
                return *(value+i+stride2_*j);
            }

            std::pair<col_element_iterator,col_element_iterator> col(size_type col = 0 )
            {
                assert(col < num_cols());
                return std::make_pair(col_element_iterator(value+stride2_*col), col_element_iterator(value+stride2_*col+num_rows_));
            }

            std::pair<const_col_element_iterator,const_col_element_iterator> col(size_type col = 0 ) const
            {
                assert(col < num_cols());
                return std::make_pair(col_element_iterator(value+stride2_*col), col_element_iterator(value+stride2_*col+num_rows_));
            }

        private:
            T* const value;
            std::size_t num_rows_, num_cols_, stride1_, stride2_;
        };

        template<typename T>
        size_t num_cols(const submatrix_view<T>& matrix) {
            return matrix.num_cols();
        }

        template<typename T>
        size_t num_rows(const submatrix_view<T>& matrix) {
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
        typename real_type<T>::type norm_square(const submatrix_view<T>& M){
            using alps::numeric::real;
            typename real_type<T>::type ret(0);
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

namespace boost { namespace numeric { namespace bindings { namespace detail {

    template <typename T, typename Id, typename Enable>
    struct adaptor< ::alps::numeric::submatrix_view<T>, Id, Enable>
    {
        typedef typename copy_const< Id, T >::type              value_type;
        // TODO: fix the types of size and stride -> currently it's a workaround, since std::size_t causes problems with boost::numeric::bindings
        //typedef typename ::alps::numeric::matrix<T,Alloc>::size_type         size_type;
        //typedef typename ::alps::numeric::matrix<T,Alloc>::difference_type   difference_type;
        typedef std::ptrdiff_t  size_type;
        typedef std::ptrdiff_t  difference_type;

        typedef mpl::map<
            mpl::pair< tag::value_type,      value_type >,
            mpl::pair< tag::entity,          tag::matrix >,
            mpl::pair< tag::size_type<1>,    size_type >,
            mpl::pair< tag::size_type<2>,    size_type >,
            mpl::pair< tag::data_structure,  tag::linear_array >,
            mpl::pair< tag::data_order,      tag::column_major >,
            mpl::pair< tag::data_side,       tag::upper >,
            mpl::pair< tag::stride_type<1>,  tag::contiguous >,
            mpl::pair< tag::stride_type<2>,  difference_type >
        > property_map;

        static size_type size1( const Id& id ) {
            return id.num_rows();
        }

        static size_type size2( const Id& id ) {
            return id.num_cols();
        }

        static value_type* begin_value( Id& id ) {
            return id.get_value();
        }

        //static value_type* end_value( Id& id ) {
            //return &(*(id.col(id.num_cols()-1).second-1));
        //}

        static difference_type stride1( const Id& id ) {
            return id.stride1();
        }

        static difference_type stride2( const Id& id ) {
            return id.stride2();
        }

    };

 }}}}

template<class T, class InputIterator>
void swap_cols_rows(alps::numeric::matrix<T>& mat, InputIterator first, InputIterator end) {
    for (InputIterator it=first; it!=end; ++it) {
        mat.swap_cols(it->first, it->second);
        mat.swap_rows(it->first, it->second);
    }
}

double mymod(double x, double beta);

template<class T> alps::numeric::matrix<T>
mygemm(const alps::numeric::matrix<T>& A, const alps::numeric::matrix<T>& B) {
    alps::numeric::matrix<T> AB(alps::numeric::num_rows(A), alps::numeric::num_cols(B), 0.0);
    assert(A.num_cols()==B.num_rows());
    alps::numeric::gemm(A, B, AB);
    return AB;
}

template<typename MatrixA, typename MatrixB, typename MatrixC>
inline void mygemm( const typename boost::numeric::bindings::value_type< MatrixA >::type alpha,
       const MatrixA& a, const MatrixB& b,
       const typename boost::numeric::bindings::value_type< MatrixA >::type beta,
       MatrixC& c ) {
    assert(a.num_rows()==c.num_rows());
    assert(a.num_cols()==b.num_rows());
    assert(b.num_cols()==c.num_cols());
    assert(a.num_rows()>0);
    assert(a.num_cols()>0);
    assert(b.num_rows()>0);
    assert(b.num_cols()>0);
    boost::numeric::bindings::blas::gemm(alpha, a, b, beta, c);
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
T mysign(T x) {
    return x/std::abs(x);
}

template<typename T>
bool my_equal(T x, T y, double eps=1E-8) {
    return std::abs(x-y)/std::max(std::abs(x),std::abs(y))<eps;
}

template<typename T>
bool my_rdiff(T x, T y) {
    return std::abs(x-y)/std::max(std::abs(x),std::abs(y));
}


#endif //IMPSOLVER_UTIL_H

