#include <boost/random.hpp>

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

#include"legendre.h"
#include"util.h"
#include"fastupdate_formula.h"

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
