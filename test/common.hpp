#include <boost/random.hpp>

#include "../src/matrix.hpp"

template<class T>
void randomize_matrix(alps::numeric::matrix<T>& mat, size_t seed=100) {
    boost::random::mt19937 gen;
    boost::random::uniform_01<double> dist;
    gen.seed(seed);

    for (int j=0; j<num_cols(mat); ++j) {
        for (int i=0; i<num_rows(mat); ++i) {
            mat(i,j) = dist(gen);
        }
    }
}
