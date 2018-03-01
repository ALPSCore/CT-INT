#pragma once

#include <alps/fastupdate/resizable_matrix.hpp>

namespace alps {
    namespace numeric {
        template<typename T>
        using matrix = alps::fastupdate::ResizableMatrix<T>;

        template<typename T>
        using submatrix_view = Eigen::Block<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >;

        /*
        template<typename T>
        int num_rows(const matrix<T>& m) {
          m.size1();
        }

        template<typename T>
        int num_cols(const matrix<T>& m) {
          m.size2();
        }
         */
    }
}