#pragma once

#include <alps/fastupdate/resizable_matrix.hpp>

namespace alps {
    namespace numeric {
        template<typename T>
        using matrix = alps::fastupdate::ResizableMatrix<T>;

        template<typename T>
        using submatrix_view = Eigen::Block<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >;
    }
}