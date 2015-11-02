export ALPS_ROOT=/opt/alps
export ALPS_SRC_DIR=/Users/hiroshi/work/src/alps
export PATH=$ALPS_ROOT/bin:$PATH
export ALPS_DIR=$ALPS_ROOT/share/alps
export PYTHONPATH=$ALPS_ROOT/lib:$PYTHONPATH
export BLAS_LAPACK_LIBRARIES="mkl_intel_lp64;mkl_core;mkl_sequential;pthread;m"
export BLAS_LAPACK_DIRS="/opt/intel/composer_xe_2013_sp1.1.103/mkl/lib"
export BLAS_LAPACK_INCLUDE_DIRS="/opt/intel/composer_xe_2013_sp1.1.103/mkl/include"
export CXX=openmpicxx
#cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-m64 -std=c++11 -stdlib=libstdc++"  -DCMAKE_INSTALL_PREFIX=/opt/ct-int /Users/hiroshi/git/interaction_expansion3
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-stdlib=libstdc++"  -DCMAKE_INSTALL_PREFIX=/opt/ct-int /Users/hiroshi/git/interaction_expansion3
