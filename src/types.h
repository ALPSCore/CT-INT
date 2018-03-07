#pragma once

#include <sys/types.h>
#include <complex>
#include <boost/cstdint.hpp>

namespace alps {
    namespace ctint {

///enum for spin up and spin down
        enum  {up=0, down=1} ;
///addressing type for site indices (cluster)
        typedef unsigned int site_t;
///addressing type for spin indices
        typedef unsigned int nop_t;
///addressing type for a combined object of c^dagger c
        typedef unsigned int spin_t;
///type of imaginary time values
        typedef double itime_t;
///addressing type of imaginary time indices (discretized)
        typedef unsigned int itime_index_t;
///addressing type of matsubara frequency
        typedef unsigned int frequency_t;

        typedef boost::uint_fast64_t my_uint64;
    }
}
