#pragma once

#include <alps/accumulators.hpp>

namespace alps {
    namespace ctint {
        typedef alps::accumulators::NoBinningAccumulator<double> SimpleRealObservable;
        typedef alps::accumulators::NoBinningAccumulator<std::vector<double> > SimpleRealVectorObservable;
    }
}
