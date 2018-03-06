#pragma once

#include <alps/accumulators.hpp>

namespace alps {
    namespace ctint {
        using SimpleRealVectorObservable = alps::accumulators::NoBinningAccumulator<std::vector<double> >;
        using SimpleRealObservable = alps::accumulators::NoBinningAccumulator<double>;
    }
}
