#pragma once

#include <alps/accumulators.hpp>

using SimpleRealVectorObservable = alps::accumulators::NoBinningAccumulator<std::vector<double> >;
using SimpleRealObservable = alps::accumulators::NoBinningAccumulator<double>;
