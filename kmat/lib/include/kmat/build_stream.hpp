#pragma once

#include "kmat/error.hpp"
#include "kmat/matrix.hpp"

namespace kmat {

/// Memory-bounded v2 build from sorted `.kset` files (multiway merge + hash shards).
Error build_matrix_from_presence_sets_streaming(const BuildOptions& opts);

}  // namespace kmat
