#pragma once

#include "kmat/runtime.hpp"

#include <algorithm>
#include <thread>
#include <vector>

namespace kmat {

template <typename Fn>
void parallel_for(std::size_t begin, std::size_t end, std::size_t num_threads, Fn&& fn) {
  if (begin >= end) {
    return;
  }
  const std::size_t n = end - begin;
  if (num_threads <= 1 || n == 1) {
    for (std::size_t i = begin; i < end; ++i) {
      fn(i);
    }
    return;
  }

  const std::size_t workers = std::min(num_threads, n);
  std::vector<std::thread> threads;
  threads.reserve(workers);

  const std::size_t chunk = (n + workers - 1) / workers;
  for (std::size_t w = 0; w < workers; ++w) {
    const std::size_t lo = begin + w * chunk;
    if (lo >= end) {
      break;
    }
    const std::size_t hi = std::min(end, lo + chunk);
    threads.emplace_back([lo, hi, &fn]() {
      for (std::size_t i = lo; i < hi; ++i) {
        fn(i);
      }
    });
  }
  for (std::thread& t : threads) {
    t.join();
  }
}

}  // namespace kmat
