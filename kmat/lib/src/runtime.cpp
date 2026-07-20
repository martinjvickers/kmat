#include "kmat/runtime.hpp"

#include <algorithm>
#include <cctype>
#include <mutex>
#include <thread>

namespace kmat {

namespace {

std::mutex g_runtime_mu;
RuntimeConfig g_runtime{};
bool g_runtime_inited = false;

void ensure_runtime_inited() {
  if (!g_runtime_inited) {
    g_runtime = resolve_runtime(RuntimeProfile::Laptop, 0);
    g_runtime_inited = true;
  }
}

std::string lower_copy(std::string s) {
  for (char& c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

}  // namespace

RuntimeConfig resolve_runtime(RuntimeProfile profile, std::size_t num_threads_override) {
  RuntimeConfig cfg;
  cfg.profile = profile;
  const unsigned hw = std::max(1u, std::thread::hardware_concurrency());

  if (profile == RuntimeProfile::Hpc) {
    cfg.num_threads = hw;
    cfg.io_buffer_bytes = 8u << 20;  // 8 MiB
  } else {
    // Laptop: keep headroom for interactive use / shared machines.
    cfg.num_threads = std::max(1u, hw / 2);
    cfg.io_buffer_bytes = 1u << 20;  // 1 MiB
  }

  if (num_threads_override > 0) {
    cfg.num_threads = num_threads_override;
  }
  return cfg;
}

void set_runtime_config(const RuntimeConfig& cfg) {
  std::lock_guard<std::mutex> lock(g_runtime_mu);
  g_runtime = cfg;
  if (g_runtime.num_threads == 0) {
    g_runtime = resolve_runtime(g_runtime.profile, 0);
  }
  g_runtime_inited = true;
}

RuntimeConfig runtime_config() {
  std::lock_guard<std::mutex> lock(g_runtime_mu);
  ensure_runtime_inited();
  return g_runtime;
}

std::size_t effective_threads(const RuntimeConfig& cfg) {
  if (cfg.num_threads > 0) {
    return cfg.num_threads;
  }
  return resolve_runtime(cfg.profile, 0).num_threads;
}

RuntimeProfile parse_runtime_profile(const std::string& name, bool& ok) {
  const std::string n = lower_copy(name);
  ok = true;
  if (n == "laptop" || n == "local") {
    return RuntimeProfile::Laptop;
  }
  if (n == "hpc" || n == "cluster") {
    return RuntimeProfile::Hpc;
  }
  ok = false;
  return RuntimeProfile::Laptop;
}

const char* runtime_profile_name(RuntimeProfile profile) {
  switch (profile) {
    case RuntimeProfile::Hpc:
      return "hpc";
    case RuntimeProfile::Laptop:
    default:
      return "laptop";
  }
}

}  // namespace kmat
