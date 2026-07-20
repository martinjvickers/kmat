#include "kmat/log.hpp"

#include <iostream>
#include <mutex>

namespace kmat {

namespace {

std::mutex g_log_mutex;
LogLevel g_log_level = LogLevel::Info;

void log_at(LogLevel level, std::string_view tag, std::string_view message) {
  if (level < g_log_level) {
    return;
  }
  std::lock_guard lock(g_log_mutex);
  std::cerr << "[kmat:" << tag << "] " << message << '\n';
}

}  // namespace

void set_log_level(LogLevel level) { g_log_level = level; }

LogLevel log_level() { return g_log_level; }

void log_debug(std::string_view message) { log_at(LogLevel::Debug, "debug", message); }

void log_info(std::string_view message) { log_at(LogLevel::Info, "info", message); }

void log_warn(std::string_view message) { log_at(LogLevel::Warn, "warn", message); }

void log_error(std::string_view message) { log_at(LogLevel::Error, "error", message); }

}  // namespace kmat
