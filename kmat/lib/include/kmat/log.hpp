#pragma once

#include <string>
#include <string_view>

namespace kmat {

enum class LogLevel { Debug, Info, Warn, Error };

void set_log_level(LogLevel level);
LogLevel log_level();

void log_debug(std::string_view message);
void log_info(std::string_view message);
void log_warn(std::string_view message);
void log_error(std::string_view message);

}  // namespace kmat
