#pragma once

#include <string>
#include <unordered_map>
#include <variant>

namespace eqlib {

using Settings = std::unordered_map<std::string, std::variant<int, double, std::string>>;

std::string get_or_default(const Settings& options, const std::string key, const std::string default_value)
{
    const auto& it = options.find(key);

    if (it == options.end()) {
        return default_value;
    }

    return std::get<std::string>(it->second);
}

int get_or_default(const Settings& options, const std::string key, const int default_value)
{
    const auto& it = options.find(key);

    if (it == options.end()) {
        return default_value;
    }

    return std::get<int>(it->second);
}

double get_or_default(const Settings& options, const std::string key, const double default_value)
{
    const auto& it = options.find(key);

    if (it == options.end()) {
        return default_value;
    }

    if (std::holds_alternative<int>(it->second)) {
        return std::get<int>(it->second);
    }

    return std::get<double>(it->second);
}

} // namespace eqlib