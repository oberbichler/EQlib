#pragma once

#include "Define.h"

#include <pybind11/pybind11.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <stack>

namespace EQlib {

class Log
{
private:    // methods
    static std::shared_ptr<spdlog::logger>
    create()
    {
        auto logger = spdlog::stdout_color_mt("console");

        logger->set_pattern("%H:%M:%S.%e %^[%L]%$ %v");

        return logger;
    }

private:    // variables
    static const inline std::shared_ptr<spdlog::logger> s_console = create();

    std::stack<int> m_info_levels;

public:     // constructors
    Log(py::dict options = py::dict())
    {
        // options

        const auto info_level = get_or_default<int>(options, "info_level", 1);

        push_info_level(info_level);
    }

public:     // methods
    void push_info_level(const int value)
    {
        m_info_levels.push(value);
    }

    int pop_info_level()
    {
        const auto value = info_level();

        m_info_levels.pop();

        if (m_info_levels.empty()) {
            push_info_level(1);
        }

        return value;
    }

    int info_level()
    {
        if (!m_info_levels.empty()) {
            return m_info_levels.top();
        } else {
            return 1;
        }
    }

    template<class... TArgs>
    void debug(TArgs&&... args)
    {
        s_console->debug(std::forward<TArgs>(args)...);
    }

    template<class... TArgs>
    void info(TArgs&&... args)
    {
        s_console->info(std::forward<TArgs>(args)...);
    }

    template<class... TArgs>
    void info(const int level, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        info(std::forward<TArgs>(args)...);
    }

    template<class... TArgs>
    void error(TArgs&&... args)
    {
        s_console->error(std::forward<TArgs>(args)...);
    }

    template<class... TArgs>
    void warn(TArgs&&... args)
    {
        s_console->warn(std::forward<TArgs>(args)...);
    }

    template<class... TArgs>
    void critical(TArgs&&... args)
    {
        s_console->critical(std::forward<TArgs>(args)...);
    }
};

} // namespace EQlib