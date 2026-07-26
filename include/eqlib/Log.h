#pragma once

#include "Define.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace eqlib {

class Log {
private: // types
    using Type = eqlib::Log;

private: // methods
    static Pointer<spdlog::logger> create()
    {
        auto logger = spdlog::stdout_color_mt("console");

        logger->set_pattern("%H:%M:%S.%e  %v");

        return logger;
    }

private: // variables
    static const inline Pointer<spdlog::logger> s_console = create();

    static inline int s_info_level = 0;

public: // methods
    static int info_level()
    {
        return s_info_level;
    }

    static void set_info_level(const int value)
    {
        s_info_level = value;
    }

    template <class... TArgs>
    static void debug(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        s_console->debug(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void task_begin(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        info(1, "\033[1;32m> {}\033[0m", fmt::format(format, std::forward<TArgs>(args)...));
    }

    template <class... TArgs>
    static void task_end(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        info(1, "{}", fmt::format(format, std::forward<TArgs>(args)...));
    }

    template <class... TArgs>
    static void task_info(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        info(2, "\033[37m{}\033[0m", fmt::format(format, std::forward<TArgs>(args)...));
    }

    template <class... TArgs>
    static void task_step(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        info(3, "\033[33m{}\033[0m", fmt::format(format, std::forward<TArgs>(args)...));
    }

    template <class... TArgs>
    static void info(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        s_console->info(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void info(const int level, fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        s_console->info(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void error(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        s_console->error(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void error(const int level, fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        s_console->error(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void warn(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        s_console->warn("\033[35m{}\033[0m", fmt::format(format, std::forward<TArgs>(args)...));
    }

    template <class... TArgs>
    static void warn(const int level, fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        warn(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void critical(fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        s_console->critical(format, std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void critical(const int level, fmt::format_string<TArgs...> format, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        s_console->critical(format, std::forward<TArgs>(args)...);
    }
};

} // namespace eqlib
