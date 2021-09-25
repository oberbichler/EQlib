#pragma once

#include <Eigen/Core>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <cstddef>
#include <limits>
#include <memory>

#include <mkl.h>

#define EQLIB_VERSION_MAJOR 1
#define EQLIB_VERSION_MINOR 0
#define EQLIB_VERSION_PATCH 0.dev1

#define EQLIB_STRINGIFY(x) #x
#define EQLIB_TOSTRING(x)  EQLIB_STRINGIFY(x)
#define EQLIB_VERSION                                                          \
    (EQLIB_TOSTRING(EQLIB_VERSION_MAJOR) "."                                   \
     EQLIB_TOSTRING(EQLIB_VERSION_MINOR) "."                                   \
     EQLIB_TOSTRING(EQLIB_VERSION_PATCH))

namespace eqlib {

const double LOWEST = std::numeric_limits<float>::lowest();
const double HIGHEST = std::numeric_limits<float>::max();

using index = std::ptrdiff_t;

template <typename T>
inline index len(const T& container) noexcept
{
    return static_cast<index>(container.size());
}

template <typename T>
using Unique = std::unique_ptr<T>;

template <typename T>
using Pointer = std::shared_ptr<T>;

template <typename T>
using Weak = std::weak_ptr<T>;

template <typename T, typename... TArgs>
Unique<T> new_(TArgs&&... args)
{
    return Unique<T>(new T(std::forward<TArgs>(args)...));
}

template <typename T>
using RobinSet = tsl::robin_set<T>;

template <typename K, typename T>
using RobinMap = tsl::robin_map<K, T>;

using Vector3 = Eigen::Matrix<double, 1, 3, Eigen::RowMajor>;

using Vector = Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using Ref = Eigen::Ref<T>;

template <typename T>
using Map = Eigen::Map<T, Eigen::RowMajor>;

struct Logger {
    static Pointer<spdlog::logger> create()
    {
        auto logger = spdlog::stdout_color_mt("console");

        logger->set_pattern("%H:%M:%S.%e  %v");

        return logger;
    }

    static const inline Pointer<spdlog::logger> s_console = create();

    static inline int s_info_level = 0;

    static int info_level()
    {
        return s_info_level;
    }

    static void set_info_level(const int value)
    {
        s_info_level = value;
    }

    template <class... TArgs>
    static void debug(TArgs&&... args)
    {
        s_console->debug(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void task_begin(std::string text, TArgs&&... args)
    {
        std::string message = "\u001b[1;32m> " + text + "\u001b[0m";
        info(1, message.c_str(), std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void task_end(std::string text, TArgs&&... args)
    {
        info(1, text.c_str(), std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void task_info(std::string text, TArgs&&... args)
    {
        std::string message = "\033[37m" + text + "\033[0m";
        info(2, message.c_str(), std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void task_step(std::string text, TArgs&&... args)
    {
        std::string message = "\033[33m" + text + "\033[0m";
        info(3, message.c_str(), std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void info(TArgs&&... args)
    {
        s_console->info(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void info(const int level, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        info(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void error(TArgs&&... args)
    {
        s_console->error(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void error(const int level, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        s_console->error(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void warn(std::string text, TArgs&&... args)
    {
        std::string message = "\033[35m" + text + "\033[0m";
        s_console->warn(message.c_str(), std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void warn(const int level, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        warn(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void critical(TArgs&&... args)
    {
        s_console->critical(std::forward<TArgs>(args)...);
    }

    template <class... TArgs>
    static void critical(const int level, TArgs&&... args)
    {
        if (level > info_level()) {
            return;
        }

        s_console->critical(std::forward<TArgs>(args)...);
    }
};

enum Request {
    F = 1 << 0,
    G = 1 << 1,
    Df = 1 << 2,
    Dg = 1 << 3,
    Hf = 1 << 4,
    Hg = 1 << 5
};

} // namespace eqlib