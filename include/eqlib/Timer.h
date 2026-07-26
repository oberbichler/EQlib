#pragma once

#include <chrono>

namespace eqlib {

class Timer {
private: // types
    using Type = Timer;
    using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using Duration = std::chrono::duration<double>;

private: // variables
    Time m_start;

private: // methods
    static Time now() noexcept
    {
        return std::chrono::high_resolution_clock::now();
    }

public: // constructors
    Timer() noexcept
        : m_start(now())
    {
    }

public: // methods
    void start() noexcept
    {
        m_start = now();
    }

    double elapsed() const noexcept
    {
        const Duration duration = now() - m_start;
        return duration.count();
    }
};

} // namespace eqlib
