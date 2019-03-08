#pragma once

#include <chrono>
#include <stack>

namespace EQlib {

class Timer
{
private:
    using TimePoint = std::chrono::time_point<
        std::chrono::high_resolution_clock>;
    using Duration = std::chrono::duration<double>;

private:
    TimePoint m_start;

public:
    Timer()
    : m_start(std::chrono::high_resolution_clock::now())
    { }

    void start()
    {
        m_start = std::chrono::high_resolution_clock::now();
    }

    double ellapsed() const
    {
        const auto now = std::chrono::high_resolution_clock::now();
        const Duration duration = now - m_start;
        return duration.count();
    }
};

} // namespace EQlib