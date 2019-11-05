#pragma once

#include <chrono>

namespace eqlib {

class Timer
{
private:    // types
    using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using Duration = std::chrono::duration<double>;

private:    // variables
    Time m_start;

private:    // methods
    static Time now() noexcept
    {
        return std::chrono::high_resolution_clock::now();
    }

public:     // constructors
    Timer() noexcept
    : m_start(now())
    {
    }

public:     // methods
    void start() noexcept
    {
        m_start = now();
    }

    double ellapsed() const noexcept
    {
        const Duration duration = now() - m_start;
        return duration.count();
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace py::literals;

        using Type = Timer;

        py::class_<Type>(m, "Timer")
            .def(py::init<>())
            .def("start", &Type::start)
            .def_property_readonly("ellapsed", &Type::ellapsed)
        ;
    }
};

} // namespace eqlib