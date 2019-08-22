#pragma once

#include <chrono>
#include <stack>

namespace EQlib {

class Timer
{
private:
    using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using Duration = std::chrono::duration<double>;

private:
    Time m_start;

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

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = EQlib::Timer;

        py::class_<Type>(m, "Timer")
            .def(py::init<>())
            .def("start", &Type::start)
            .def_property_readonly("ellapsed", &Type::ellapsed)
        ;
    }
};

} // namespace EQlib