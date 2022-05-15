#include "../common.h"

#include <eqlib/objectives/spring.h>

void bind_spring(py::module_ &m)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace eqlib::objectives;
    using namespace pybind11::literals;

    auto cls = bind_objective<Spring>(m, "Spring");
    
    cls
        .def(py::init<Pointer<Node>>(), "node"_a);
}