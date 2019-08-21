#pragma once

#include "Define.h"
#include "Parameter.h"

#include <memory>
#include <string>
#include <unordered_map>

namespace EQlib {

class Node
{
private:    // variables
    Pointer<Parameter> m_x;
    Pointer<Parameter> m_y;
    Pointer<Parameter> m_z;

    std::unordered_map<std::string, Parameter> m_parameters;

public:     // constructors
    Node(const double x, const double y, const double z) noexcept
    : m_x(std::make_shared<Parameter>(x))
    , m_y(std::make_shared<Parameter>(y))
    , m_z(std::make_shared<Parameter>(z))
    { }

    Node() noexcept : Node(0, 0, 0)
    { }

public:     // getters and setters
    Pointer<Parameter> x() noexcept
    {
        return m_x;
    }

    Pointer<Parameter> y() noexcept
    {
        return m_y;
    }

    Pointer<Parameter> z() noexcept
    {
        return m_z;
    }

    Vector3D ref_location() const noexcept
    {
        return Vector3D(m_x->ref_value(), m_y->ref_value(), m_z->ref_value());
    }

    void set_ref_location(Vector3D value) noexcept
    {
        m_x->set_ref_value(value[0]);
        m_y->set_ref_value(value[1]);
        m_z->set_ref_value(value[2]);
    }

    Vector3D act_location() const noexcept
    {
        return Vector3D(m_x->act_value(), m_y->act_value(), m_z->act_value());
    }

    void set_act_location(Vector3D value) noexcept
    {
        m_x->set_act_value(value[0]);
        m_y->set_act_value(value[1]);
        m_z->set_act_value(value[2]);
    }

    Vector3D displacements() const noexcept
    {
        return Vector3D(m_x->delta(), m_y->delta(), m_z->delta());
    }

    void set_displacements(Vector3D value) noexcept
    {
        m_x->set_delta(value[0]);
        m_y->set_delta(value[1]);
        m_z->set_delta(value[2]);
    }

    Vector3D forces() const noexcept
    {
        return Vector3D(m_x->target(), m_y->target(), m_z->target());
    }

    void set_forces(Vector3D value) noexcept
    {
        m_x->set_target(value[0]);
        m_y->set_target(value[1]);
        m_z->set_target(value[2]);
    }

public:     // operators
    Parameter& operator[](const std::string& name) noexcept
    {
        if (name == "x") {
            return *m_x;
        } else if (name == "y") {
            return *m_y;
        } else if (name == "z") {
            return *m_z;
        }

        return m_parameters[name];
    }

public:     // methods
    bool has_parameter(const std::string& name) noexcept
    {
        return m_parameters.find(name) != m_parameters.end();
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = EQlib::Node;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Holder>(m, "Node", py::dynamic_attr())
            // constructors
            .def(py::init<>())
            .def(py::init<double, double, double>(), "x"_a=0, "y"_a=0, "z"_a=0)
            // readonly properties
            .def_property_readonly("x", &Type::x)
            .def_property_readonly("y", &Type::y)
            .def_property_readonly("z", &Type::z)
            // properties
            .def_property("ref_location", &Type::ref_location,
                &Type::set_ref_location)
            .def_property("act_location", &Type::act_location,
                &Type::set_act_location)
            .def_property("displacements", &Type::displacements,
                &Type::set_displacements)
            .def_property("forces", &Type::forces, &Type::set_forces)
            // methods
            .def("has_parameter", &Type::has_parameter, "name"_a)
            // operators
            .def("__getitem__", &Type::operator[],
                py::return_value_policy::reference_internal)
        ;
    }
};

} // namespace EQlib