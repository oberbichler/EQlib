#pragma once

#include "Define.h"
#include "Parameter.h"

namespace EQlib {

class Node {
private:    // variables
    Parameter m_x;
    Parameter m_y;
    Parameter m_z;

public:     // constructors
    Node(const double x, const double y, const double z)
    : m_x(x, 0)
    , m_y(y, 0)
    , m_z(z, 0)
    { }

    Node() : Node(0, 0, 0) { }

public:     // getters and setters
    Parameter& x() {
        return m_x;
    }

    Parameter& y() {
        return m_y;
    }

    Parameter& z() {
        return m_z;
    }

    Vector3D ref_location() const {
        return Vector3D(m_x.ref_value(), m_y.ref_value(), m_z.ref_value());
    }

    void set_ref_location(Vector3D value) {
        m_x.set_ref_value(value(0));
        m_y.set_ref_value(value(1));
        m_z.set_ref_value(value(2));
    }

    Vector3D act_location() const {
        return Vector3D(m_x.act_value(), m_y.act_value(), m_z.act_value());
    }

    void set_act_location(Vector3D value) {
        m_x.set_act_value(value(0));
        m_y.set_act_value(value(1));
        m_z.set_act_value(value(2));
    }

    Vector3D displacements() const {
        return Vector3D(m_x.delta(), m_y.delta(), m_z.delta());
    }

    void set_displacements(Vector3D value) {
        m_x.set_delta(value(0));
        m_y.set_delta(value(1));
        m_z.set_delta(value(2));
    }

    Vector3D forces() const {
        return Vector3D(m_x.target(), m_y.target(), m_z.target());
    }

    void set_forces(Vector3D value) {
        m_x.set_target(value(0));
        m_y.set_target(value(1));
        m_z.set_target(value(2));
    }
};

} // namespace EQlib