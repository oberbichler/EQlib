#pragma once

#include "../node.h"
#include "../objective.h"

#include <string>
#include <vector>

namespace eqlib {
namespace objectives {

struct Spring : public Objective {
    Spring(const Pointer<Node>& node)
        : Objective(3)
    {
        m_parameters[0] = node->x();
        m_parameters[1] = node->y();
        m_parameters[2] = node->z();
    }

    // compute

    double compute(Ref<Vector> df, Ref<Matrix> hm, Request request)
    {
        df(0) = 1;
        df(1) = 2;
        df(2) = 3;

        return 1.0;
    }
};

} // namespace objectives
} // namespace eqlib