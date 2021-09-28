#pragma once

#include "../node.h"
#include "../objective.h"

#include <string>
#include <vector>

namespace eqlib {
namespace objectives {

struct Spring : public Objective {
    Spring(const Pointer<Node>& node)
        : Objective()
    {
        add_variable(node->x());
        add_variable(node->y());
        add_variable(node->z());
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