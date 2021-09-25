#pragma once

#include "../common.h"
#include "../node.h"
#include "../objective.h"
#include "../variable.h"

#include <vector>
#include <string>

namespace eqlib
{

    struct MembraneTriangle : public Objective
    {
        MembraneTriangle(const std::vector<Pointer<Node>> &nodes) : Objective()
        {
            const auto node_a = nodes[0];
            const auto node_b = nodes[1];
            const auto node_c = nodes[2];

            add_parameter(node_a->x());
            add_parameter(node_a->y());
            add_parameter(node_a->z());
            add_parameter(node_b->x());
            add_parameter(node_b->y());
            add_parameter(node_b->z());
            add_parameter(node_c->x());
            add_parameter(node_c->y());
            add_parameter(node_c->z());
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

} // namespace eqlib