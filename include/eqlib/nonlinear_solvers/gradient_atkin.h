#pragma once

#include "../common.h"

namespace eqlib {

struct Momentum
{
    Momentum(const double alpha, const double beta, const Ref<Vector> x)
    {
        this->alpha = alpha;
        this->beta = beta;
        this->v = Vector::Zeros(len(x));
        this->g = Vector(len(g));
    }

    double alpha;
    double beta;
    Vector v;
    Vector g;

    void step(Vector(Vector) df, Ref<Vector> x)
    {
        g = df(x);
        v = beta * v - alpha * g;
        x += v;
    }
}

struct GradientAtkin
{

}

};

} // namespace eqlib