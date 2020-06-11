#pragma once

#ifdef EQLIB_USE_MKL
#include "PardisoLDLT.h"
#else
#include "SimplicialLDLT.h"
#endif

namespace eqlib {
#ifdef EQLIB_USE_MKL
    using SparseLDLT = PardisoLDLT;
#else
    using SparseLDLT = SimplicialLDLT;
#endif
} // namespace eqlib