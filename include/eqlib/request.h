#pragma once

namespace eqlib {

enum Request {
    F = 0b00000001,
    Df = 0b00000010,
    Hf = 0b00000100,
    G = 0b00010000,
    Dg = 0b00100000,
    Hg = 0b01000000,
    Hm = 0b01000100,
};

} // namespace eqlib