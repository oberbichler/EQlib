#pragma once

namespace eqlib {

enum Request {
    F   = 0b00000001,
    DF  = 0b00000010,
    HF  = 0b00000100,
    G   = 0b00010000,
    DG  = 0b00100000,
    HG  = 0b01000000,
    HM  = 0b01000100,
    R   = 0b10000000,
    All = 0b11111111,
};

} // namespace eqlib