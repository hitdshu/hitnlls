#pragma once

#include "matrix/matrixx.h"

namespace hitnlls {
namespace matrix {

class Matrixsym : public Matrixx<unsigned char> {
public:
    Matrixxi OrderingByMinHeap();
};

} // namespace matrix
} // namespace hitnlls