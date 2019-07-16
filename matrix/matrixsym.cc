#include "matrix/matrixsym.h"

namespace hitnlls {
namespace matrix {

Matrixxi Matrixsym::OrderingByMinHeap() {
    Matrixxi ordering(Rows());
    for (int idx = 0; idx < Rows(); ++idx) {
        ordering[idx] = idx;
    }
    return ordering;
}

} // namespace matrix
} // namespace hitnlls