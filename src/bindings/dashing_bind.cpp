#include <pybind11/pybind11.h>
#include "hll_bind.h"

PYBIND11_MODULE(dashing, m) {
    m.doc() = "Python bindings for Dashing's HLL sketching"; 
    bind_hll(m);
} 