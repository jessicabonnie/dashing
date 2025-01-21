#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../sketchcorehll.cpp"  // Adjust path as needed

namespace py = pybind11;

void bind_hll(py::module &m) {
    py::class_<sketch::hll>(m, "HLL")
        .def(py::init<uint8_t>())
        // Core operations
        .def("add", &sketch::hll::add)
        .def("estimate", &sketch::hll::estimate)
        // Additional utility methods
        .def("clear", &sketch::hll::clear)
        .def("merge", &sketch::hll::merge)
        
        // Properties
        .def_property_readonly("size", &sketch::hll::size)
        .def_property_readonly("precision", &sketch::hll::get_precision)
        
        // Optional: Add numpy array interface
        .def("get_registers", [](const sketch::hll &self) {
            auto data = self.get_core();
            auto size = self.size();
            return py::array_t<uint8_t>(
                {size},                // shape
                {sizeof(uint8_t)},     // stride
                data,                  // data pointer
                py::cast(self)         // owner object
            );
        })
        // Add k-mer processing
        .def("add_sequence", [](sketch::hll &self, const std::string &seq, size_t k) {
            for(size_t i = 0; i <= seq.length() - k; ++i) {
                uint64_t hash = sketch::FgHash()(seq.substr(i, k));
                self.add(hash);
            }
        });
} 