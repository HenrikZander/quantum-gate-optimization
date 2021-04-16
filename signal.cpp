#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j)
{
    return i + j;
}

PYBIND11_MODULE(signal, m)
{
    m.doc() = "Different signals that will be applied to the tunnable bus in the quantum circuit."; // optional module docstring

    m.def("add", &add, "A function which adds two numbers", py::arg("i"), py::arg("j"));
}