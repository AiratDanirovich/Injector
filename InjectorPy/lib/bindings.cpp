#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(InjectorGPN, m) {
    py::class_<Wrapper>(m, "Wrapper")
        .def(py::init<
            RealType,     // q
            RealType,     // density
            RealType,     // capacity
            RealType,     // viscosity
            RealType,     // rMin
            RealType,     // rMax
            size_t,       // r_nodes_nmbr
            RealType,     // zTop
            const VR&,    // thickness
            const VR&,    // conductivity
            const VR&,    // porosity
            const VR&,    // is_permeable
            const VR&,    // solid_density
            const VR&,    // solid_specific_heatcapacity
            RealType,     // t_start
            const VR&     // time_intervals
        >(),
        py::arg("q"),
        py::arg("density"),
        py::arg("capacity"),
        py::arg("viscosity"),
        py::arg("rMin"),
        py::arg("rMax"),
        py::arg("r_nodes_nmbr"),
        py::arg("zTop"),
        py::arg("thickness"),
        py::arg("conductivity"),
        py::arg("porosity"),
        py::arg("is_permeable"),
        py::arg("solid_density"),
        py::arg("solid_specific_heatcapacity"),
        py::arg("t_start"),
        py::arg("time_intervals")
        )
        .def("get_times", &Wrapper::get_times)
        .def("get_temps", &Wrapper::get_temps);
}
