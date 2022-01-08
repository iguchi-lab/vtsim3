#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <stdio.h>
#include <math.h>                   //数学　ヘッダーファイルの読み込み
#include <iostream>                 //出力　ヘッダーファイルの読み込み
#include <iomanip>
#include <vector>
#include <tuple>

#include "vtsimc.h"

using namespace::std;
namespace py = pybind11;

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, 
      vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> calc(InputData inp){
    VTSim calc(inp);
    calc.calc();
    return calc.result();
}

PYBIND11_MODULE(vtsimc, m) {
    m.def("calc", &calc, "");
    m.attr("SOLVE_LU")  = SOLVE_LU;
    m.attr("SOLVE_SOR") = SOLVE_SOR;
    py::class_<InputData>(m, "InputData")
        .def(py::init<>())
        .def_readwrite("sts",           &InputData::sts)
        .def_readwrite("length",        &InputData::length)
        .def_readwrite("t_step",        &InputData::t_step)
        .def_readwrite("nodes",         &InputData::nodes)
        .def_readwrite("v_nets",        &InputData::v_nets)
        .def_readwrite("t_nets",        &InputData::t_nets)
        .def_readwrite("sn_P_set",      &InputData::sn_P_set)
        .def_readwrite("sn_C_set",      &InputData::sn_C_set)
        .def_readwrite("sn_T_set",      &InputData::sn_T_set)
        .def_readwrite("sn_h_sr_set",   &InputData::sn_h_sr_set)
        .def_readwrite("sn_h_inp_set",  &InputData::sn_h_inp_set)
        .def_readwrite("sn_v_set",      &InputData::sn_v_set)
        .def_readwrite("sn_capa_set",   &InputData::sn_capa_set)
        .def_readwrite("sn_m_set",      &InputData::sn_m_set)
        .def_readwrite("sn_beta_set",   &InputData::sn_beta_set)
        .def_readwrite("vn_simple_set", &InputData::vn_simple_set)
        .def_readwrite("vn_gap_set",    &InputData::vn_gap_set)
        .def_readwrite("vn_fix_set",    &InputData::vn_fix_set)
        .def_readwrite("vn_fan_set",    &InputData::vn_fan_set)
        .def_readwrite("vn_eta_set",    &InputData::vn_eta_set)
        .def_readwrite("tn_simple_set", &InputData::tn_simple_set)
        .def_readwrite("tn_aircon_set", &InputData::tn_aircon_set)
        .def_readwrite("tn_solar_set",  &InputData::tn_solar_set)
        .def_readwrite("tn_ground_set", &InputData::tn_ground_set);
}