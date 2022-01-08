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
    VTSim calc;
    calc.set_inp(inp);
    calc.calc();
    return calc.result();
}

PYBIND11_MODULE(vtsimc, m) {
    m.def("calc", &calc, "");

    m.attr("SN_NONE")    = SN_NONE;
    m.attr("SN_CALC")    = SN_CALC;
    m.attr("SN_FIX")     = SN_FIX;
    m.attr("SN_DLY")     = SN_DLY;
    
    m.attr("VN_SIMPLE")  = VN_SIMPLE;
    m.attr("VN_GAP")     = VN_GAP;
    m.attr("VN_FIX")     = VN_FIX;
    m.attr("VN_AIRCON")  = VN_AIRCON;
    m.attr("VN_FAN")     = VN_FAN;
    
    m.attr("TN_SIMPLE")  = TN_SIMPLE;
    m.attr("TN_AIRCON")  = TN_AIRCON;
    m.attr("TN_SOLAR")   = TN_SOLAR;
    m.attr("TN_GROUND")  = TN_GROUND;
    m.attr("TN_HEATER")  = TN_HEATER;
    
    m.attr("AC_AUTO")    = AC_AUTO;
    m.attr("AC_HEATING") = AC_HEATING;
    m.attr("AC_COOLING") = AC_COOLING;
    m.attr("AC_STOP")    = AC_STOP;

    py::class_<CalcStatus>(m, "CalcStatus")
        .def(py::init<>())
        .def_read_write("length",    &CalcStatus::length)
        .def_read_write("t_step",    &CalcStatus::t_step)
        .def_read_write("solve",     &CalcStatus::solve)
        .def_read_write("step_p",    &CalcStatus::step_p) 
        .def_read_write("vent_err",  &CalcStatus::vent_err) 
        .def_read_write("step_t",    &CalcStatus::step_t) 
        .def_read_write("thrm_err",  &CalcStatus::thrm_err)
        .def_read_write("conv_err",  &CalcStatus::conv_err) 
        .def_read_write("sor_ratio", &CalcStatus::sor_ratio) 
        .def_read_write("sor_err",   &CalcStatus::sor_err);
    
    py::class_<InputData>(m, "InputData")
        .def(py::init<>())
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

    py::class_<VTSim>(m, "VTSim")
        .def(py::init<>())
        .def("set_calc_status", &VTSim::set_calc_status, "")
        .def("set_inp",         &VTSim::set_inp, "")
        .def("calc",            &VTSim::calc,    "")
        .def("result",          &VTSim::result,  "")
        .def_readwrite("sn",    &VTSim::sn, "")
        .def_readwrite("vn",    &VTSim::vn, "")
        .def_readwrite("tn",    &VTSim::tn, "");
}