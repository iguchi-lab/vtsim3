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

    m.attr("SOLVE_LU")   = SOLVE_LU;
    m.attr("SOLVE_SOR")  = SOLVE_SOR;
    m.attr("STEP_P")     = STEP_P;          //偏微分時の圧力変化
    m.attr("VENT_ERR")   = VENT_ERR;        //換気回路網の許容残差
    m.atrr("STEP_T")     = STEP_T;          //偏微分時の温度変化
    m.atrr("THRM_ERR")   = THRM_ERR;        //熱回路網の許容残差
    m.atrr("CONV_ERR")   = CONV_ERR;        //収束許容誤差
    m.atrr("SOR_RATIO")  = SOR_RATIO;       //SOR法の緩和係数
    m.atrr("SOR_ERR")    = SOR_ERR;         //SOR法の許容残差
    m.atrr("SN_NONE")    = SN_NONE;         //計算しない
    m.atrr("SN_CALC")    = SN_CALC;         //計算する
    m.atrr("SN_FIX")     = SN_FIX;          //固定値（計算には利用するが、更新しない）
    m.atrr("SN_DLY")     = SN_DLY;          //遅延（熱容量計算用）
    m.atrr("VN_SIMPLE")  = VN_SIMPLE;       //換気回路網：単純開口
    m.atrr("VN_GAP")     = VN_GAP;          //換気回路網：隙間
    m.atrr("VN_FIX")     = VN_FIX;          //換気回路網：風量固定
    m.atrr("VN_AIRCON")  = VN_AIRCON;       //換気回路網：エアコン=風量固定、換気による熱移動=0
    m.atrr("VN_FAN")     = VN_FAN;          //換気回路網：送風ファン、PQ特性
    m.atrr("TN_SIMPLE")  = TN_SIMPLE;       //熱回路網：単純熱回路
    m.atrr("TN_AIRCON")  = TN_AIRCON;       //熱回路網：エアコン、熱量収支付け替え
    m.atrr("TN_SOLAR")   = TN_SOLAR;        //熱回路網：日射取得
    m.atrr("TN_GROUND")  = TN_GROUND;       //熱回路網：地盤
    m.atrr("TN_HEATER")  = TN_HEATER;       //熱回路網：発熱
    m.atrr("AC_AUTO")    = AC_AUTO;         //エアコン：自動
    m.atrr("AC_HEATING") = AC_HEATING;      //エアコン：暖房
    m.atrr("AC_COOLING") = AC_COOLING;      //エアコン：冷房
    m.atrr("AC_STOP")    = AC_STOP;         //エアコン：停止

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