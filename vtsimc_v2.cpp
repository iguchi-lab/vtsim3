#include <stdio.h>
#include <math.h>                   //数学　ヘッダーファイルの読み込み
#include <iostream>                 //出力　ヘッダーファイルの読み込み
#include <iomanip>
#include <vector>
#include <tuple>

#include "vtsimc.h"

using namespace std;

void calc_01(void){
    InputData inp;
    inp.length        = 2;
    vector<double> vol(inp.length, 100.0 / 3600.0),
                   alpha(inp.length, 0.6),
                   area(inp.length, 1.0),
                   h(inp.length, 0.0);

    inp.sts           = {SOLVE_LU, STEP_P, VENT_ERR, STEP_T, THRM_ERR, CONV_ERR, SOR_RATIO, SOR_ERR};
    inp.nodes         = {{SN_FIX, SN_NONE, SN_NONE},
                         {SN_CALC, SN_NONE, SN_NONE},
                         {SN_CALC, SN_NONE, SN_NONE}};
    inp.v_nets        = {{0, 1, VN_FIX,    h, h}, 
                         {1, 2, VN_SIMPLE, h, h}, 
                         {2, 0, VN_SIMPLE, h, h}};                                                          
    inp.vn_simple_set = {{1, alpha, area},
                         {2, alpha, area}}; 
    inp.vn_fix_set    = {{0, vol}};

    VTSim calc(inp);
    calc.calc();
}

void calc_02(void){
    InputData inp;
    inp.length = 2;
    vector<double> t1(inp.length, 20.0),              
                   t2(inp.length,  0.0),
                   vol1(inp.length, 1000.0 / 3600.0), 
                   vol2(inp.length,  500.0 / 3600.0),
                   h_sr(inp.length, 4000.0),             
                   h_inp(inp.length, 2000.0),
                   cdtc1(inp.length, 100.0 * 0.87),
                   cdtc2(inp.length, 200.0 * 0.87),
                   ms(inp.length, 1.0),
                   h(inp.length, 0.0);
    vector<int>    ac_mode(inp.length, AC_COOLING);

    inp.sts           = {SOLVE_LU, STEP_P, VENT_ERR, STEP_T, THRM_ERR, CONV_ERR, SOR_RATIO, SOR_ERR};
    inp.nodes         = {{SN_NONE, SN_NONE, SN_CALC},                                                     //0: AC-out
                         {SN_NONE, SN_NONE, SN_CALC},                                                     //1: Room1
                         {SN_NONE, SN_NONE, SN_CALC},                                                     //2: Room2
                         {SN_NONE, SN_NONE, SN_FIX},                                                      //3: AC-in
                         {SN_NONE, SN_NONE, SN_FIX},                                                      //4: Outside
                         {SN_NONE, SN_NONE, SN_NONE},                                                     //5: solar
                         {SN_NONE, SN_NONE, SN_NONE}};                                                    //6: heater  

    inp.v_nets        = {{0, 1, VN_FIX,    h, h}, 
                         {0, 2, VN_FIX,    h, h},
                         {1, 3, VN_FIX,    h, h}, 
                         {2, 3, VN_FIX,    h, h},
                         {3, 0, VN_FIX,    h, h}};
    inp.t_nets        = {{1, 4, TN_SIMPLE}, 
                         {2, 4, TN_SIMPLE}, 
                         {3, 0, TN_AIRCON},
                         {2, 5, TN_SOLAR},
                         {2, 6, TN_HEATER}};
    inp.sn_T_set      = {{3, t1},
                         {4, t2}};
    inp.sn_h_sr_set   = {{5, h_sr}};
    inp.sn_h_inp_set  = {{6, h_inp}};
    inp.vn_fix_set    = {{0, vol2},
                         {1, vol2},
                         {2, vol2},
                         {3, vol2},
                         {4, vol1}};
    inp.tn_simple_set = {{0, cdtc1},
                         {1, cdtc2}};
    inp.tn_aircon_set = {{2, ac_mode, t1}};
    inp.tn_solar_set  = {{3, ms}};

    VTSim calc(inp);
    calc.calc();
}

void calc_03(void){
    InputData inp;
    inp.length        = 3;
    vector<double> alpha(inp.length, 0.6),
                   area(inp.length, 0.21),
                   a1(inp.length, 0.005516666666666667),
                   a2(inp.length, 0.018633333333333335),
                   a3(inp.length, 0.005516666666666667),
                   a4(inp.length, 0.006208333333333334),
                   a5(inp.length, 0.004141666666666667),
                   n(inp.length, 1.5),
                   cdtc01(inp.length, 15.292200000000001),
                   cdtc02(inp.length, 51.6516),
                   cdtc03(inp.length, 15.292200000000001),
                   cdtc04(inp.length, 17.209500000000002),
                   cdtc05(inp.length, 11.4807),
                   cdtc06(inp.length, 1.7745),
                   cdtc07(inp.length, 3.549),
                   cdtc08(inp.length, 2.4843),
                   cdtc09(inp.length, 2.4843),
                   cdtc10(inp.length, 1.2441),
                   h(inp.length, 0.0),
                   g_area(inp.length, 6.62),
                   rg(inp.length, 1.85);

    inp.sts           = {SOLVE_SOR, STEP_P, VENT_ERR, STEP_T, THRM_ERR, CONV_ERR, SOR_RATIO, 0.01};
    inp.nodes         = {{2, 0, 2}, {2, 0, 2}, {2, 0, 2}, 
                         {1, 0, 1}, {1, 0, 1}, {1, 0, 1}, {1, 0, 1}, {1, 0, 1}, 
                         {0, 0, 2}, {0, 0, 2}};

    inp.sn_T_set      = {{0, {21.55, 21.5,  21.4}}, 
                         {1, {21.4 , 21.45, 21.4}}, 
                         {2, {20.22, 20.12, 19.98}}, 
                         {8, {14.9,  14.6,  14.3}}, 
                         {9, {17.03, 17.03, 17.03}}};
    inp.v_nets        = {{3, 4, 0, h, h}, 
                         {4, 5, 0, h, h},
                         {4, 6, 0, h, h},
                         {6, 7, 0, h, h},
                         {3, 0, 1, h, h},
                         {4, 0, 1, h, h},
                         {5, 0, 1, h, h},
                         {6, 0, 1, h, h},
                         {7, 0, 1, h, h},
                         {0, 1, 2, h, h},
                         {1, 3, 2, h, h},
                         {0, 2, 2, h, h},
                         {2, 3, 2, h, h},
                         {4, 0, 4, h, h},
                         {5, 0, 4, h, h},
                         {6, 0, 4, h, h},
                         {6, 0, 4, h, h},
                         {7, 0, 4, h, h}};
    inp.t_nets        = {{3, 0, 0}, {4, 0, 0}, {5, 0, 0}, {6, 0, 0}, {7, 0, 0}, 
                         {3, 8, 0}, {4, 8, 0}, {5, 8, 0}, {6, 8, 0}, {7, 8, 0}, 
                         {3, 9, 3}, {4, 9, 3}, {5, 9, 3}, {6, 9, 3}, {7, 9, 3}};
    inp.vn_simple_set = {{0, alpha, area}, 
                         {1, alpha, area}, 
                         {2, alpha, area}, 
                         {3, alpha, area}};
    inp.vn_gap_set    = {{4, a1, n}, 
                         {5, a2, n}, 
                         {6, a3, n}, 
                         {7, a4, n}, 
                         {8, a5, n}};
    inp.vn_fix_set    = {{9,  {0.15178725, 0.15178725, 0.15178725}}, 
                         {10, {0.15178725, 0.15178725, 0.15178725}}, 
                         {11, {0.05555556, 0.05555556, 0.05555556}}, 
                         {12, {0.05555556, 0.05555556, 0.05555556}}};
    inp.vn_fan_set    = {{13, {0.04166667, 0.04166667, 0.04166667}, {80., 80., 80.}, {0., 0., 0.}, {80.,  80.,  80.}}, 
                         {14, {0.04166667, 0.04166667, 0.04166667}, {80., 80., 80.}, {0., 0., 0.}, {80.,  80.,  80.}}, 
                         {15, {0.04166667, 0.04166667, 0.04166667}, {80., 80., 80.}, {0., 0., 0.}, {80.,  80.,  80.}}, 
                         {16, {0.04166667, 0.04166667, 0.04166667}, {80., 80., 80.}, {0., 0., 0.}, {80.,  80.,  80.}}, 
                         {17, {0.04166667, 0.04166667, 0.04166667}, {80., 80., 80.}, {0., 0., 0.}, {80.,  80.,  80.}}};
    inp.tn_simple_set = {{0, cdtc01}, 
                         {1, cdtc02}, 
                         {2, cdtc03}, 
                         {3, cdtc04}, 
                         {4, cdtc05}, 
                         {5, cdtc06}, 
                         {6, cdtc07}, 
                         {7, cdtc08}, 
                         {8, cdtc09}, 
                         {9, cdtc10}};
    inp.tn_ground_set = {{10, g_area, rg, 0.0255, 
                         {0.99999618, 0.99998474, 0.99993896, 0.99975593, 0.99902388, 0.99610162, 0.98449151, 0.93941063, 0.77880078, 0.36787944}, 
                         {-3.7393300e-08,  7.3357800e-07, -9.9402100e-06,  5.8429900e-04, 3.8443500e-04,  
                           6.1867900e-04,  2.7429050e-03,  8.7726000e-04,  1.1610039e-02,  1.5655690e-03}},
                         {11, g_area, rg, 0.0255, 
                         {0.99999618, 0.99998474, 0.99993896, 0.99975593, 0.99902388, 0.99610162, 0.98449151, 0.93941063, 0.77880078, 0.36787944}, 
                         {-3.7393300e-08,  7.3357800e-07, -9.9402100e-06,  5.8429900e-04, 3.8443500e-04,  
                           6.1867900e-04,  2.7429050e-03,  8.7726000e-04,  1.1610039e-02,  1.5655690e-03}},
                         {12, g_area, rg, 0.0255, 
                         {0.99999618, 0.99998474, 0.99993896, 0.99975593, 0.99902388, 0.99610162, 0.98449151, 0.93941063, 0.77880078, 0.36787944}, 
                         {-3.7393300e-08,  7.3357800e-07, -9.9402100e-06,  5.8429900e-04, 3.8443500e-04,  
                           6.1867900e-04,  2.7429050e-03,  8.7726000e-04,  1.1610039e-02,  1.5655690e-03}},
                         {13, g_area, rg, 0.0255, 
                         {0.99999618, 0.99998474, 0.99993896, 0.99975593, 0.99902388, 0.99610162, 0.98449151, 0.93941063, 0.77880078, 0.36787944}, 
                         {-3.7393300e-08,  7.3357800e-07, -9.9402100e-06,  5.8429900e-04, 3.8443500e-04,  
                           6.1867900e-04,  2.7429050e-03,  8.7726000e-04,  1.1610039e-02,  1.5655690e-03}},
                         {14, g_area, rg, 0.0255, 
                         {0.99999618, 0.99998474, 0.99993896, 0.99975593, 0.99902388, 0.99610162, 0.98449151, 0.93941063, 0.77880078, 0.36787944}, 
                         {-3.7393300e-08,  7.3357800e-07, -9.9402100e-06,  5.8429900e-04, 3.8443500e-04,  
                           6.1867900e-04,  2.7429050e-03,  8.7726000e-04,  1.1610039e-02,  1.5655690e-03}}};
    VTSim calc(inp);
    calc.calc();
}

int get_xyz(int x, int y, int z){
    return x * 15 + y * 5 + z;
}

void calc_04(void){
    InputData inp;
    int X = 3, Y = 3, Z = 5;
    double cx[3]   = {1.52, 1.52, 1.52}, cy[3] = {1.02, 1.60, 1.02}, cz[5] = {0.48, 0.48, 0.48, 0.48, 0.48};
    double hcz[5]  = {0.24, 0.72, 1.20, 1.68, 2.16};
    double htz[5]  = {0.48, 0.96, 1.44, 1.92};
    double cdtc[6] = {0.4, 0.5, 0.94, 1.88, 3.49, 5.59};

    inp.length = 4;
    vector<double> Ti(inp.length, 20.0), 
                   To(inp.length, 6.0), 
                   H_inp(inp.length, 462.3),
                   alpha(inp.length, 1.0);

    inp.t_step = 0.1;

    inp.sts = {SOLVE_SOR, STEP_P, 1e-3, STEP_T, THRM_ERR, CONV_ERR, 0.8, 1.0e-6};

    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int z = 0; z < Z; z++){
                inp.nodes.push_back({SN_CALC, SN_NONE, SN_CALC});
                inp.sn_T_set.push_back({get_xyz(x, y, z), Ti});
            }
        }
    }


    inp.nodes.push_back({SN_NONE, SN_NONE, SN_FIX});
    inp.sn_T_set.push_back({45, To});
    inp.nodes.push_back({SN_NONE, SN_NONE, SN_FIX});
    inp.sn_T_set.push_back({46, Ti});
    inp.nodes.push_back({SN_NONE, SN_NONE, SN_NONE});
    inp.sn_h_inp_set.push_back({47, H_inp});

    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int z = 0; z < Z; z++){
                inp.nodes.push_back({SN_NONE, SN_NONE, SN_DLY});
                inp.sn_T_set.push_back({get_xyz(x, y, z) + 48, Ti});
                inp.sn_capa_set.push_back({get_xyz(x, y, z) + 48, get_xyz(x, y, z)});
            }
        }
    }

    int v_n = 0, t_n = 0;
    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int z = 0; z < Z; z++){
                if(x != X - 1){
                    vector<double> h(inp.length, hcz[z]);
                    inp.v_nets.push_back({get_xyz(x, y, z), get_xyz(x + 1, y, z), VN_SIMPLE, h, h});
                    double area = cy[y] * cz[z];
                    inp.vn_simple_set.push_back({v_n, alpha, {area, area, area, area}});
                    v_n++;
                }
                if(y != Y - 1){
                    vector<double> h(inp.length, hcz[z]);   
                    inp.v_nets.push_back({get_xyz(x, y, z), get_xyz(x, y + 1, z), VN_SIMPLE, h, h});
                    double area = cx[x] * cz[z];
                    inp.vn_simple_set.push_back({v_n, alpha, {area, area, area, area}});
                    v_n++;
                }
                if(z != Z - 1){
                    vector<double> h(inp.length, htz[z]);   
                    inp.v_nets.push_back({get_xyz(x, y, z), get_xyz(x, y, z + 1), VN_SIMPLE, h, h});
                    double area = cx[x] * cy[y];
                    inp.vn_simple_set.push_back({v_n, alpha, {area, area, area, area}});
                    v_n++;
                }
                if(z == 0 || z == Z - 1){  
                    inp.t_nets.push_back({get_xyz(x, y, z), 45, TN_SIMPLE});
                    double c = cdtc[0];
                    inp.tn_simple_set.push_back({t_n, {c, c, c, c}});
                    t_n++;
                }
                if(y == 0 || y == Y - 1){   
                    inp.t_nets.push_back({get_xyz(x, y, z), 45, TN_SIMPLE});
                    double c = cdtc[0];
                    inp.tn_simple_set.push_back({t_n, {0, 0, 0, 0}});
                    t_n++;
                }
                if(x == 0){
                    if(y == 1 && z < 4){
                        inp.t_nets.push_back({get_xyz(x, y, z), 54, TN_SIMPLE});
                        double c = cdtc[5];
                        inp.tn_simple_set.push_back({t_n, {c, c, c, c}});
                        t_n++;
                    }
                    else{
                        inp.t_nets.push_back({get_xyz(x, y, z), 45, TN_SIMPLE});
                        double c = cdtc[0];
                        inp.tn_simple_set.push_back({t_n, {c, c, c, c}});
                        t_n++;
                    }
                }
                if(x == X - 1){   
                    inp.t_nets.push_back({get_xyz(x, y, z), 46, TN_SIMPLE});
                    double c = cdtc[0];
                    inp.tn_simple_set.push_back({t_n, {c, c, c, c}});
                    t_n++;
                }
                inp.t_nets.push_back({get_xyz(x, y, z), get_xyz(x, y, z) + 48, TN_SIMPLE});
                double c = cx[x] * cy[y] * cz[z] * 1.205 * 1006 / inp.t_step;
                inp.tn_simple_set.push_back({t_n, {c, c, c, c}});
                t_n++;                
            }
        }
    }
    inp.t_nets.push_back({get_xyz(2, 1, 0), 47, TN_HEATER});

    VTSim calc(inp);
    calc.calc();
}

void calc_05(void){
    InputData inp;
    inp.length        = 3;

    vector<double>    c0(inp.length, 40000.0), 
                       m(inp.length, 28470.0),
                      v1(inp.length,  13.66),
                      v2(inp.length,   6.00),
                      v3(inp.length,  40.25),
                      v4(inp.length, 126.49),
                      v5(inp.length,  16.56),
                      v6(inp.length,   9.94),
                      v7(inp.length, 100.46),
                   beta0(inp.length, 0.00005),
                   vol00(inp.length,  320.0 / 3600),
                   vol01(inp.length, 1206.0 / 3600),
                   vol02(inp.length,  155.0 / 3600),
                   vol03(inp.length,  313.0 / 3600),
                   vol04(inp.length,  108.0 / 3600),
                   vol05(inp.length,  176.0 / 3600),
                   vol06(inp.length,  454.0 / 3600),
                   vol07(inp.length,  155.0 / 3600),
                   vol08(inp.length,  313.0 / 3600),
                   vol09(inp.length,  108.0 / 3600),
                   vol10(inp.length,  100.0 / 3600),
                   vol11(inp.length,   76.0 / 3600),
                   vol12(inp.length,  220.0 / 3600),
                   vol13(inp.length,  886.0 / 3600),              
                    eta0(inp.length,    1.0),
                    eta1(inp.length,    0.9),
                       h(inp.length,    0.0);

    inp.t_step        = 10;
    inp.sts           = {SOLVE_LU, STEP_P, VENT_ERR, STEP_T, THRM_ERR, CONV_ERR, SOR_RATIO, SOR_ERR};
    inp.nodes         = {{2, 2, 0}, {2, 1, 0}, {2, 1, 0}, {2, 1, 0}, {2, 1, 0}, {2, 1, 0}, {2, 1, 0}, {2, 1, 0}};   
    inp.v_nets        = {{0, 1, 2, h, h},
                         {1, 2, 2, h, h},
                         {2, 3, 2, h, h},
                         {2, 4, 2, h, h},
                         {2, 5, 2, h, h},
                         {2, 6, 2, h, h},
                         {2, 7, 2, h, h},
                         {3, 7, 2, h, h},
                         {4, 7, 2, h, h},
                         {5, 7, 2, h, h},                     
                         {6, 0, 2, h, h},
                         {6, 7, 2, h, h},
                         {7, 0, 2, h, h},
                         {7, 2, 2, h, h}};
    inp.sn_C_set      = {{0, c0},
                         {1, c0},
                         {2, c0},
                         {3, c0},
                         {4, c0},
                         {5, c0},
                         {6, c0},
                         {7, c0}};
    inp.sn_m_set      = {{7, m}};
    inp.sn_v_set      = {{1, v1},
                         {2, v2},
                         {3, v3},
                         {4, v4},
                         {5, v5},
                         {6, v6},
                         {7, v7}};
    inp.sn_beta_set   = {{1, beta0},
                         {2, beta0},
                         {3, beta0},
                         {4, beta0},
                         {5, beta0},
                         {6, beta0},
                         {7, beta0}};
    inp.vn_fix_set    = {{0,  vol00},
                         {1,  vol01},
                         {2,  vol02},
                         {3,  vol03},
                         {4,  vol04},
                         {5,  vol05},
                         {6,  vol06},
                         {7,  vol07},
                         {8,  vol08},
                         {9,  vol09},
                         {10, vol10},
                         {11, vol11},
                         {12, vol12},
                         {13, vol13}};
    inp.vn_eta_set    = {{0,  eta0},
                         {1,  eta1},
                         {2,  eta0},
                         {3,  eta0},
                         {4,  eta0},
                         {5,  eta0},
                         {6,  eta0},
                         {7,  eta0},
                         {8,  eta0},
                         {9,  eta0},
                         {10, eta0},
                         {11, eta0},
                         {12, eta0},
                         {13, eta0}};
    VTSim calc(inp);
    calc.calc();
}

void calc_06(void){
    InputData inp;

    inp.length = 5;

    vector<double> h(inp.length,    0.0);
    vector<double> pre_tmp(inp.length, 26.0),
                   t_ex(inp.length, 15.2);
    vector<double> v1(inp.length, 339.57 / 3600),
                   v2(inp.length,  44.07 / 3600),
                   v3(inp.length,  53.46 / 3600),
                   v4(inp.length, 172.26 / 3600),
                   v5(inp.length, 127.64 / 3600),
                   v6(inp.length,  42.90 / 3600),
                   v7(inp.length, 781.90 / 3600);
    vector<double> cdtc01(inp.length, 58.82),
                   cdtc02(inp.length,  9.69),
                   cdtc03(inp.length, 22.69),
                   cdtc04(inp.length, 27.23),
                   cdtc05(inp.length, 10.30),
                   cdtc06(inp.length, 19.15);
    vector<double> h_sr(inp.length, 0.0);
    vector<double> ms(inp.length, 2.699);
    vector<int>    ac_mode(inp.length, AC_COOLING);

    inp.sts    = {SOLVE_LU, STEP_P, VENT_ERR, STEP_T, THRM_ERR, CONV_ERR, SOR_RATIO, SOR_ERR};
    inp.nodes  = {{SN_NONE, SN_NONE, SN_CALC},          //0  建築ダクト
                  {SN_NONE, SN_NONE, SN_CALC},          //1  LDK
                  {SN_NONE, SN_NONE, SN_CALC},          //2  洗面所
                  {SN_NONE, SN_NONE, SN_CALC},          //3  1階ホール
                  {SN_NONE, SN_NONE, SN_CALC},          //4  洋室大
                  {SN_NONE, SN_NONE, SN_CALC},          //5  洋室小
                  {SN_NONE, SN_NONE, SN_FIX},           //6  2Fホール
                  {SN_NONE, SN_NONE, SN_FIX},           //7  外部
                  {SN_NONE, SN_NONE, SN_NONE}};         //8　日射    
    inp.v_nets = {{0, 1, VN_FIX, h, h},         //0
                  {0, 2, VN_FIX, h, h},         //1
                  {0, 3, VN_FIX, h, h},         //2
                  {0, 4, VN_FIX, h, h},         //3
                  {0, 5, VN_FIX, h, h},         //4
                  {0, 6, VN_FIX, h, h},         //5
                  {1, 6, VN_FIX, h, h},         //6
                  {2, 6, VN_FIX, h, h},         //7
                  {3, 6, VN_FIX, h, h},         //8
                  {4, 6, VN_FIX, h, h},         //9
                  {5, 6, VN_FIX, h, h},         //10
                  {6, 0, VN_AIRCON, h, h}};        //11    
    inp.t_nets = {{1, 7, TN_SIMPLE},            //0
                  {2, 7, TN_SIMPLE},            //1
                  {3, 7, TN_SIMPLE},            //2
                  {4, 7, TN_SIMPLE},            //3
                  {5, 7, TN_SIMPLE},            //4
                  {6, 7, TN_SIMPLE},            //5
                  {1, 8, TN_SOLAR},             //6
                  {6, 0, TN_AIRCON}};           //7

    inp.sn_T_set    = {{6, pre_tmp},
                       {7, t_ex}};
    inp.sn_h_sr_set = {{8, h_sr}};

    inp.vn_fix_set  = {{0,  v1},
                       {1,  v2},
                       {2,  v3},
                       {3,  v4},
                       {4,  v5},
                       {5,  v6},
                       {6,  v1},
                       {7,  v2},
                       {8,  v3},
                       {9,  v4},
                       {10, v5},
                       {11, v7}};
    inp.tn_simple_set = {{0,  cdtc01}, 
                         {1,  cdtc02}, 
                         {2,  cdtc03}, 
                         {3,  cdtc04}, 
                         {4,  cdtc05}, 
                         {5,  cdtc06}};
    inp.tn_solar_set  = {{6, ms}};
    inp.tn_aircon_set = {{7, ac_mode, pre_tmp}};

    VTSim calc(inp);
    calc.calc();

    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, 
          vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> res; 
    vector<vector<double>>   r_p, r_c, r_t, r_qv, r_qt1, r_qt2;
        
    res = calc.result();
    r_t = get<2>(res);

    cout << endl;
    for(int i = 0; i < inp.length; i++){
        for(int j = 0; j < r_t.size(); j ++)    cout << r_t[j][i] << ",";
        cout << endl;
    }

}

int main(void){
    cout << endl << "calc 1" << endl;
    calc_01();
    cout << endl << "calc 2" << endl;
    calc_02();
    cout << endl << "calc 3" << endl;
    calc_03();
    cout << endl << "calc 4" << endl;
    calc_04();
    cout << endl << "calc 5" << endl;
    calc_05();
    cout << endl << "calc 6" << endl;
    calc_06();
    return 0;
}