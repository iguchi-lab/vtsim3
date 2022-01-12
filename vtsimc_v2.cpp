#include <stdio.h>
#include <math.h>                   //数学　ヘッダーファイルの読み込み
#include <iostream>                 //出力　ヘッダーファイルの読み込み
#include <iomanip>
#include <vector>
#include <tuple>
#include <map>

#include "vtsimc.h"

using namespace std;

void calc_01(void){
    CalcStatus sts;
    VTSim calc;

    sts.length = 2;
    calc.sts = sts;

    vector<double> vol(sts.length, 100.0 / 3600.0),
                   alpha(sts.length, 0.6),
                   area(sts.length, 1.0),
                   h(sts.length, 0.0);

    calc.sn_add(0, {SN_FIX, SN_NONE, SN_NONE});
    calc.sn_add(1, {SN_CALC, SN_NONE, SN_NONE});
    calc.sn_add(2, {SN_CALC, SN_NONE, SN_NONE});

    calc.vn_add(0, 0, 1, VN_FIX, h, h);
    calc.vn_add(1, 1, 2, VN_SIMPLE, h, h);
    calc.vn_add(2, 2, 0, VN_SIMPLE, h, h);

    calc.vn[1].alpha = alpha;
    calc.vn[1].area  = area;

    calc.vn[2].alpha = alpha;
    calc.vn[2].area  = area;

    calc.vn[0].qv = vol;

    calc.calc();
}

void calc_02(void){
    CalcStatus sts;
    VTSim calc;

    sts.length = 2;
    calc.sts = sts;

    vector<double> t1(sts.length, 20.0),              
                   t2(sts.length,  0.0),
                   vol1(sts.length, 1000.0 / 3600.0), 
                   vol2(sts.length,  500.0 / 3600.0),
                   h_sr(sts.length, 4000.0),             
                   h_inp(sts.length, 2000.0),
                   cdtc1(sts.length, 100.0 * 0.87),
                   cdtc2(sts.length, 200.0 * 0.87),
                   ms(sts.length, 1.0),
                   h(sts.length, 0.0);
    vector<int>    ac_mode(sts.length, AC_COOLING);

    calc.sn_add(0, {SN_NONE, SN_NONE, SN_CALC});                                                     //0: AC-out
    calc.sn_add(1, {SN_NONE, SN_NONE, SN_CALC});                                                     //1: Room1
    calc.sn_add(2, {SN_NONE, SN_NONE, SN_CALC});                                                     //2: Room2
    calc.sn_add(3, {SN_NONE, SN_NONE, SN_FIX});                                                      //3: AC-in
    calc.sn_add(4, {SN_NONE, SN_NONE, SN_FIX});                                                      //4: Outside
    calc.sn_add(5, {SN_NONE, SN_NONE, SN_NONE});                                                     //5: solar
    calc.sn_add(6, {SN_NONE, SN_NONE, SN_NONE});                                                     //6: heater  

    calc.vn_add(0, 0, 1, VN_FIX,    h, h); 
    calc.vn_add(1, 0, 2, VN_FIX,    h, h);
    calc.vn_add(2, 1, 3, VN_FIX,    h, h);
    calc.vn_add(3, 2, 3, VN_FIX,    h, h);
    calc.vn_add(4, 3, 0, VN_FIX,    h, h);
    
    calc.tn_add(0, 1, 4, TN_SIMPLE); 
    calc.tn_add(1, 2, 4, TN_SIMPLE); 
    calc.tn_add(2, 3, 0, TN_AIRCON);
    calc.tn_add(3, 2, 5, TN_SOLAR);
    calc.tn_add(4, 2, 6, TN_HEATER);

    calc.sn[3].t = t1;
    calc.sn[4].t = t2;

    calc.sn[5].h_sr  = h_sr;
    calc.sn[6].h_inp = h_inp;
    
    calc.vn[0].qv = vol2;
    calc.vn[1].qv = vol2;
    calc.vn[2].qv = vol2;
    calc.vn[3].qv = vol2;
    calc.vn[4].qv = vol2;
    
    calc.tn[0].cdtc = cdtc1;
    calc.tn[1].cdtc = cdtc2;

    calc.tn[2].ac_mode = ac_mode;
    calc.tn[2].pre_tmp = t1;

    calc.tn[3].ms = ms;

    calc.calc();
}
/*
void calc_03(void){
    CalcStatus sts;
    InputData inp;
    VTSim calc;

    sts.length  = 3;
    sts.sor_err = 0.01;
    calc.set_calc_status(sts);

    vector<double> alpha(sts.length, 0.6),
                   area(sts.length, 0.21),
                   a1(sts.length, 0.005516666666666667),
                   a2(sts.length, 0.018633333333333335),
                   a3(sts.length, 0.005516666666666667),
                   a4(sts.length, 0.006208333333333334),
                   a5(sts.length, 0.004141666666666667),
                   n(sts.length, 1.5),
                   cdtc01(sts.length, 15.292200000000001),
                   cdtc02(sts.length, 51.6516),
                   cdtc03(sts.length, 15.292200000000001),
                   cdtc04(sts.length, 17.209500000000002),
                   cdtc05(sts.length, 11.4807),
                   cdtc06(sts.length, 1.7745),
                   cdtc07(sts.length, 3.549),
                   cdtc08(sts.length, 2.4843),
                   cdtc09(sts.length, 2.4843),
                   cdtc10(sts.length, 1.2441),
                   h(sts.length, 0.0),
                   g_area(sts.length, 6.62),
                   rg(sts.length, 1.85);


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

    calc.set_inp(inp);
    calc.calc();
}
*/
int get_xyz(int x, int y, int z){
    return x * 15 + y * 5 + z;
}

void calc_04(void){
    CalcStatus sts;
    VTSim calc;

    int X = 3, Y = 3, Z = 5;
    double cx[3]   = {1.52, 1.52, 1.52}, cy[3] = {1.02, 1.60, 1.02}, cz[5] = {0.48, 0.48, 0.48, 0.48, 0.48};
    double hcz[5]  = {0.24, 0.72, 1.20, 1.68, 2.16};
    double htz[5]  = {0.48, 0.96, 1.44, 1.92};
    double cdtc[6] = {0.4, 0.5, 0.94, 1.88, 3.49, 5.59};

    sts.length    = 2;
    sts.t_step    = 0.1;

    sts.step_p    = 1.0e-6;
    sts.vent_err  = 1.0e-1;
    
    sts.step_t    = 1.0e-6;
    sts.thrm_err  = 1.0e-1;

    sts.conv_err  = 1.0e-3;
    sts.sor_ratio = 0.8;
    sts.sor_err   = 1.0e-6;
    calc.sts = sts;

    vector<double> Ti(sts.length, 20.0), 
                   To(sts.length, 6.0), 
                   H_inp(sts.length, 462.3),
                   alpha(sts.length, 1.0);

    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int z = 0; z < Z; z++){
                calc.sn_add(get_xyz(x, y, z), {SN_CALC, SN_NONE, SN_CALC});
                calc.sn[get_xyz(x, y, z)].t = Ti;
            }
        }
    }

    calc.sn_add(45, {SN_NONE, SN_NONE, SN_FIX});
    calc.sn[45].t = To;
    calc.sn_add(46, {SN_NONE, SN_NONE, SN_FIX});
    calc.sn[46].t = Ti;
    calc.sn_add(47, {SN_NONE, SN_NONE, SN_NONE});
    calc.sn[47].h_inp = H_inp;

    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int z = 0; z < Z; z++){
                calc.sn_add(get_xyz(x, y, z) + 48, {SN_NONE, SN_NONE, SN_DLY});
                calc.sn[get_xyz(x, y, z) + 48].t = Ti;
                calc.sn[get_xyz(x, y, z) + 48].s_i = get_xyz(x, y, z);
            }
        }
    }

    int v_n = 0, t_n = 0;
    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int z = 0; z < Z; z++){
                if(x != X - 1){
                    vector<double> h(sts.length, hcz[z]);
                    calc.vn_add(v_n, get_xyz(x, y, z), get_xyz(x + 1, y, z), VN_SIMPLE, h, h);
                    double area = cy[y] * cz[z];
                    calc.vn[v_n].alpha = alpha;
                    calc.vn[v_n].area = {area, area, area, area};
                    v_n++;
                }
                if(y != Y - 1){
                    vector<double> h(sts.length, hcz[z]);   
                    calc.vn_add(v_n, get_xyz(x, y, z), get_xyz(x, y + 1, z), VN_SIMPLE, h, h);
                    double area = cx[x] * cz[z];
                    calc.vn[v_n].alpha = alpha;
                    calc.vn[v_n].area = {area, area, area, area};
                    v_n++;
                }
                if(z != Z - 1){
                    vector<double> h(sts.length, htz[z]);   
                    calc.vn_add(v_n, get_xyz(x, y, z), get_xyz(x, y, z + 1), VN_SIMPLE, h, h);
                    double area = cx[x] * cy[y];
                    calc.vn[v_n].alpha = alpha;
                    calc.vn[v_n].area = {area, area, area, area};
                    v_n++;
                }
                if(z == 0 || z == Z - 1){  
                    calc.tn_add(t_n, get_xyz(x, y, z), 45, TN_SIMPLE);
                    double c = cdtc[0];
                    calc.tn[t_n].cdtc = {c, c, c, c};
                    t_n++;
                }
                if(y == 0 || y == Y - 1){   
                    calc.tn_add(t_n, get_xyz(x, y, z), 45, TN_SIMPLE);
                    double c = cdtc[0];
                    calc.tn[t_n].cdtc = {c, c, c, c};
                    t_n++;
                }
                if(x == 0){
                    if(y == 1 && z < 4){
                        calc.tn_add(t_n, get_xyz(x, y, z), 54, TN_SIMPLE);
                        double c = cdtc[5];
                        calc.tn[t_n].cdtc = {c, c, c, c};
                        t_n++;
                    }
                    else{
                        calc.tn_add(t_n, get_xyz(x, y, z), 45, TN_SIMPLE);
                        double c = cdtc[0];
                        calc.tn[t_n].cdtc = {c, c, c, c};
                        t_n++;
                    }
                }
                if(x == X - 1){   
                    calc.tn_add(t_n, get_xyz(x, y, z), 46, TN_SIMPLE);
                    double c = cdtc[0];
                    calc.tn[t_n].cdtc = {c, c, c, c};
                    t_n++;
                }
                calc.tn_add(t_n, get_xyz(x, y, z), get_xyz(x, y, z) + 48, TN_SIMPLE);
                double c = cx[x] * cy[y] * cz[z] * 1.205 * 1006 / sts.t_step;
                calc.tn[t_n].cdtc = {c, c, c, c};
                t_n++;                
            }
        }
    }
    calc.tn_add(t_n, get_xyz(2, 1, 0), 47, TN_HEATER);

    calc.calc();
}

/*
void calc_05(void){
    CalcStatus sts;
    InputData inp;
    VTSim calc;

    sts.length        = 3;
    sts.t_step        = 10;
    calc.set_calc_status(sts);

    vector<double>    c0(sts.length, 40000.0), 
                       m(sts.length, 28470.0),
                      v1(sts.length,  13.66),
                      v2(sts.length,   6.00),
                      v3(sts.length,  40.25),
                      v4(sts.length, 126.49),
                      v5(sts.length,  16.56),
                      v6(sts.length,   9.94),
                      v7(sts.length, 100.46),
                   beta0(sts.length, 0.00005),
                   vol00(sts.length,  320.0 / 3600),
                   vol01(sts.length, 1206.0 / 3600),
                   vol02(sts.length,  155.0 / 3600),
                   vol03(sts.length,  313.0 / 3600),
                   vol04(sts.length,  108.0 / 3600),
                   vol05(sts.length,  176.0 / 3600),
                   vol06(sts.length,  454.0 / 3600),
                   vol07(sts.length,  155.0 / 3600),
                   vol08(sts.length,  313.0 / 3600),
                   vol09(sts.length,  108.0 / 3600),
                   vol10(sts.length,  100.0 / 3600),
                   vol11(sts.length,   76.0 / 3600),
                   vol12(sts.length,  220.0 / 3600),
                   vol13(sts.length,  886.0 / 3600),              
                    eta0(sts.length,    1.0),
                    eta1(sts.length,    0.9),
                       h(sts.length,    0.0);

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

    calc.set_inp(inp);
    calc.calc();
}

void calc_06(void){
    CalcStatus sts;
    InputData inp;

    VTSim calc;
    sts.length        = 5;
    calc.set_calc_status(sts);

    vector<double> h(sts.length,    0.0);
    vector<double> pre_tmp(sts.length, 26.0),
                   t_ex(sts.length, 15.2);
    vector<double> v1(sts.length, 339.57 / 3600),
                   v2(sts.length,  44.07 / 3600),
                   v3(sts.length,  53.46 / 3600),
                   v4(sts.length, 172.26 / 3600),
                   v5(sts.length, 127.64 / 3600),
                   v6(sts.length,  42.90 / 3600),
                   v7(sts.length, 781.90 / 3600);
    vector<double> cdtc01(sts.length, 58.82),
                   cdtc02(sts.length,  9.69),
                   cdtc03(sts.length, 22.69),
                   cdtc04(sts.length, 27.23),
                   cdtc05(sts.length, 10.30),
                   cdtc06(sts.length, 19.15);
    vector<double> h_sr(sts.length, 0.0);
    vector<double> ms(sts.length, 2.699);
    vector<int>    ac_mode(sts.length, AC_COOLING);

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

    calc.set_inp(inp);
    calc.calc();

    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, 
          vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> res; 
    vector<vector<double>>   r_p, r_c, r_t, r_qv, r_qt1, r_qt2;
        
    res = calc.result();
    r_t = get<2>(res);

    cout << endl;
    for(int i = 0; i < sts.length; i++){
        for(int j = 0; j < r_t.size(); j ++)    cout << r_t[j][i] << ",";
        cout << endl;
    }

}
*/
int main(void){
    cout << endl << "calc 1" << endl;
    calc_01();
    cout << endl << "calc 2" << endl;
    calc_02();
    /*
    cout << endl << "calc 3" << endl;
    calc_03();
    */
    cout << endl << "calc 4" << endl;
    calc_04();
    /*
    cout << endl << "calc 5" << endl;
    calc_05();
    cout << endl << "calc 6" << endl;
    calc_06();
    */
    return 0;
}