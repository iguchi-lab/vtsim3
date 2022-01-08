###############################################################################
# import
###############################################################################

import numpy as np
import pandas as pd
import time

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

import vtsimc as vt

###############################################################################
# define const
###############################################################################
OPT_DF:    int  = 0              #DataFrameを出力
OPT_CSV:   int  = 1              #上記に加え、csvファイルを出力
OPT_GRAPH: int  = 2              #上記に加えグラフを描画

SN_NONE:    int = vt.SN_NONE
SN_CALC:    int = vt.SN_CALC
SN_FIX:     int = vt.SN_FIX
SN_DLY:     int = vt.SN_DLY
    
VN_SIMPLE:  int = vt.VN_SIMPLE
VN_GAP:     int = vt.VN_GAP
VN_FIX:     int = vt.VN_FIX
VN_AIRCON:  int = vt.VN_AIRCON
VN_FAN:     int = vt.VN_FAN
    
TN_SIMPLE:  int = vt.TN_SIMPLE
TN_AIRCON:  int = vt.TN_AIRCON
TN_SOLAR:   int = vt.TN_SOLAR
TN_GROUND:  int = vt.TN_GROUND
TN_HEATER:  int = vt.TN_HEATER
    
AC_AUTO:    int = vt.AC_AUTO
AC_HEATING: int = vt.AC_HEATING
AC_COOLING: int = vt.AC_COOLING
AC_STOP:    int = vt.AC_STOP

###############################################################################
# define lambda
###############################################################################

node = lambda name, v_flag, c_flag, t_flag: {'name':   name, 
                                             'v_flag': v_flag, 'c_flag': c_flag, 't_flag': t_flag}      #ノードの設定
net  = lambda name1, name2, tp:             {'name1': name1, 'name2': name2, 'type': tp}                #ネットワークの設定
r_df = lambda fn:                           pd.read_csv(fn, index_col = 0, 
                                                        parse_dates = True).fillna(method = 'bfill')\
                                                                           .fillna(method = 'ffill')     #csvファイルの読み込み
nc   = lambda id, v:                        np.array([v] * len(id))                                     #idの長さ分の値value
nd   = lambda df, cl:                       np.array(df[cl])                                            #dfの列clを設定
ix   = lambda length:                       pd.date_range(datetime(2021, 1, 1, 0, 0, 0), 
                                                          datetime(2021, 1, 1, 0, 0, 0) + timedelta(seconds = length), 
                                                          freq='1s')                                    #長さlength、1s毎の時刻
d_node  = lambda name:                      name + '_c'                                                 #遅延ノードの名前作成



calc = vt.VTSim()
sts  = vt.CalcStatus()
inp  = vt.InputData()

node = {}

###############################################################################
# define function
###############################################################################

def run_calc(ix, sn, **kwargs):                                                                            #はじめに呼び出される関数
    
    set_calc_status(ix, **kwargs)
    set_node_net(sn, **kwargs)
    
    print('Start vtsim calc.')
    s_time = time.time()
    calc.calc()
    e_time = time.time() - s_time    
    print('Finish vtsim calc.')
    print("calc time = {0}".format(e_time * 1000) + "[ms]")
    
    opt = kwargs['output'] if 'output' in kwargs else OPT_GRAPH                                     #出力フラグ
    return output_calc(calc.result(), ix, opt)

def set_calc_status(ix, **kwargs):
    sts.length = len(ix)
    sts.t_step = (ix[1] - ix[0]).seconds + (ix[1] - ix[0]).microseconds / 1000000                          #t_stepの読み込み    

    if 'solve'     in kwargs:   sts.solve     = kwargs['solve']
    if 'step_p'    in kwargs:   sts.step_pe   = kwargs['step_p']
    if 'vent_err'  in kwargs:   sts.vent_err  = kwargs['vent_err']
    if 'step_t'    in kwargs:   sts.step_t    = kwargs['step_t']
    if 'thrm_err'  in kwargs:   sts.thrm_err  = kwargs['thrm_err']
    if 'conv_err'  in kwargs:   sts.conv_err  = kwargs['conv_err']
    if 'sor_ratio' in kwargs:   sts.sor_ratio = kwargs['sor_ratio']
    if 'sor_err'   in kwargs:   sts.sor_err   = kwargs['sor_err']

    calc.set_calc_status(sts)

def set_node_net(sn, **kwargs):
    vn         = kwargs['vn']     if 'vn'  in kwargs else []                                               #vnの読み込み
    tn         = kwargs['tn']     if 'tn'  in kwargs else []                                               #tnの読み込み

    nodes, v_nets, t_nets                                         = [], [], []
    sn_P_set, sn_C_set, sn_T_set, sn_h_sr_set, sn_h_inp_set       = [], [], [], [], []
    sn_v_set, sn_capa_set, sn_m_set, sn_beta_set                  = [], [], [], []
    vn_simple_set, vn_gap_set, vn_fix_set, vn_fan_set, vn_eta_set = [], [], [], [], []
    tn_simple_set, tn_aircon_set, tn_solar_set, tn_ground_set     = [], [], [], []

    #リストかnp.ndarrayでなければlength分の長さのリストにする
    to_list_f = lambda v:   [float(v)] * sts.length if type(v) != list and type(v) != np.ndarray else v  
    to_list_i = lambda v:   [int(v)]   * sts.length if type(v) != list and type(v) != np.ndarray else v  


    for i, n in enumerate(sn):                                                                              #sn
        node[n['name']] = i                                                                                 #ノード番号
        
        v_flag = n['v_flag'] if 'v_flag' in n else SN_NONE
        c_flag = n['c_flag'] if 'c_flag' in n else SN_NONE
        t_flag = n['t_flag'] if 't_flag' in n else SN_NONE

        calc.sn[i] = vt.Node(sts.length, i, [v_flag, c_flag, t_flag])
        if 'p'     in n:    sn[i].p     = to_list_f(n['p'])
        if 'c'     in n:    sn[i].c     = to_list_f(n['c'])
        if 't'     in n:    sn[i].t     = to_list_f(n['t'])
        if 'h_sr'  in n:    sn[i].h_sr  = to_list_f(n['h_sr'])
        if 'h_inp' in n:    sn[i].h_inp = to_list_f(n['h_inp'])
        if 'v'     in n:    sn[i].v     = to_list_f(n['v'])
        if 'm'     in n:    sn[i].m     = to_list_f(n['m'])
        if 'beta'  in n:    sn[i].beta  = to_list_f(n['beta'])

        """"
        nodes.append([v_flag, c_flag, t_flag])                                                              #計算フラグ

        if 'p' in n:            sn_P_set.append([i, to_list_f(n['p'],     sts.length)])                     #圧力、行列で設定可能                                                 
        if 'c' in n:            sn_C_set.append([i, to_list_f(n['c'],     sts.length)])                     #濃度、行列で設定可能
        if 't' in n:            sn_T_set.append([i, to_list_f(n['t'],     sts.length)])                     #温度、行列で設定可能
        if 'h_sr' in n:      sn_h_sr_set.append([i, to_list_f(n['h_sr'],  sts.length)])                     #日射量、行列で設定可能
        if 'h_inp' in n:    sn_h_inp_set.append([i, to_list_f(n['h_inp'], sts.length)])                     #発熱、行列で設定可能
        if 'v' in n:            sn_v_set.append([i, to_list_f(n['v'],     sts.length)])                     #気積、行列で設定可能
        if 'm' in n:            sn_m_set.append([i, to_list_f(n['m'],     sts.length)])                     #発生量、行列で設定可能
        if 'beta' in n:      sn_beta_set.append([i, to_list_f(n['beta'],  sts.length)])                     #濃度減少率、行列で設定可能
        """

    for i, nt in enumerate(vn):                                                                             #vn
        h1 = to_list_f(nt['h1'], sts.length) if 'h1' in nt else to_list_f(0.0, sts.length)                  #高さ1、行列設定不可
        h2 = to_list_f(nt['h2'], sts.length) if 'h2' in nt else to_list_f(0.0, sts.length)                  #高さ2、行列設定不可
        
        vn_type = nt['type'] if 'type' in nt else VN_FIX
        v_nets.append([node[nt['name1']], node[nt['name2']], vn_type, h1, h2])                           #ネットワークタイプ＆高さ
        
        if vn_type == VN_FIX:           vn_fix_set.append([i, to_list_f(nt['vol'],   sts.length)])       #風量固定値、行列で設定可能
        if vn_type == VN_AIRCON:        vn_fix_set.append([i, to_list_f(nt['vol'],   sts.length)])       #風量固定値、行列で設定可能
        if vn_type == VN_SIMPLE:     vn_simple_set.append([i, to_list_f(nt['alpha'], sts.length), 
                                                              to_list_f(nt['area'],  sts.length)])       #単純開口、行列で設定可能
        if vn_type == VN_GAP:           vn_gap_set.append([i, to_list_f(nt['a'],     sts.length), 
                                                              to_list_f(nt['n'],     sts.length)])       #隙間、行列で設定可能
        if vn_type == VN_FAN:           vn_fan_set.append([i, to_list_f(nt['qmax'],  sts.length), 
                                                              to_list_f(nt['pmax'],  sts.length), 
                                                              to_list_f(nt['q1'],    sts.length),
                                                              to_list_f(nt['p1'],    sts.length)])       #ファン、行列で設定可能

        if 'eta' in nt:                    vn_eta_set.append([i, to_list_f(nt['eta'],   sts.length)])               
        else:                              vn_eta_set.append([i, to_list_f(0.0,         sts.length)])       #粉じん除去率、行列で設定可能

    for i, nt in enumerate(tn):                                                                             #tn
        tn_type = nt['type'] if 'type' in nt else TN_SIMPLE
        t_nets.append([node[nt['name1']], node[nt['name2']], tn_type])                                   #ネットワークタイプ
        
        if tn_type == TN_SIMPLE:     tn_simple_set.append([i, to_list_f(nt['cdtc'],    sts.length)])     #コンダクタンス、行列設定可能
        if tn_type == TN_AIRCON:     tn_aircon_set.append([i, to_list_i(nt['ac_mode'], sts.length), 
                                                              to_list_f(nt['pre_tmp'], sts.length)])     #エアコン運転モード
        if tn_type == TN_SOLAR:       tn_solar_set.append([i, to_list_f(nt['ms'],      sts.length)])     #日射熱取得率、行列設定可能
        if tn_type == TN_GROUND:     tn_ground_set.append([i, to_list_f(nt['area'],    sts.length),           
                                                              to_list_f(nt['rg'],      sts.length), 
                                                              nt['phi_0'], nt['cof_r'], nt['cof_phi']])  #地盤熱応答、行列設定不可（面積と断熱性能はOK）
        
    for i, n in enumerate([n for n in sn if 'capa' in n]):                                                  #熱容量の設定のあるノード
        node[d_node(n['name'])] = len(sn) + i                                                               #時間遅れノードのノード番号
        nodes.append([SN_NONE, SN_NONE, SN_DLY])                                                            #計算フラグ
        sn_capa_set.append([node[d_node(n['name'])], node[n['name']]])                                      #熱容量の設定
        if 't' in n:    sn_T_set.append([len(sn) + i, to_list_f(n['t'], sts.length)])
        t_nets.append([node[n['name']], node[d_node(n['name'])], TN_SIMPLE])                                #ネットワークの設定
        tn_simple_set.append([len(tn) + i, to_list_f(n['capa'] / sts.t_step, sts.length)])                  #コンダクタンス（熱容量）

    inp.nodes, inp.v_nets, inp.t_nets                                                 = nodes, v_nets, t_nets
    inp.sn_P_set, inp.sn_C_set, inp.sn_T_set, inp.sn_h_sr_set, inp.sn_h_inp_set       = sn_P_set, sn_C_set, sn_T_set, sn_h_sr_set, sn_h_inp_set  
    inp.sn_v_set, inp.sn_capa_set, inp.sn_m_set, inp.sn_beta_set                      = sn_v_set, sn_capa_set, sn_m_set, sn_beta_set
    inp.vn_simple_set, inp.vn_gap_set, inp.vn_fix_set, inp.vn_fan_set, inp.vn_eta_set = vn_simple_set, vn_gap_set, vn_fix_set, vn_fan_set, vn_eta_set
    inp.tn_simple_set, inp.tn_aircon_set, inp.tn_solar_set, inp.tn_ground_set         = tn_simple_set, tn_aircon_set, tn_solar_set, tn_ground_set      

    calc.set_inp(inp)

def output_calc(res, ix, opt):
    print('Create pd.DataFrames')

    node_swap = {v: k for k, v in node.items()}
    n_columns = [node_swap[i] for i in range(len(inp.nodes))]                                                                           #出力用カラムの作成（ノード）
    v_columns = [str(i) + " " + node_swap[inp.v_nets[i][0]] + "->" + node_swap[inp.v_nets[i][1]] for i in range(len(inp.v_nets))]       #出力用カラムの作成（換気回路網）
    t_columns = [str(i) + " " + node_swap[inp.t_nets[i][0]] + "->" + node_swap[inp.t_nets[i][1]] for i in range(len(inp.t_nets))]       #出力用カラムの作成（熱回路網）
    
    dat_list  = [{'df': pd.DataFrame(), 'columns': n_columns, 'fn': 'vent_p.csv',   'title': '圧力',  'unit': '[Pa]'},
                 {'df': pd.DataFrame(), 'columns': n_columns, 'fn': 'vent_c.csv',   'title': '濃度',  'unit': '[個/L]'},
                 {'df': pd.DataFrame(), 'columns': n_columns, 'fn': 'them_t.csv',   'title': '温度',  'unit': '[℃]'},
                 {'df': pd.DataFrame(), 'columns': v_columns, 'fn': 'vent_qv.csv',  'title': '風量',  'unit': '[m3/s]'},
                 {'df': pd.DataFrame(), 'columns': v_columns, 'fn': 'thrm_qt1.csv', 'title': '熱量1', 'unit': '[W]'},
                 {'df': pd.DataFrame(), 'columns': t_columns, 'fn': 'thrm_qt2.csv', 'title': '熱量2', 'unit': '[W]'}]
    
    for i, d in enumerate(dat_list):
        if len(d) != 0: d['df'] = pd.DataFrame(np.array(res[i]).T,  index = ix, columns = d['columns'])

    if opt > 0:
        print('Output csv files.')
        for d in dat_list:      d['df'].to_csv(d['fn'], encoding = 'utf_8_sig')

    if opt > 1:
        print('Draw Graphs.')
        fig = plt.figure(facecolor = 'w', figsize = (18, len(dat_list) * 4))
        fig.subplots_adjust(wspace = -0.1, hspace=0.9)

        for i, d in enumerate(dat_list):
            a = fig.add_subplot(len(dat_list), 1, i + 1)
            for cl in d['df'].columns:
                a.plot(d['df'][cl], linewidth = 1.0, label = cl)
            a.legend(ncol = 5, bbox_to_anchor = (0, 1.05, 1, 0), 
                     loc = 'lower right', borderaxespad = 0, facecolor = 'w', edgecolor = 'k')
            a.set_title(d['title'], loc='left')
            a.set_ylabel(d['unit'])

    return dat_list[0]['df'], dat_list[1]['df'], dat_list[2]['df'], dat_list[3]['df'], dat_list[4]['df'], dat_list[5]['df']