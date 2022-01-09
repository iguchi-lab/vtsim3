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

#リストかnp.ndarrayでなければlength分の長さのリストにする
to_list_f = lambda v:   [float(v)] * sts.length if type(v) != list and type(v) != np.ndarray else v  
to_list_i = lambda v:   [int(v)]   * sts.length if type(v) != list and type(v) != np.ndarray else v 

calc = vt.VTSim()
sts  = vt.CalcStatus()

node = {}

###############################################################################
# define function
###############################################################################

def run_calc(ix, sn, **kwargs):                                                     #はじめに呼び出される関数
    
    set_calc_status(ix, **kwargs)
    
    set_sim_node(sn, **kwargs)
    set_vent_net(**kwargs)
    set_thrm_net(**kwargs)

    add_capa(sn)

    print('Start vtsim calc.')
    s_time = time.time()
    calc.calc()
    e_time = time.time() - s_time    
    print('Finish vtsim calc.')
    print("calc time = {0}".format(e_time * 1000) + "[ms]")
    
    opt = kwargs['output'] if 'output' in kwargs else OPT_GRAPH                     #出力フラグ
    return output_calc(calc.result(), ix, opt)

def set_calc_status(ix, **kwargs):
    sts.length = len(ix)
    sts.t_step = (ix[1] - ix[0]).seconds + (ix[1] - ix[0]).microseconds / 1000000   #t_stepの読み込み    

    if 'solve'     in kwargs:   sts.solve     = kwargs['solve']
    if 'step_p'    in kwargs:   sts.step_pe   = kwargs['step_p']
    if 'vent_err'  in kwargs:   sts.vent_err  = kwargs['vent_err']
    if 'step_t'    in kwargs:   sts.step_t    = kwargs['step_t']
    if 'thrm_err'  in kwargs:   sts.thrm_err  = kwargs['thrm_err']
    if 'conv_err'  in kwargs:   sts.conv_err  = kwargs['conv_err']
    if 'sor_ratio' in kwargs:   sts.sor_ratio = kwargs['sor_ratio']
    if 'sor_err'   in kwargs:   sts.sor_err   = kwargs['sor_err']

    calc.sts = sts

def set_sim_node(sn, **kwargs):
    v_idc, c_idc, t_idc = [], [], []

    for i, n in enumerate(sn):                                                      #sn
        node[n['name']] = i                                                         #ノード番号
        
        v_flag = n['v_flag'] if 'v_flag' in n else SN_NONE
        c_flag = n['c_flag'] if 'c_flag' in n else SN_NONE
        t_flag = n['t_flag'] if 't_flag' in n else SN_NONE                          #計算フラグ

        calc.sn_add(i, [v_flag, c_flag, t_flag])

        if v_flag == SN_CALC:   v_idc.append(i)
        if c_flag == SN_CALC:   c_idc.append(i)
        if t_flag == SN_CALC:   t_idc.append(i)

        if 'p'     in n:    calc.sn[i].p     = to_list_f(n['p'])                    #圧力、行列で設定可能
        if 'c'     in n:    calc.sn[i].c     = to_list_f(n['c'])                    #濃度、行列で設定可能
        if 't'     in n:    calc.sn[i].t     = to_list_f(n['t'])                    #温度、行列で設定可能
        if 'h_sr'  in n:    calc.sn[i].h_sr  = to_list_f(n['h_sr'])                 #日射量、行列で設定可能
        if 'h_inp' in n:    calc.sn[i].h_inp = to_list_f(n['h_inp'])                #発熱、行列で設定可能
        if 'v'     in n:    calc.sn[i].v     = to_list_f(n['v'])                    #気積、行列で設定可能
        if 'm'     in n:    calc.sn[i].m     = to_list_f(n['m'])                    #発生量、行列で設定可能
        if 'beta'  in n:    calc.sn[i].beta  = to_list_f(n['beta'])                 #濃度減少率、行列で設定可能

    calc.v_idc = v_idc
    calc.c_idc = c_idc
    calc.t_idc = t_idc  

def set_vent_net(**kwargs):
    vn         = kwargs['vn']     if 'vn'  in kwargs else []                        #vnの読み込み
    for i, nt in enumerate(vn):
        h1 = to_list_f(nt['h1']) if 'h1' in nt else to_list_f(0.0)                  #高さ1、行列設定不可
        h2 = to_list_f(nt['h2']) if 'h2' in nt else to_list_f(0.0)                  #高さ2、行列設定不可
        vn_type = nt['type'] if 'type' in nt else VN_FIX

        calc.vn_add(i, node[nt['name1']], node[nt['name2']], vn_type, h1, h2)
        
        if (vn_type == VN_FIX) or (vn_type == VN_AIRCON):       
            calc.vn[i].qv = to_list_f(nt['vol'])                                       #風量固定値、行列で設定可能
        if vn_type == VN_SIMPLE:    
            calc.vn[i].alpha = to_list_f(nt['alpha'])
            calc.vn[i].area  = to_list_f(nt['area'])                                     #単純開口、行列で設定可能
        if vn_type == VN_GAP:           
            calc.vn[i].a     = to_list_f(nt['a'])
            calc.vn[i].n     = to_list_f(nt['n'])                                        #隙間、行列で設定可能
        if vn_type == VN_FAN:           
            calc.vn[i].qmax = to_list_f(nt['qmax']) 
            calc.vn[i].pmax = to_list_f(nt['pmax']) 
            calc.vn[i].q1   = to_list_f(nt['q1'])
            calc.vn[i].p1   = to_list_f(nt['p1'])                                       #ファン、行列で設定可能
        if vn_type == VN_AIRCON:
            calc.i_vn_ac = i
        calc.vn[i].eta = to_list_f(nt['eta']) if 'eta' in nt else to_list_f(0.0)        #粉じん除去率、行列で設定可能

def set_thrm_net(**kwargs):
    tn         = kwargs['tn']     if 'tn'  in kwargs else []                            #tnの読み込み
    for i, nt in enumerate(tn):                                                         #tn
        tn_type = nt['type'] if 'type' in nt else TN_SIMPLE

        calc.tn_add(i, node[nt['name1']], node[nt['name2']], tn_type)
        
        if tn_type == TN_SIMPLE:     
            calc.tn[i].cdtc = to_list_f(nt['cdtc'])                                     #コンダクタンス、行列設定可能
        if tn_type == TN_AIRCON:     
            calc.tn[i].ac_mode = to_list_i(nt['ac_mode']) 
            calc.tn[i].pre_tmp = to_list_f(nt['pre_tmp'])                               #エアコン運転モード
            calc.i_tn_ac = i
        if tn_type == TN_SOLAR:       
            calc.tn[i].ms      = to_list_f(nt['ms'])                                    #日射熱取得率、行列設定可能
        if tn_type == TN_GROUND:     
            calc.tn[i].area    = to_list_f(nt['area'])           
            calc.tn[i].rg      = to_list_f(nt['rg']) 
            calc.tn[i].phi_0   = nt['phi_0']
            calc.tn[i].cof_r   = nt['cof_r']
            calc.tn[i].cof_phi = nt['cof_phi']                                          #地盤熱応答、行列設定不可（面積と断熱性能はOK）

def add_capa(sn):  
    for i, n in enumerate([n for n in calc.sn if 'capa' in n]):                         #熱容量の設定のあるノード
        node[d_node(n['name'])] = len(sn) + i                                           #時間遅れノードのノード番号
        
        calc.sn_add(len(calc.sn) + i, [SN_NONE, SN_NONE, SN_DLY])                       #計算フラグ
        if 't' in n:    sn[len(calc.sn) + i] = to_list_f(n['t'])

        calc.tn_add(len(calc.tn) + i, d_node(n['name']), node[n['name']], TN_SIMPLE)    #熱容量の設定
        calc.tn[i].cdtc = to_list_f(n['capa'] / sts.t_step)                             #コンダクタンス（熱容量）            

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