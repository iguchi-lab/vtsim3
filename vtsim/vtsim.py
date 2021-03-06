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

SOLVE_LU:   int = vt.SOLVE_LU
SOLVE_SOR:  int = vt.SOLVE_SOR

OPT_DF:     int  = 0              #DataFrameを出力
OPT_CSV:    int  = 1              #上記に加え、csvファイルを出力
OPT_GRAPH:  int  = 2              #上記に加えグラフを描画

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
read_csv = lambda fn:                       pd.read_csv(fn, index_col = 0, 
                                                        parse_dates = True).fillna(method = 'bfill')\
                                                                           .fillna(method = 'ffill')     #csvファイルの読み込み

nc   = lambda id, v:                        np.array([v] * len(id))                                     #idの長さ分の値value
nd   = lambda df, cl:                       np.array(df[cl])                                            #dfの列clを設定
ix   = lambda length:                       pd.date_range(datetime(2021, 1, 1, 0, 0, 0), 
                                                          datetime(2021, 1, 1, 0, 0, 0) + timedelta(seconds = length), 
                                                          freq='1s')                                    #長さlength、1s毎の時刻
d_node  = lambda name:                      name + '_c'                                                 #遅延ノードの名前作成

calc = vt.VTSim()

df_p   = pd.DataFrame()
df_c   = pd.DataFrame()
df_t   = pd.DataFrame()
df_qv  = pd.DataFrame()
df_qt1 = pd.DataFrame()
df_qt2 = pd.DataFrame()

###############################################################################
# define function
###############################################################################

def to_list_f(v):
    if   type(v) == list:                   return(v)
    elif type(v) == np.ndarray:             return(v)
    elif type(v) == pd.core.series.Series:  return(np.array(v))
    else:                                   return[float(v)] * calc.sts.length

def to_list_i(v):
    if   type(v) == list:                   return(v)
    elif type(v) == np.ndarray:             return(v)
    elif type(v) == pd.core.series.Series:  return(np.array(v))
    else:                                   return[int(v)] * calc.sts.length

def run_calc(ix, sn, **kwargs):                                                     #はじめに呼び出される関数    
    
    sn_c, vn_c, tn_c = [], [], []
    
    print('Set calc status.')
    set_calc_status(ix, **kwargs)
    
    print('Set SimNode.')
    sn_c = set_sim_node(sn)
    
    vn = kwargs['vn']     if 'vn'  in kwargs else []
    print('Set VentNet.')
    vn_c = set_vent_net(vn)

    tn = kwargs['tn']     if 'tn'  in kwargs else []    
    print('Set ThrmNet.')
    sn_c, tn_c = set_thrm_net(sn, tn, sn_c)

    print('ready')
    print('sts     ', calc.sts)
    print('sn      ', calc.sn)
    print('node    ', calc.node)
    print('vn      ', calc.vn)
    print('tn      ', calc.tn)
    print('v_idc   ', calc.v_idc)
    print('c_idc   ', calc.c_idc)
    print('t_idc   ', calc.t_idc)
    print('i_vn_ac ', calc.i_vn_ac)
    print('i_tn_ac ', calc.i_tn_ac)

    print('Start vtsim calc.')
    s_time = time.time()
    calc.calc()
    e_time = time.time() - s_time    
    print('Finish vtsim calc.')
    print("calc time = {0}".format(e_time * 1000) + "[ms]")
    
    opt = kwargs['output'] if 'output' in kwargs else OPT_GRAPH                     #出力フラグ

    return output_calc(ix, opt, calc.result(), sn_c, vn_c, tn_c)

def set_calc_status(ix, **kwargs):
    sts  = vt.CalcStatus()

    sts.length = len(ix)
    sts.t_step = (ix[1] - ix[0]).seconds + (ix[1] - ix[0]).microseconds / 1000000   #t_stepの読み込み    

    if 'solve'     in kwargs:   sts.solve     = kwargs['solve']
    if 'step_p'    in kwargs:   sts.step_p    = kwargs['step_p']
    if 'vent_err'  in kwargs:   sts.vent_err  = kwargs['vent_err']
    if 'step_t'    in kwargs:   sts.step_t    = kwargs['step_t']
    if 'thrm_err'  in kwargs:   sts.thrm_err  = kwargs['thrm_err']
    if 'conv_err'  in kwargs:   sts.conv_err  = kwargs['conv_err']
    if 'sor_ratio' in kwargs:   sts.sor_ratio = kwargs['sor_ratio']
    if 'sor_err'   in kwargs:   sts.sor_err   = kwargs['sor_err']

    calc.setup(sts)

def set_sim_node(sn):
    sn_c = []
    for i, n in enumerate(sn):                                                      #sn
        calc.set_node(n["name"], i)                                                 #ノード番号
        v_flag = n['v_flag'] if 'v_flag' in n else SN_NONE
        c_flag = n['c_flag'] if 'c_flag' in n else SN_NONE
        t_flag = n['t_flag'] if 't_flag' in n else SN_NONE                          #計算フラグ

        calc.sn_add(i, [v_flag, c_flag, t_flag])
        sn_c.append(n['name'])

        if 'p'     in n:    calc.sn[i].p     = to_list_f(n['p'])                    #圧力、行列で設定可能
        if 'c'     in n:    calc.sn[i].c     = to_list_f(n['c'])                    #濃度、行列で設定可能
        if 't'     in n:    calc.sn[i].t     = to_list_f(n['t'])                    #温度、行列で設定可能
        if 'h_sr'  in n:    calc.sn[i].h_sr  = to_list_f(n['h_sr'])                 #日射量、行列で設定可能
        if 'h_inp' in n:    calc.sn[i].h_inp = to_list_f(n['h_inp'])                #発熱、行列で設定可能
        if 'v'     in n:    calc.sn[i].v     = to_list_f(n['v'])                    #気積、行列で設定可能
        if 'm'     in n:    calc.sn[i].m     = to_list_f(n['m'])                    #発生量、行列で設定可能
        if 'beta'  in n:    calc.sn[i].beta  = to_list_f(n['beta'])                 #濃度減少率、行列で設定可能

    return sn_c

def set_vent_net(vn):
    vn_c = []
    for i, nt in enumerate(vn):
        h1 = to_list_f(nt['h1']) if 'h1' in nt else to_list_f(0.0)                  #高さ1、行列設定不可
        h2 = to_list_f(nt['h2']) if 'h2' in nt else to_list_f(0.0)                  #高さ2、行列設定不可
        vn_type = nt['type'] if 'type' in nt else VN_FIX

        calc.vn_add(i, calc.node[nt['name1']], calc.node[nt['name2']], vn_type, h1, h2)
        vn_c.append(str(i) + " " + nt['name1'] + "->" + nt['name2'])

        if (vn_type == VN_FIX) or (vn_type == VN_AIRCON):       
            calc.vn[i].qv = to_list_f(nt['vol'])                                    #風量固定値、行列で設定可能
        if vn_type == VN_SIMPLE:    
            calc.vn[i].alpha = to_list_f(nt['alpha'])
            calc.vn[i].area  = to_list_f(nt['area'])                                #単純開口、行列で設定可能
        if vn_type == VN_GAP:           
            calc.vn[i].a     = to_list_f(nt['a'])
            calc.vn[i].n     = to_list_f(nt['n'])                                   #隙間、行列で設定可能
        if vn_type == VN_FAN:           
            calc.vn[i].q_max = to_list_f(nt['qmax']) 
            calc.vn[i].p_max = to_list_f(nt['pmax']) 
            calc.vn[i].q1   = to_list_f(nt['q1'])
            calc.vn[i].p1   = to_list_f(nt['p1'])                                   #ファン、行列で設定可能
        if vn_type == VN_AIRCON:
            calc.i_vn_ac = i
        calc.vn[i].eta = to_list_f(nt['eta']) if 'eta' in nt else to_list_f(0.0)    #粉じん除去率、行列で設定可能
    return vn_c

def set_thrm_net(sn, tn, sn_c):
    tn_c = []
    for i, nt in enumerate(tn):                                                     #tn
        tn_type = nt['type'] if 'type' in nt else TN_SIMPLE

        calc.tn_add(i, calc.node[nt['name1']], calc.node[nt['name2']], tn_type)
        tn_c.append(str(i) + " " + nt['name1'] + "->" + nt['name2'])

        if tn_type == TN_SIMPLE:     
            calc.tn[i].cdtc = to_list_f(nt['cdtc'])                                 #コンダクタンス、行列設定可能
        if tn_type == TN_AIRCON:     
            calc.tn[i].ac_mode = to_list_i(nt['ac_mode']) 
            calc.tn[i].pre_tmp = to_list_f(nt['pre_tmp'])                           #エアコン運転モード
            calc.i_tn_ac = i
        if tn_type == TN_SOLAR:       
            calc.tn[i].ms      = to_list_f(nt['ms'])                                #日射熱取得率、行列設定可能
        if tn_type == TN_GROUND:     
            calc.tn[i].area    = to_list_f(nt['area'])           
            calc.tn[i].rg      = to_list_f(nt['rg']) 
            calc.tn[i].phi_0   = nt['phi_0']
            calc.tn[i].cof_r   = nt['cof_r']
            calc.tn[i].cof_phi = nt['cof_phi']                                      #地盤熱応答、行列設定不可（面積と断熱性能はOK）

    print('Add Capacity.')
    for i, n in enumerate([n for n in sn if 'capa' in n]):                                              #熱容量の設定のあるノード
        calc.set_node(d_node(n['name']), len(sn) + i)                                                   #時間遅れノードのノード番号
        sn_c.append(d_node(n['name']))

        calc.sn_add(len(sn) + i, [SN_NONE, SN_NONE, SN_DLY])                                            #計算フラグ
        calc.sn[len(sn) + i].s_i = calc.node[n['name']]
        if 't' in n:    calc.sn[len(sn) + i].t = to_list_f(n['t'])

        calc.tn_add(len(tn) + i, calc.node[n['name']], calc.node[d_node(n['name'])], TN_SIMPLE)         #熱容量の設定
        tn_c.append(str(len(tn) + i) + " " + n['name'] + "->" + d_node(n['name']))
        calc.tn[len(tn) + i].cdtc = to_list_f(n['capa'] / calc.sts.t_step)                              #コンダクタンス（熱容量）            

    return sn_c, tn_c

def output_calc(ix, opt, res, sn_c, vn_c, tn_c):
    print('Create pd.DataFrames')
    
    dat_list  = [{'df': pd.DataFrame(), 'columns': sn_c, 'fn': 'vent_p.csv',   'title': '圧力',  'unit': '[Pa]'},
                 {'df': pd.DataFrame(), 'columns': sn_c, 'fn': 'vent_c.csv',   'title': '濃度',  'unit': '[個/L]'},
                 {'df': pd.DataFrame(), 'columns': sn_c, 'fn': 'them_t.csv',   'title': '温度',  'unit': '[℃]'},
                 {'df': pd.DataFrame(), 'columns': vn_c, 'fn': 'vent_qv.csv',  'title': '風量',  'unit': '[m3/s]'},
                 {'df': pd.DataFrame(), 'columns': vn_c, 'fn': 'thrm_qt1.csv', 'title': '熱量1', 'unit': '[W]'},
                 {'df': pd.DataFrame(), 'columns': tn_c, 'fn': 'thrm_qt2.csv', 'title': '熱量2', 'unit': '[W]'}]
    
    for i, d in enumerate(dat_list):
        if len(res[i]) != 0: d['df'] = pd.DataFrame(np.array(res[i]).T,  index = ix, columns = d['columns'])

    if opt >= OPT_CSV:
        print('Output csv files.')
        for d in dat_list:      d['df'].to_csv(d['fn'], encoding = 'utf_8_sig')

    if opt >= OPT_GRAPH:
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

    global df_p, df_c, df_t, df_qv, df_qt1, df_qt2
    df_p   = dat_list[0]['df']
    df_c   = dat_list[1]['df']
    df_t   = dat_list[2]['df'] 
    df_qv  = dat_list[3]['df']
    df_qt1 = dat_list[4]['df']
    df_qt2 = dat_list[5]['df']

    return dat_list[0]['df'], dat_list[1]['df'], dat_list[2]['df'], dat_list[3]['df'], dat_list[4]['df'], dat_list[5]['df']