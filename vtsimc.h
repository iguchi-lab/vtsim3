#include "vtsim_set.h"              //計算の設定　ヘッダーファイルの読み込み
#include "node.h"                   //ノード　ヘッダーファイルの読み込み
#include "vent_net.h"               //換気回路網　ヘッダーファイルの読み込み
#include "thrm_net.h"               //熱回路網　ヘッダーファイルの読み込み
#include "mymath.h"
#include <iostream>
#include <fstream>

//#define DEBUG_ON

#ifdef  DEBUG_ON
#define LOG_PRINT(...)     ofs << __FILE__ << " (" << __LINE__ << ") " << __func__ << ":" << __VA_ARGS__
#define LOG_CONTENTS(...)  ofs << __VA_ARGS__
const char *logfileName = "log.txt";
ofstream ofs(logfileName);
#else
#define LOG_PRINT(...)
#define LOG_CONTENTS(...)
#endif

using namespace::std;

class CalcStatus{
public:
    long length      = 0;
    double t_step    = 3600;
    int solve        = SOLVE_LU;
    double step_p    = STEP_P, 
           vent_err  = VENT_ERR, 
           step_t    = STEP_T, 
           thrm_err  = THRM_ERR, 
           conv_err  = CONV_ERR, 
           sor_ratio = SOR_RATIO, 
           sor_err   = SOR_ERR;
};

class VTSim{
public:
    CalcStatus sts;
    vector<Node> sn;                                        //ノード
    vector<Vent_Net> vn;                                    //換気回路網
    vector<Thrm_Net> tn;                                    //熱回路網
    vector<int> v_idc, c_idc, t_idc, ac_idc;
    int i_vn_ac = -1, i_tn_ac = -1;

    void setup(sts_){
        sts = sts_;
        sn.clear();
        vn.clear();
        tn.clear();
        v_idc.clear();
        c_idc.clear();
        t_idc.clear();
        ac_idc.clear();
        i_vn_ac = -1;
        i_tn_ac - -1;
    }

    void sn_add(int i, tuple<int, int, int> flag){
        sn.push_back(Node(sts.length, i, flag));
    }

    void vn_add(int i, int i1, int i2, int vn_type, vector<double> h1, vector<double> h2){
        vn.push_back(Vent_Net(sts.length, i, i1, i2, vn_type, h1, h2));
    }

    void tn_add(int i, int i1, int i2, int tn_type){
        tn.push_back(Thrm_Net(sts.length, i, i1, i2, tn_type));
    }

    void change_sn_t_flag(int i, int flag_){
        get<2>(sn[i].flag) = flag_;
        t_idc.clear();
        for(int i = 0; i < sn.size(); i++)   
            if(get<2>(sn[i].flag) == SN_CALC)     t_idc.push_back(sn[i].i);
    }

    void set_aircon_status(long ts, int flag_){
        if(tn[i_tn_ac].aircon_on[ts] == AC_OFF && flag_ == AC_ON){
            tn[i_tn_ac].aircon_on[ts] = AC_ON;
            change_sn_t_flag(tn[i_tn_ac].i1, SN_FIX);
            vn[i_vn_ac].vn_type = VN_AIRCON;
        }
        else if(tn[i_tn_ac].aircon_on[ts] == AC_ON && flag_ == AC_OFF){
            tn[i_tn_ac].aircon_on[ts] = AC_OFF;
            change_sn_t_flag(tn[i_tn_ac].i1, SN_CALC);  
            vn[i_vn_ac].vn_type = VN_FIX;
        }
    }

    vector<double> qv_sum(vector<double> p, long ts){
        vector<double> qvsum(sn.size(), 0.0);                                                                             //風量収支の初期化                                 
        for(int i = 0; i < vn.size(); i++){
            double rgh1 = get_rho(sn[vn[i].i1].t[ts]) * G * vn[i].h1[ts];
            double rgh2 = get_rho(sn[vn[i].i2].t[ts]) * G * vn[i].h2[ts];
            double qv   = vn[i].get_qv((p[vn[i].i1] - rgh1) - (p[vn[i].i2] - rgh2), ts);
            qvsum[vn[i].i1] -= qv;                                                                                        //風量収支の計算
            qvsum[vn[i].i2] += qv;                                                                                        //風量収支の計算      
        }
        return qvsum;
    }

    double calc_vent(long ts){
        vector<double>          p0(sn.size()), qvsum_0(sn.size()), qvsum_d(sn.size());
        vector<vector<double>>  a(v_idc.size(), vector<double>(v_idc.size()));
        vector<double>          b(v_idc.size()), dp(v_idc.size());
        double                  rmse;

        for(int i = 0; i < sn.size(); i++)    p0[sn[i].i] = sn[i].p[ts];                                                //圧力の初期化       
        do{
            rmse = 0.0;
            qvsum_0 = qv_sum(p0, ts);
            
            for(int j = 0; j < v_idc.size(); j++){
                p0[v_idc[j]] += sts.step_p;                                                                             //ダミー圧力の作成       
                qvsum_d = qv_sum(p0, ts);                                                                               //ダ三－風量収支の計算
                for(int i = 0; i < v_idc.size(); i++)  a[i][j] = (qvsum_d[v_idc[i]] - qvsum_0[v_idc[i]]) / sts.step_p;  //aの計算
                b[j] = -qvsum_0[v_idc[j]];
                p0[v_idc[j]] -= sts.step_p; 
            }

            if(sts.solve == SOLVE_SOR)  dp = SOR(a, b, v_idc.size(), sts.sor_ratio, sts.sor_err);                                   //SOR法による計算     
            else                        dp = LU(a, b, v_idc.size()); 

            for(int i = 0; i < v_idc.size(); i++)   p0[v_idc[i]] += dp[i];                                              //圧力の更新
            qvsum_0 = qv_sum(p0, ts);    
            for(int i = 0; i < v_idc.size(); i++)   rmse += pow(qvsum_0[v_idc[i]], 2.0) / v_idc.size();

        }while(sts.vent_err < sqrt(rmse));
        for(int i = 0; i < v_idc.size(); i++)   sn[v_idc[i]].p[ts] = p0[v_idc[i]];                                      //圧力の計算                                                               
        return rmse;
    }

    vector<double> qt_sum(vector<double> t, long ts){
        vector<double> qtsum(sn.size(), 0.0);           
       
        for(int i = 0; i < vn.size(); i++){                                                                             //移流に伴う熱移動
            if(vn[i].qv[ts] > 0)            qtsum[vn[i].i2] += vn[i].get_qt(t[vn[i].i1] - t[vn[i].i2], ts);
            else                            qtsum[vn[i].i1] += vn[i].get_qt(t[vn[i].i1] - t[vn[i].i2], ts);
        }

        for(int i = 0; i < tn.size(); i++){                                                                             //貫流、日射、発熱による熱移動
            switch(tn[i].tn_type){
                case TN_SIMPLE:
                case TN_GROUND:
                    qtsum[tn[i].i1] -= tn[i].get_qt(t[tn[i].i1] - t[tn[i].i2], ts);     
                    qtsum[tn[i].i2] += tn[i].get_qt(t[tn[i].i1] - t[tn[i].i2], ts);
                    break;
                case TN_SOLAR:
                    qtsum[tn[i].i1] += tn[i].get_qt_s(sn[tn[i].i2].h_sr[ts], ts);
                    break;
                case TN_HEATER:
                    qtsum[tn[i].i1] += tn[i].get_qt_h(sn[tn[i].i2].h_inp[ts], ts);
                    break;
            }    
        }

        if(i_tn_ac != -1 && tn[i_tn_ac].aircon_on[ts] == AC_ON){                                                        //エアコンによる温度調整
            qtsum[tn[i_tn_ac].i2] -= qtsum[tn[i_tn_ac].i1];
            qtsum[tn[i_tn_ac].i1]  = 0.0;
        }
        return qtsum;
    }

    double calc_thrm(long ts){
        vector<double>          t0(sn.size()), qtsum_0(sn.size()), qtsum_d(sn.size());
        vector<vector<double>>  a(t_idc.size() + 1, vector<double>(t_idc.size() + 1));                                  //エアコン分も余計に確保
        vector<double>          b(t_idc.size() + 1), dt(t_idc.size() + 1);                                              //エアコン分も余計に確保
        double                  rmse;        
        int                     itr = 0;

        for(int i = 0; i < sn.size(); i++)    t0[sn[i].i] = sn[i].t[ts];                                                //温度の初期化
        
        do{
            if(i_tn_ac != -1){
                set_aircon_status(ts, AC_ON);                                                                           //全部ONになってない？？
                if(itr == 0)    set_aircon_status(ts, AC_OFF);
                else{
                    switch(tn[i_tn_ac].ac_mode[ts]){
                        case AC_AUTO:       set_aircon_status(ts, AC_ON);
                                            break;
                        case AC_COOLING:    if(t0[tn[i_tn_ac].i1] >= tn[i_tn_ac].pre_tmp[ts])   set_aircon_status(ts, AC_ON);
                                            break;
                        case AC_HEATING:    if(t0[tn[i_tn_ac].i1] <= tn[i_tn_ac].pre_tmp[ts])   set_aircon_status(ts, AC_ON);
                                            break;
                    }
                }
            }

            rmse = 0.0;
            qtsum_0 = qt_sum(t0, ts);                                                                                   //熱量収支の計算

            for(int j = 0; j < t_idc.size(); j++){
                t0[t_idc[j]] += sts.step_t;                                                                             //ダミー温度の作成
                qtsum_d = qt_sum(t0, ts);                                                                               //ダミー熱量の計算
                for(int i = 0; i < t_idc.size(); i++)  a[i][j] = (qtsum_d[t_idc[i]] - qtsum_0[t_idc[i]]) / sts.step_t;  //aの計算
                b[j] = -qtsum_0[t_idc[j]];                                                                              //bの計算
                t0[t_idc[j]] -= sts.step_t;                                                                             //ダミー温度を戻す
            }

            if(sts.solve == SOLVE_SOR)  dt = SOR(a, b, t_idc.size(), sts.sor_ratio, sts.sor_err);                                   //SOR法による計算
            else                        dt = LU(a, b, t_idc.size());

            for(int i = 0; i < t_idc.size(); i++)   t0[t_idc[i]] += dt[i];                                              //温度の更新
            qtsum_0 = qt_sum(t0, ts);
            for(int i = 0; i < t_idc.size(); i++)   rmse += pow(qtsum_0[t_idc[i]], 2.0) / t_idc.size(); 

            if(i_tn_ac != -1){
                switch(tn[i_tn_ac].ac_mode[ts]){
                    case AC_AUTO:       if(abs(t0[tn[i_tn_ac].i1] - tn[i_tn_ac].pre_tmp[ts]) > sts.thrm_err){
                                            t0[tn[i_tn_ac].i1] = tn[i_tn_ac].pre_tmp[ts];
                                            rmse = 999;
                                        }
                                        break;
                    case AC_HEATING:    if(t0[tn[i_tn_ac].i1] < tn[i_tn_ac].pre_tmp[ts]){
                                            t0[tn[i_tn_ac].i1] = tn[i_tn_ac].pre_tmp[ts];
                                            rmse = 999;
                                        }                   
                                        break;
                    case AC_COOLING:    if(t0[tn[i_tn_ac].i1] > tn[i_tn_ac].pre_tmp[ts]){
                                            t0[tn[i_tn_ac].i1] = tn[i_tn_ac].pre_tmp[ts];
                                            rmse = 999;
                                        }
                                        break;
                }       
            }
            itr++;
        }while(sts.thrm_err < sqrt(rmse));

        for(int i = 0; i < t_idc.size(); i++)   sn[t_idc[i]].t[ts] = t0[t_idc[i]];                                          //温度の計算
        return rmse;
    }

    void calc_qv(long ts1, long ts2){
        for(int i = 0; i < vn.size(); i++){    
            double rgh1 = get_rho(sn[vn[i].i1].t[ts2]) * G * vn[i].h1[ts2];
            double rgh2 = get_rho(sn[vn[i].i2].t[ts2]) * G * vn[i].h2[ts2];
            vn[i].qv[ts1] = vn[i].get_qv((sn[vn[i].i1].p[ts2] - rgh1) - (sn[vn[i].i2].p[ts2] - rgh2), ts2);                 //風量の計算
        }
    }

    void calc_qt(long ts1, long ts2){
        for(int i = 0; i < tn.size(); i++){    
            switch(tn[i].tn_type){
                case TN_SIMPLE: 
                case TN_GROUND:         tn[i].qt[ts1] =   tn[i].get_qt(sn[tn[i].i1].t[ts2] - sn[tn[i].i2].t[ts2], ts2);
                                        break;
                case TN_SOLAR:          tn[i].qt[ts1] = - tn[i].get_qt_s(sn[tn[i].i2].h_sr[ts2], ts2);
                                        break;
                case TN_HEATER:         tn[i].qt[ts1] = - tn[i].get_qt_h(sn[tn[i].i2].h_inp[ts2], ts2);
                                        break;
            }
        }
        for(int i = 0; i < vn.size(); i++)
            vn[i].qt[ts1] = vn[i].get_qt(sn[vn[i].i1].t[ts2] - sn[vn[i].i2].t[ts2], ts2);                                   //熱量の計算
    }

    int calc(){
        vector<double>  pre_p(sn.size(), 0.0), pre_t(sn.size(), 0.0);
        double          delta_p, delta_t; 

        LOG_PRINT("start calc");     

        for(long ts = 1; ts < sts.length; ts++){
            if(v_idc.size() > 0){
                for(int i = 0; i < sn.size(); i++)
                    if(get<0>(sn[i].flag) == SN_CALC)   sn[i].p[ts] = sn[i].p[ts - 1];
                calc_qv(ts, ts - 1);
            }

            if(t_idc.size() > 0){    
                for(int i = 0; i < sn.size(); i++){
                    if(get<2>(sn[i].flag) == SN_CALC)   sn[i].t[ts] = sn[i].t[ts - 1];
                    if(get<2>(sn[i].flag) == SN_DLY)    sn[i].t[ts] = sn[sn[i].s_i].t[ts - 1];
                }
                calc_qt(ts, ts - 1);
                for(int i = 0; i < tn.size(); i++)
                    if(tn[i].tn_type == TN_GROUND)      tn[i].refresh(sn[tn[i].i1].t[ts], sn[tn[i].i2].t[ts], ts);   //地盤のリフレッシュ
            }

            do{
                delta_p = 0.0;
                delta_t = 0.0;

                for(int i = 0; i < sn.size(); i++){
                    if(get<0>(sn[i].flag) == SN_CALC)   pre_p[i] = sn[i].p[ts];
                    if(get<2>(sn[i].flag) == SN_CALC)   pre_t[i] = sn[i].t[ts];
                }
                LOG_PRINT("ts = " << ts << " calc_vent" << endl);
                calc_vent(ts);
                LOG_PRINT("ts = " << ts << " calc_thrm" << endl);
                calc_thrm(ts);
                
                for(int i = 0; i < sn.size(); i++){
                    if(get<0>(sn[i].flag) == SN_CALC)   delta_p += pow(sn[i].p[ts] - pre_p[i], 2.0) / v_idc.size();
                    if(get<2>(sn[i].flag) == SN_CALC)   delta_t += pow(sn[i].t[ts] - pre_t[i], 2.0) / t_idc.size();
                }
                LOG_PRINT("ts = " << ts << " delta_p = " << delta_p << ", delta_t = " << delta_t << " >>>>> " << sqrt((delta_p + delta_t)/ 2) << endl);

            }while(sts.conv_err < sqrt((delta_p + delta_t)/ 2));
            
            if(v_idc.size() > 0){
                calc_qv(ts, ts);
                LOG_CONTENTS(endl << "p ");
                for(int i = 0; i < sn.size(); i++)    LOG_CONTENTS(sn[i].i << ": " << sn[i].p[ts] << "Pa, ");
                LOG_CONTENTS("qv ");
                for(int i = 0; i < vn.size(); i++)    LOG_CONTENTS(vn[i].i << ": " << vn[i].qv[ts] * 3600 << "m3/h, ");
                LOG_CONTENTS(endl << endl);    
            }

            if(t_idc.size()> 0){
                calc_qt(ts, ts);
                LOG_CONTENTS(endl << "t ");
                for(int i = 0; i < sn.size(); i++)    LOG_CONTENTS(sn[i].i << ": " << sn[i].t[ts] << "deg, ");
                LOG_CONTENTS(endl << "qt1 ");
                for(int i = 0; i < vn.size(); i++)    LOG_CONTENTS(vn[i].i << ": " << vn[i].qt[ts] << "W, ");
                LOG_CONTENTS(endl << "qt2 ");
                for(int i = 0; i < tn.size(); i++)    LOG_CONTENTS(tn[i].i << ": " << tn[i].qt[ts] << "W, ");
                LOG_CONTENTS(endl << endl);
            }

            if(c_idc.size() > 0){
                for(int i = 0; i < c_idc.size(); i++){
                    sn[c_idc[i]].c[ts] =  sn[c_idc[i]].c[ts - 1] * exp(-sn[c_idc[i]].beta[ts] * sts.t_step);
                    sn[c_idc[i]].c[ts] += sn[c_idc[i]].m[ts] / (sn[c_idc[i]].beta[ts] * sn[c_idc[i]].v[ts]) * (1 - exp(-sn[c_idc[i]].beta[ts] * sts.t_step));
                }
                for(int i = 0; i < vn.size(); i++){
                    if(vn[i].qv[ts] > 0 && get<1>(sn[vn[i].i2].flag) == SN_CALC)    
                        sn[vn[i].i2].c[ts] += vn[i].qv[ts] * (sn[vn[i].i1].c[ts] * (1 - vn[i].eta[ts]) - sn[vn[i].i2].c[ts]) * sts.t_step / sn[vn[i].i2].v[ts];
                    if(vn[i].qv[ts] < 0 && get<1>(sn[vn[i].i1].flag) == SN_CALC)    
                        sn[vn[i].i1].c[ts] += vn[i].qv[ts] * (sn[vn[i].i2].c[ts] * (1 - vn[i].eta[ts]) - sn[vn[i].i1].c[ts]) * sts.t_step / sn[vn[i].i1].v[ts];
                }
                for(int i = 0; i < sn.size(); i++)    LOG_CONTENTS("c" << sn[i].i << ": " << sn[i].c[ts] << "num/L, ");
                LOG_CONTENTS(endl << endl);
            }
        }
        return 0;
    }

    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, 
          vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> result(){

        vector<vector<double>>   r_p, r_c, r_t, r_qv, r_qt1, r_qt2;
    
        for(Node n: sn){
            r_p.push_back(n.p);
            r_c.push_back(n.c);
            r_t.push_back(n.t);
        }
        for(Vent_Net nt: vn){
            r_qv.push_back(nt.qv);
            r_qt1.push_back(nt.qt);
        }
        for(Thrm_Net nt: tn)   
            r_qt2.push_back(nt.qt);

        return forward_as_tuple(r_p, r_c, r_t, r_qv, r_qt1, r_qt2);
    }

};
