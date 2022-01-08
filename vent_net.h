using namespace std;

class Vent_Net{
public:
    int i, i1, i2, vn_type;                                         //入口出口のノード番号、換気回路網の種類
    vector<double> h1, h2;                                          //入口出口の高さ
    vector<double> alpha, area;                                     //単純開口　開口率、面積
    vector<double> a, n;                                            //隙間　　　隙間量、隙間特性値
    vector<double> qv, qt;                                          //風量、移流に伴う熱量
    vector<double> eta;                                             //粉塵除去率
    vector<double> q_max, p_max, q1, p1;                            //送風ファンの風量・圧力

    Vent_Net(long length, int i_, tuple<int, int, int, vector<double>, vector<double>> v_nets){
        i       = i_;
        forward_as_tuple(i1, i2, vn_type, h1, h2) = v_nets; 
        qv.assign(length, 0.0);                                     //風量を0に初期化
        qt.assign(length, 0.0);                                     //熱量を0に初期化
    }

    void set_Simple(vector<double> alpha_, vector<double> area_){
        alpha = alpha_;                                             //開口率
        area  = area_;                                              //面積
    }

    void set_Gap(vector<double> a_, vector<double> n_){
        a = a_;                                                     //隙間量
        n = n_;                                                     //隙間特性値
    }

    void set_Fix(vector<double> vol_){                                             
        qv = vol_;                                                  //風量固定値
    }

    void set_Fan(vector<double> q_max_, vector<double> p_max_, 
                 vector<double> q1_,    vector<double> p1_){
        q_max = q_max_;                                             //最大風量
        p_max = p_max_;                                             //最大静圧
        q1 = q1_;                                                   //中間の風量
        p1 = p1_;                                                   //中間の静圧
    }

    void set_Eta(vector<double> eta_){
        eta = eta_;
    }

    double get_qv(double dp, long ts){
        switch(vn_type){
            case VN_SIMPLE:                                                         //単純開口
                if(dp >= 0)                                   return  alpha[ts] * area[ts] * sqrt( 2.0 * dp / RHO20);
                else                                          return -alpha[ts] * area[ts] * sqrt(-2.0 * dp / RHO20);
            case VN_GAP:                                                            //隙間
                if(dp >= 0)                                   return  a[ts] * pow( dp, 1 / n[ts]);
                else                                          return -a[ts] * pow(-dp, 1 / n[ts]);
            case VN_FIX:
            case VN_AIRCON:                                   return qv[ts];      //固定・エアコン 
            case VN_FAN:                                                          //送風ファン
                if (-dp <= 0.0)                               return q_max[ts];
                else if((0 < -dp) && (-dp <= p1[ts]))         return q_max[ts] + (q1[ts] - q_max[ts]) / p1[ts] * (-dp);
                else if((p1[ts] < -dp) && (-dp <= p_max[ts])) return q1[ts] - q1[ts] * (p_max[ts] - p1[ts]) * (-dp - q1[ts]);
                else                                          return 0.0;
            default:                                          return 0.0;
        }
    }

    double get_qt(double dt, long ts){
        switch(vn_type){
            case VN_AIRCON:                                 return 0.0;           //エアコンの場合のみ0
            default:                                        return RHO20 * AIR_CP * qv[ts] * dt;
        }
    }
};