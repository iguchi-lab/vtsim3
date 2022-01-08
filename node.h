using namespace std;

class Node{
public:
    int i, s_i;
    tuple<int, int, int> flag;                              //換気、濃度、熱計算フラグ                
    vector<double> p, c, t;                                 //圧力、濃度、温度
    vector<double> m;                                       //発生量
    vector<double> h_sr, h_inp;                             //日射量、発熱量
    vector<double> v;                                       //気積     
    vector<double> beta;                                    //沈着率

    Node(long length, int i_, tuple<int, int, int> SN_){
        i = i_;
        flag = SN_;
        p.assign(length, 0.0);
        c.assign(length, 0.0);
        t.assign(length, 20.0);
        if(get<1>(flag) == SN_CALC)   m.assign(length, 0.0);
        if(get<2>(flag) == SN_CALC)   h_inp.assign(length, 0.0);
    }

    void set_P(vector<double> p_)           {p     = p_;}    //圧力固定値
    void set_C(vector<double> c_)           {c     = c_;}    //濃度固定値
    void set_T(vector<double> t_)           {t     = t_;}    //温度固定値
    void set_h_sr(vector<double> h_sr_)     {h_sr  = h_sr_;}
    void set_h_inp(vector<double> h_inp_)   {h_inp = h_inp_;}
    void set_v(vector<double> v_)           {v     = v_;}
    void set_capa(int s_i_)                 {s_i = s_i_;}
    void set_m(vector<double> m_)           {m     = m_;}
    void set_beta(vector<double> beta_)     {beta  = beta_;}
};