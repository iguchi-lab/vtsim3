import numpy as np

# case0：省エネ法の計算で使用される係数
PHI_0_0 = 0.0255                                                                                                        #吸熱応答係数0      0.0255
COF_r_0 =   np.array([ 0.999996185, 0.99998474,   0.999938962, 0.99975593,  0.999023877, 
                       0.996101618, 0.984491515,  0.939410628, 0.778800783, 0.367879441])                               #公比
COF_phi_0 = np.array([-3.73933E-08, 7.33578E-07, -9.94021E-06, 0.000584299, 0.000384435, 
                       0.000618679, 0.002742905,  0.00087726,  0.011610039, 0.001565569])                               #φ
                                                                                                #吸熱応答係数0      0.0255

# case1：断熱材50mm + コンクリート120mm ＋ 地盤3000mm
PHI_0_1 = 1.713017E+00  
COF_r_1 =   np.array([0.999999619, 0.999998474, 0.999993897, 0.999975586, 0.999902349, 
                      0.999609451, 0.99843872,  0.993769491, 0.975309912, 0.904837418])                                 #公比
COF_phi_1 = np.array([-9.784586E-10, 2.722309E-08, -5.966737E-07, 5.476697E-05, 5.780554E-05, 
                      -2.704139E-05, 6.676871E-04, -1.508605E-03, 7.168433E-03, -1.244022E-02])                         #φ
                      
PHI_0_1_10min   = 1.71301691744404
COF_r_1_10min   = np.array([0.999999619, 0.999998474, 0.999993897, 0.999975586, 0.999902349,
                            0.999609451, 0.99843872,  0.993769491, 0.975309912, 0.904837418])
COF_phi_1_10min = np.array([-9.784586E-10, 2.722309E-08, -5.966737E-07, 5.476697E-05,  5.780554E-05,
                            -2.704139E-05, 6.676871E-04, -1.508605E-03, 7.168433E-03, -1.244022E-02])

# case2：コンクリート120mm + 断熱材50mm + 地盤3000mm
PHI_0_2 = 2.934455E-02                                                                                                  #吸熱応答係数0      0.0255
COF_r_2 =   np.array([ 0.999999921, 0.999999682, 0.999998728, 0.999994914, 0.999979655, 
                       0.999918623, 0.999674532, 0.998698764, 0.994805207, 0.979382181])                                #公比
COF_phi_2 = np.array([-4.068474E-12, 1.085816E-10, -1.975424E-09, 3.304520E-08, -5.629434E-07,
                       1.167294E-05, 9.036604E-04, -3.759300E-04, 9.467514E-03, 5.798151E-03])                          #φ

PHI_0_2_10min   = 0.0221839067807554
COF_r_2_10min   = np.array([0.999999992, 0.999999968, 0.999999873, 0.999999491, 0.999997965,
                            0.999991862, 0.999967448, 0.9998698,   0.999479302, 0.997918835])
COF_phi_2_10min = np.array([-4.068474E-13, 1.085816E-11, -1.975426E-10, 3.304535E-09, -5.629537E-08,
                             1.167380E-06, 9.039252E-05, -3.763708E-05, 9.511976E-04,  5.907680E-04])

#case3：コンクリート120mm + 地盤3000mm
PHI_0_3 = 1.263304E-01                                                                                                  #吸熱応答係数0      0.0255
COF_r_3 =   np.array([ 0.999999921, 0.999999682, 0.999998728, 0.999994914, 0.999979655, 
                       0.999918623, 0.999674532, 0.998698764, 0.994805207, 0.979382181])                                #公比
COF_phi_3 = np.array([-5.612334E-13,  1.496760E-11, -2.715127E-10, 4.488739E-09, -7.281398E-08,
                       1.195868E-06, -2.125076E-05,  6.342970E-04, 3.113842E-03, -1.747075E-03])                        #φ

PHI_0_3_10min   = 0.125444510728293                                                                                     #吸熱応答係数0      0.0255
COF_r_3_10min   = np.array([0.999999992, 0.999999968, 0.999999873, 0.999999491, 0.999997965,
                            0.999991862, 0.999967448, 0.9998698,   0.999479302, 0.997918835])                           #公比
COF_phi_3_10min = np.array([-5.612334E-14,  1.496761E-12, -2.715130E-11, 4.488759E-10, -7.281531E-09,
                             1.195956E-07, -2.125699E-06,  6.350407E-05, 3.128466E-04, -1.780078E-04])                  #φ