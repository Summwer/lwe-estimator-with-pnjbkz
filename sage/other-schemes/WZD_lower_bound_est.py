import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))
from Lower_Bound_Estimation import *



#WZD param
q = 2**38
n = 1792
m = n
D_s = {x : 1./2 for x in range(0, 2)}
D_e =  {x : 1./2 for x in range(0, 2)}
d_max = m+n+1
mu_s, s_s = average_variance(D_s)
mu_e, s_e = average_variance(D_e)
print("-----WZD")
sigma = sqrt(s_e)
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))

#Lower_Bound_Estimation
print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)
