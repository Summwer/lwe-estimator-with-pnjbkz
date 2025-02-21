import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))
from Lower_Bound_Estimation import *


#Kyber-I
n = 512
m = 512
q = 3329
d_max = m+n+1
D_s = build_centered_binomial_law(3)
D_e = D_s
mu_s, s_s = average_variance(D_s)
mu_e, s_e = average_variance(D_e)
print("-----Kyber-I")
sigma = sqrt(s_e)
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))

#Core-SVP
#m_op,beta_op=Core_SVP(d_max,n,q,sigma)
#print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )


#Lower_Bound_Estimation
print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)


#Lower_Bound_Estimation_Optimal_m
print("\n Lower_Bound_Estimation_Optimal_m ")
d_op,beta_op,d_svp_op,T_LBE_op = Lower_Bound_Estimation_Optimal_d(m,n,sigma,q,m-26)
print("Lower_Bound_Estimation_Optimal_Dimension = %d" % d_op)
print("Lower_Bound_Estimation_Optimal_m: Current Beta = %d" % beta_op)
print("Lower_Bound_Estimation_Optimal_m: d_svp = %d" % d_svp_op)
print("Lower_Bound_Estimation_Optimal_m: Hardness = %f Bits" % T_LBE_op)


#Kyber-II
n = 768
m = 768
q = 3329
d_max = m+n+1
D_s = build_centered_binomial_law(2)
D_e = D_s
mu_s, s_s = average_variance(D_s)
mu_e, s_e = average_variance(D_e)
print("-----Kyber-II")
sigma = sqrt(s_e)
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))


#Core-SVP
#m_op,beta_op=Core_SVP(d_max,n,q,sigma)
#print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )

#Lower_Bound_Estimation
print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)


#Lower_Bound_Estimation_Optimal_m
print("\n \n Lower_Bound_Estimation_Optimal_m ")
d_op,beta_op,d_svp_op,T_LBE_op = Lower_Bound_Estimation_Optimal_d(m,n,sigma,q,m-79)
print("Lower_Bound_Estimation_Optimal_Dimension = %d" % d_op)
print("Lower_Bound_Estimation_Optimal_m: Current Beta = %d" % beta_op)
print("Lower_Bound_Estimation_Optimal_m: d_svp = %d" % d_svp_op)
print("Lower_Bound_Estimation_Optimal_m: Hardness = %f Bits" % T_LBE_op)




#Kyber-III
n = 1024
m = 1024
q = 3329
d_max = m+n+1
D_s = build_centered_binomial_law(2)
D_e = D_s
mu_s, s_s = average_variance(D_s)
mu_e, s_e = average_variance(D_e)
print("-----Kyber-III")
sigma = sqrt(s_e)
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))


#Core-SVP
#m_op,beta_op=Core_SVP(d_max,n,q,sigma)
#print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )

print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)


#Lower_Bound_Estimation_Optimal_m
print("\n \n Lower_Bound_Estimation_Optimal_m ")
d_op,beta_op,d_svp_op,T_LBE_op = Lower_Bound_Estimation_Optimal_d(m,n,sigma,q,m-182)
print("Lower_Bound_Estimation_Optimal_Dimension = %d" % d_op)
print("Lower_Bound_Estimation_Optimal_m: Current Beta = %d" % beta_op)
print("Lower_Bound_Estimation_Optimal_m: d_svp = %d" % d_svp_op)
print("Lower_Bound_Estimation_Optimal_m: Hardness = %f Bits" % T_LBE_op)




#Dilithium-III
#N= 256
#l = 7
#n = l*N
#k = 8
#m = k*N
#q = 8380417
#eta_1 = 2
#eta_2 = 2
#d_max = m+n+1
#sigma = eta_1/sqrt(3)
print("-----Dilithium-III")
n = 7*256
m = 8*256
d_max = m+n+1
q = 8380417
eta = 2
D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)}
D_e = D_s
mu_e, s_e = average_variance(D_e)
mu_s, s_s = average_variance(D_s)
sigma = sqrt(s_e)
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))


#Core-SVP
#m_op,beta_op=Core_SVP(d_max,n,q,sigma)
#print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )


#Lower_Bound_Estimation
print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)

#Lower_Bound_Estimation_Optimal_m
print("\n Lower_Bound_Estimation_Optimal_m ")
d_op,beta_op,d_svp_op,T_LBE_op = Lower_Bound_Estimation_Optimal_d(m,n,sigma,q,m-408)
print("Lower_Bound_Estimation_Optimal_Dimension = %d" % d_op)
print("Lower_Bound_Estimation_Optimal_m: Current Beta = %d" % beta_op)
print("Lower_Bound_Estimation_Optimal_m: d_svp = %d" % d_svp_op)
print("Lower_Bound_Estimation_Optimal_m: Hardness = %f Bits" % T_LBE_op)



        
#Dilithium-I
n = 4*256
m = 4*256
d_max = m+n+1
q = 8380417
eta = 2
D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)}
D_e = D_s
mu_e, s_e = average_variance(D_e)
mu_s, s_s = average_variance(D_s)
print("-----Dilithium-I")
sigma = sqrt(s_e)
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))

#Core-SVP
#m_op,beta_op=Core_SVP(d_max,n,q,sigma)
#print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )


#Lower_Bound_Estimation
print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)

#Lower_Bound_Estimation_Optimal_m
print("\n Lower_Bound_Estimation_Optimal_m ")
d_op,beta_op,d_svp_op,T_LBE_op = Lower_Bound_Estimation_Optimal_d(m,n,sigma,q,m-73)
print("Lower_Bound_Estimation_Optimal_Dimension = %d" % d_op)
print("Lower_Bound_Estimation_Optimal_m: Current Beta = %d" % beta_op)
print("Lower_Bound_Estimation_Optimal_m: d_svp = %d" % d_svp_op)
print("Lower_Bound_Estimation_Optimal_m: Hardness = %f Bits" % T_LBE_op)



#Dilithium-II
#N= 256
#l = 5
#n = l*N
#k = 6
#m = k*N
#q = 8380417
#eta_1 = 4
#eta_2 = 4

#sigma = eta_1/sqrt(3)
n = 5*256
m = 6*256
d_max = m+n+1
q = 8380417
eta = 4
D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)}
D_e = D_s
mu_e, s_e = average_variance(D_e)
mu_s, s_s = average_variance(D_s)
sigma = sqrt(s_e)
print("-----Dilithium-II")
print("s_e = %.1f, s_s = %.1f" %(s_e,s_s))


##Core-SVP
#m_op,beta_op=Core_SVP(d_max,n,q,sigma)
#print("Core-SVP model \n Optimal d = %d, Optimal beta = %.1f" %(m_op+n+1,beta_op) )

#Lower_Bound_Estimation
print("\n Lower_Bound_Estimation")
d = d_max
M = sigma * sqrt(d)
V_d = q**(m/d_max)

Current_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
print(" Current beta = %d" % Current_Beta)
print("Lower_Bound_Estimation: d_svp = %d" % d_svp)
print("Lower_Bound_Estimation = %f" % T_LBE)


#Lower_Bound_Estimation_Optimal_m
print("\n Lower_Bound_Estimation_Optimal_m ")
d_op,beta_op,d_svp_op,T_LBE_op = Lower_Bound_Estimation_Optimal_d(m,n,sigma,q,m-297)
print("Lower_Bound_Estimation_Optimal_Dimension = %d" % d_op)
print("Lower_Bound_Estimation_Optimal_m: Current Beta = %d" % beta_op)
print("Lower_Bound_Estimation_Optimal_m: d_svp = %d" % d_svp_op)
print("Lower_Bound_Estimation_Optimal_m: Hardness = %f Bits" % T_LBE_op)






