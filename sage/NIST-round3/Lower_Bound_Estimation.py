
import math
from math import ceil, pi ,log2, e
from math import log ,sqrt


def delta(beta):
    return float( ( (pi*beta)**(1/beta)*(beta/(2*pi*e)) )**(1/(2*(beta-1))) )

def Core_SVP(d_max,n,q,sigma):
    beta_op=d_max
    m_op=d_max
    for m in range(n,n//2,-1):
        d = m+n+1
        for beta in range(50,d,1):
            if sigma*sqrt(beta) < delta(beta)**(2*beta-d-1)*q**(m/d):
                if beta_op >= beta:
                    beta_op=beta
                    m_op=m
                beta = d-1
    return m_op,beta_op

def minimum_sieving_dimension(M,d,delta0,V_d):
    for d_svp in range(2,d,1):
       #if M*sqrt(d_svp/d) <  V_d *( sqrt(d_svp/(2*pi*e))**(1/d) ) *( delta0**( (d_svp-d)/(d-1) )  ):
        if M*sqrt(d_svp/d) <  V_d *sqrt(d_svp/(2*pi*e)) *( delta0**( d_svp-d )  ):
            break
    return d_svp


def rhf(delta0,beta_1,d):
    return (delta0**((d-beta_1)/(d-1)) ) * (sqrt( beta_1/(2*pi*e) ) )**(1/d)


def Lower_Bound_Estimation(M,dimension,V_d):
    for beta in range(50,dimension,1):
        con = True
        d_svp=minimum_sieving_dimension(M,dimension,delta(beta),V_d)
        for beta_1 in range(beta+1,dimension,1):
            delta_1 = rhf(delta(beta),beta_1,dimension)
            d_svp_1 = minimum_sieving_dimension(M,dimension,delta_1,V_d)
            if 0.292*d_svp_1 > log2(dimension-beta_1+1)+0.292*beta_1:
                if 0.292*d_svp  > 0.292*d_svp_1 + log2(1+(dimension-beta_1+1) * 2**(0.292*(beta_1-d_svp_1))):
                    con = False
            elif (dimension-beta_1+1) * 2**(0.292*(d_svp_1-beta_1)) > 0:
                if 0.292*d_svp>log2(dimension-beta_1+1)+0.292*beta_1+log2(1+1/((dimension-beta_1+1) * 2**(0.292*(d_svp_1-beta_1)) )):
                    con = False
            elif 0.292*d_svp>log2(dimension-beta_1+1)+0.292*beta_1:
                    con = False
                #print("False Current beta=%d, beta_1=%d" %(beta,beta_1))
                #print("Left = %f, Right = %f"%(0.292*d_svp,log2(dimension-beta_1+1)+0.292*beta_1 + 0.292*d_svp_1))
                
                
        if con:
            #print("True! Current beta=%d, d_svp=%d" %(beta,d_svp))
            #print("True! beta_1=%d, d_svp_1=%d" %(beta_1,d_svp_1))
            #print("Left = %f, Right = %f"%(0.292*d_svp,log2(dimension-beta_1+1)+0.292*beta_1 + 0.292*d_svp_1))
            #print("Right_1 = %f, Right_2 = %f"%(log2(dimension-beta_1+1)+0.292*beta_1, 0.292*d_svp_1))
            return beta,d_svp,0.292*d_svp
        
def Lower_Bound_Estimation_Optimal_d(m_max,n,sigma,q,m_adps16):
    lst_d_svp = []
    lst_C_Beta = []
    lst_d = []   
    for m in range(m_max,m_adps16,-5):
        d_max = m+n+1
        d = d_max
        M = sigma * sqrt(d)
        V_d = q**(m/d)
        C_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
        lst_d_svp.append(d_svp)
        lst_C_Beta.append(C_Beta)
        lst_d.append(d_max)
        
    Index=lst_d_svp.index(min(lst_d_svp))
    d_svp_op = lst_d_svp[Index]
    beta_op = lst_C_Beta[Index]
    d_op = lst_d[Index]

    
    return d_op,beta_op,d_svp_op,0.292*d_svp_op


#Kyber-I
N= 256
k = 2
n = k*N
m = n
q = 3329
eta_1 = 3
eta_2 = 2
d_max = m+n+1
sigma = sqrt(eta_1/2)
print("-----Kyber-I")

sigma = sqrt(eta_1/2)
#Core-SVP    
m_op,beta_op=Core_SVP(d_max,n,q,sigma)
print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )


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
N= 256
k = 3
n = k*N
m = n
q = 3329
eta_1 = 2
eta_2 = 2
d_max = m+n+1
sigma = sqrt(eta_1/2)
print("-----Kyber-II")
#Core-SVP    
m_op,beta_op=Core_SVP(d_max,n,q,sigma)
print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )

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
N= 256
k = 4
n = k*N
m = n
q = 3329
eta_1 = 2
eta_2 = 2
d_max = m+n+1
sigma = sqrt(eta_1/2)
print("-----Kyber-III")
#Core-SVP    
m_op,beta_op=Core_SVP(d_max,n,q,sigma)
print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )

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
N= 256
l = 7
n = l*N
k = 8
m = k*N
q = 8380417
eta_1 = 2
eta_2 = 2
d_max = m+n+1
sigma = eta_1/sqrt(3)
print("-----Dilithium-III")


#Core-SVP
m_op,beta_op=Core_SVP(d_max,n,q,sigma)
print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )


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
N= 256
l = 4
n = l*N
k = 4
m = k*N
q = 8380417
eta_1 = 2
eta_2 = 2
d_max = m+n+1
sigma = eta_1/sqrt(3)
print("-----Dilithium-I")


#Core-SVP
m_op,beta_op=Core_SVP(d_max,n,q,sigma)
print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )


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
N= 256
l = 5
n = l*N
k = 6
m = k*N
q = 8380417
eta_1 = 4
eta_2 = 4
d_max = m+n+1
sigma = eta_1/sqrt(3)
print("-----Dilithium-II")


#Core-SVP
m_op,beta_op=Core_SVP(d_max,n,q,sigma)
print("Core-SVP model \n Optimal d = %d, Optimal beta = %d" %(m_op+n+1,beta_op) )

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






