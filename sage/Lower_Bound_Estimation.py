import math
from math import ceil, pi ,log2, e
from math import log ,sqrt
from multiprocessing import Pool
from math import factorial as fac
    
def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2 * k, x + k) / 2.**(2 * k)
    
#binomial distribution
def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom

def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k + 1):
        D[i] = centered_binomial_pdf(k, i)
    return D
    

def average_variance(D):
    mu = 0.
    s = 0.
    
    for (v, p) in D.items():
        mu += v * p
        s += v * v * p

    s -= mu * mu
    return mu, s


def compute_delta(beta):
    return float( ( (pi*beta)**(1/beta)*(beta/(2*pi*e)) )**(1/(2*(beta-1))) )

#def GH_sv_factor_squared(k):
#    return ((pi * k)**(1. / k) * k / (2. * pi * e))

#def compute_delta(k):
#    """Computes delta from the block size k. Interpolation from the following
#    data table:
#    Source : https://bitbucket.org/malb/lwe-estimator/
#    src/9302d4204b4f4f8ceec521231c4ca62027596337/estima
#    tor.py?at=master&fileviewer=file-view-default
#    :k: integer
#    estimator.py table:
#    """
#
#    small = {0: 1e20, 1: 1e20, 2: 1.021900, 3: 1.020807, 4: 1.019713, 5: 1.018620,
#             6: 1.018128, 7: 1.017636, 8: 1.017144, 9: 1.016652, 10: 1.016160,
#             11: 1.015898, 12: 1.015636, 13: 1.015374, 14: 1.015112, 15: 1.014850,
#             16: 1.014720, 17: 1.014590, 18: 1.014460, 19: 1.014330, 20: 1.014200,
#             21: 1.014044, 22: 1.013888, 23: 1.013732, 24: 1.013576, 25: 1.013420,
#             26: 1.013383, 27: 1.013347, 28: 1.013310, 29: 1.013253, 30: 1.013197,
#             31: 1.013140, 32: 1.013084, 33: 1.013027, 34: 1.012970, 35: 1.012914,
#             36: 1.012857, 37: 1.012801, 38: 1.012744, 39: 1.012687, 40: 1.012631,
#             41: 1.012574, 42: 1.012518, 43: 1.012461, 44: 1.012404, 45: 1.012348,
#             46: 1.012291, 47: 1.012235, 48: 1.012178, 49: 1.012121, 50: 1.012065}
#
#    if k != round(k):
#        x = k - floor(k)
#        d1 = compute_delta(floor(k))
#        d2 = compute_delta(floor(k) + 1)
#        return x * d2 + (1 - x) * d1
#
#    k = int(k)
#    if k < 50:
#        return small[k]
#    else:
#        delta = GH_sv_factor_squared(k)**(1. / (2. * k - 2.))
#        return delta

def Core_SVP(d_max,n,q,sigma):
    beta_op=d_max
    m_op=d_max
    for m in range(n,n//2,-1):
        d = m+n+1
        for beta in range(50,d,1):
            if sigma*sqrt(beta) <= (compute_delta(beta)**(2*beta-d-1))*(q**(m/d)):
                if beta_op > beta:
                    beta_op=beta
                    m_op=m
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
        d_svp=minimum_sieving_dimension(M,dimension,compute_delta(beta),V_d)
        for beta_1 in range(beta+1,dimension,1):
            delta_1 = rhf(compute_delta(beta),beta_1,dimension)
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
            if not con:
                break
                
                
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
#    for m in range(m_adps16+300,m_adps16,-30):
    for m in range(m_max,m_adps16,-5):
        d_max = m+n+1
        d = d_max
        M = sigma * sqrt(d)
        V_d = q**(m/d)
        C_Beta,d_svp,T_LBE = Lower_Bound_Estimation(M,d_max,V_d)
        lst_d_svp.append(d_svp)
        lst_C_Beta.append(C_Beta)
        lst_d.append(d_max)
#
#    param_set = []
#    lst_d = []
#    for m in range(m_adps16+300,m_adps16,-30):
#        d_max = m+n+1
#        d = d_max
#        M = sigma * sqrt(d)
#        V_d = q**(m/d)
#        param_set.append((M,d_max,V_d))
#        lst_d.append(d_max)
#
#    with Pool(len(lst_d)) as p:
#        result = p.map(Lower_Bound_Estimation,param_set)
    
#    lst_C_Beta = [_[0] for _ in range(len(result))]
#    lst_d_svp = [_[1] for _ in range(len(result))]
        
    Index=lst_d_svp.index(min(lst_d_svp))
    d_svp_op = lst_d_svp[Index]
    beta_op = lst_C_Beta[Index]
    d_op = lst_d[Index]

    
    return d_op,beta_op,d_svp_op,0.292*d_svp_op


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






