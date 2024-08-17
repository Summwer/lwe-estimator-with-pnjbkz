import math
from math import ceil, pi ,log2, e
from math import log ,sqrt
import matplotlib.pylab as plt
import numpy as np


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


def rhf(delta0,beta_1,d):
    return (delta0**((d-beta_1)/(d-1)) ) * (sqrt( beta_1/(2*pi*e) ) )**(1/d)


def logG1(tours,beta,d):
    return log2(5.46*tours*(d-beta+1))+0.292*beta

def d4f_delta(beta,V_d,sigma):
    return int(log(V_d*sqrt(2/(3*pi*e))/sigma) / log(delta(beta)))

def logG2(beta,d,V_d,sigma):
    return 0.292*(d-d4f_delta(beta,V_d,sigma))

def Log_G_beta(d,V_d,sigma,tours):
    log_G_op = 0.292*d
    beta_op = d
    log_G_beta = []
    for beta in range(50,d):
        if logG1(tours,beta,d)>logG2(beta,d,V_d,sigma):
            log_G_beta.append( logG1(tours,beta,d)+log(1+2**(logG2(beta,d,V_d,sigma)-logG1(tours,beta,d))) )
            #print("Case 1 %d" %beta)
            if log_G_op > log_G_beta[beta-50]:
                log_G_op = log_G_beta[beta-50]
                beta_op = beta      
        else:
            #print("Case 2 %d" %beta)
            log_G_beta.append( logG2(beta,d,V_d,sigma)+log(1+2**(logG1(tours,beta,d)-logG2(beta,d,V_d,sigma))) )
            if log_G_op > log_G_beta[beta-50]:
                log_G_op = log_G_beta[beta-50]
                beta_op = beta
    print('Optimal log T = %f, Optimal beta = %d'%(log_G_op,beta_op))
    return log_G_op,beta_op,log_G_beta

#Figure 
def Figure_Two_step(d,log_G_beta,beta_leaky,log_G_op,beta_op,NIST):
    x = [_ for _ in range(50,d)]
    plt.figure(figsize=(9, 6), dpi=1000)
#    plt.xlabel(r'Quality of Current Lattice Basis with $\beta$-reduced (%s parameter)'%NIST, size=20)
    plt.xlabel(r'$\beta$', size=20)
    plt.ylabel(r'$log(T=T_1+T_2)$', size=20)
    plt.plot(x, log_G_beta, '*k', label = "Two-step") #label=r'$log(T=T_1(\beta)+T_2(\beta))$ based on BKZ-$\beta$ reduced basis')
    logT_leak_LWE = log2(5.46*(d- beta_leaky+1)) + beta_leaky*0.292
    plt.plot([50,d],[logT_leak_LWE,logT_leak_LWE],c='k',linestyle='--',label = r"ADPS16,$\log_2{\rm T}$ = %.1f" %logT_leak_LWE)# label=r'Progressive BKZ-$\beta_0$ (ADPS16 $\beta_0$) logT=%f'%(logT_leak_LWE) )
    plt.plot([beta_op,beta_op],[log_G_op-2,log_G_beta[0]],c='r',linestyle='--', label = r"$\beta_{\rm op}$ = %d, $\log_2{\rm T}$ = %.1f" %(beta_op, log_G_op))# label='Optimal beta in Two-step Mode, logT = %f'%log_G_op)

    plt.legend(fontsize = 20)
    plt.grid()
    plt.ylim(None,max(log_G_beta)+ 5)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig("Log T with different reduction basis (%s).png"%NIST, bbox_inches='tight')
    plt.close()

#Figure Zoom in
def Figure_Two_step_Zoom_in(d,log_G_beta,beta_0,log_G_op,beta_op,xleft,xright,Zoom_in_yfactor,NIST):
    x = [_ for _ in range(beta_op-xleft,beta_op+xright)]
    plt.figure(figsize=(9, 6), dpi=1000)
#    plt.xlabel(r'Quality of Current Lattice Basis with $\beta$-reduced (%s parameter)'%NIST, size=20)
    plt.xlabel(r'$\beta$', size=20)
    plt.ylabel(r'$log(T=T_1+T_2)$', size=20)
    logT_leak_LWE = log2(5.46*(d- beta_0+1)) + beta_0*0.292
#    flag = True
    small_x = [_ for _ in x if log_G_beta[_-50] <= logT_leak_LWE]
    large_x = [_ for _ in x if log_G_beta[_-50] > logT_leak_LWE]
    small_T = [log_G_beta[_-50] for _ in x if log_G_beta[_-50] <= logT_leak_LWE]
    large_T = [log_G_beta[_-50] for _ in x if log_G_beta[_-50] > logT_leak_LWE]
#    print(large_x)
    for i in range(len(large_x)-1):
        if(large_x[i+1]- large_x[i] > 1):
            break
    plt.plot(small_x, small_T,color='b', marker='o', markerfacecolor='none', markersize='6', label = "Two-step")
    plt.plot(large_x[:i+1]+[large_x[i]+1], large_T[:i+1]+[log_G_beta[large_x[i]+1-50]],color='k', marker='o', markerfacecolor='none', markersize='6')
    plt.plot( [large_x[i+1]-1] + large_x[i+1:], [log_G_beta[large_x[i+1]-1-50]]+large_T[i+1:],color='k', marker='o', markerfacecolor='none', markersize='6')
#    for _ in x:
#        if logT_leak_LWE < log_G_beta[_-50]:
#            if flag:
#                plt.plot(_, log_G_beta[_-50],color='k', marker='o', markerfacecolor='none', markersize='6', label = "Two-step") #label=r'$log(T=T_1(\beta)+T_2(\beta))$ based on BKZ-$\beta$ reduced basis')
#                flag = False
#            else:
#                plt.plot(_, log_G_beta[_-50],color='k', marker='o', markerfacecolor='none', markersize='6')
#        else:
#            plt.plot(_, log_G_beta[_-50],color='b', marker='o', markerfacecolor='none', markersize='6')
    
    logT_leak_LWE = log2(5.46*(d- beta_0+1)) + beta_0*0.292
    plt.plot([beta_op-xleft,beta_op+xright],[logT_leak_LWE,logT_leak_LWE],c='k',linestyle='--', label = r"ADPS16,$\log_2{\rm T}$ = %.1f" %logT_leak_LWE)#, label=r'Progressive BKZ-$\beta_0$ (ADPS16 $\beta_0$) logT=%f'%(logT_leak_LWE))
    plt.plot([beta_op,beta_op],[log_G_op,log_G_op+Zoom_in_yfactor],c='r',linestyle='-.', label = r"$\beta_{\rm op}$ = %d, $\log_2{\rm T}$ = %.1f" %(beta_op, log_G_op))# label='Optimal beta in Two-step Mode, logT = %f'%log_G_op)
    plt.plot([beta_op],[log_G_op], color = 'red',marker="o", markersize='6', zorder = 2)#markerfacecolor='none'
    plt.legend(fontsize = 20)
    plt.ylim(None,max([log_G_beta[_] for _ in x]))
    plt.grid()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig("Log T with different reduction basis Zoom in (%s).png"%NIST, bbox_inches='tight')
    plt.close()
    
    
# def Figure_Two_step_Zoom_in(d,log_G_beta,beta_0,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST):
#     x = [_ for _ in range(beta_op-Zoom_in_xfactor,beta_op+Zoom_in_xfactor)]
#     plt.figure(figsize=(15, 10), dpi=100)
#     plt.xlabel('Quality of Current Lattice Basis with beta-reduced (%s parameter)'%NIST, size=20)
#     plt.ylabel(r'$log(T=T_1+T_2)$', size=20)
#     plt.plot(x, log_G_beta[beta_op-50-Zoom_in_xfactor:beta_op-50+Zoom_in_xfactor],color='k', marker='o', markerfacecolor='none', markersize='6', label='Log T with different reduction basis')
#     logT_leak_LWE = log2(5.46*(d- beta_0+1)) + beta_0*0.292
#     plt.plot([beta_op-Zoom_in_xfactor,beta_op+Zoom_in_xfactor],[logT_leak_LWE,logT_leak_LWE],c='k',linestyle='--', label='Leaky-LWE logT=%f'%(logT_leak_LWE))
#     plt.plot([beta_op,beta_op],[log_G_op,log_G_op+Zoom_in_yfactor],c='r',linestyle='-.', label='Optimal beta in Two-step Mode, logT = %f'%log_G_op)
#     plt.legend(fontsize = 20)
#     plt.grid()
#     plt.savefig("Log T with different reduction basis Zoom in (%s).png"%NIST)
#     plt.close()
    
    
    
# Heuristic need tours times progressive BKZ-beta to achive BKZ-beta reduced
tours = 1
 
#Kyber-I
N = 256
k = 2
n = k*N
m = n
q = 3329
eta_1 = 3
eta_2 = 2
d_max = m+n+1
sigma = sqrt(eta_1/2)
print("-----Kyber-I")
NIST = "Kyber-I"

sigma = sqrt(eta_1/2)

Kyber_beta_leaky = 413
Core_SVP_beta = 403
#
d = d_max
V_d = q**(m/d_max)

log_G_op,beta_op,log_G_beta = Log_G_beta(d,V_d,sigma,tours)


#Figure
#Figure_Two_step(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,NIST)
Figure_Two_step(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,NIST)


#Figure Zoom in
#Zoom_in_xfactor=40
Zoom_in_yfactor=12
#Figure_Two_step_Zoom_in(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST)
Figure_Two_step_Zoom_in(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,40,40,Zoom_in_yfactor,NIST)
            

    
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
NIST = "Kyber-II"

sigma = sqrt(eta_1/2)

Kyber_beta_leaky = 637
Core_SVP_beta = 625


#
d = d_max
V_d = q**(m/d_max)

log_G_op,beta_op,log_G_beta = Log_G_beta(d,V_d,sigma,tours)


#Figure
#Figure_Two_step(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,NIST)
Figure_Two_step(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,NIST)

#Figure Zoom in
#Zoom_in_xfactor=60
Zoom_in_yfactor=16
#Figure_Two_step_Zoom_in(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST)
Figure_Two_step_Zoom_in(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,50,50,Zoom_in_yfactor,NIST)



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
NIST = "Kyber-III"

sigma = sqrt(eta_1/2)

Kyber_beta_leaky = 894
Core_SVP_beta = 877

#
d = d_max
V_d = q**(m/d_max)

log_G_op,beta_op,log_G_beta = Log_G_beta(d,V_d,sigma,tours)


#Figure
#Figure_Two_step(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,NIST)
Figure_Two_step(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,NIST)

#Figure Zoom in
#Zoom_in_xfactor=70
Zoom_in_yfactor=19
#Figure_Two_step_Zoom_in(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST)
Figure_Two_step_Zoom_in(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,60,60,Zoom_in_yfactor,NIST)


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
NIST = "Dilithium-I"


beta_leaky = 433
Core_SVP_beta = 423

#
d = d_max
V_d = q**(m/d_max)

log_G_op,beta_op,log_G_beta = Log_G_beta(d,V_d,sigma,tours)


#Figure
#Figure_Two_step(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,NIST)
Figure_Two_step(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,NIST)

#Figure Zoom in
#Zoom_in_xfactor=25
Zoom_in_yfactor=40
#Figure_Two_step_Zoom_in(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST)
Figure_Two_step_Zoom_in(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,20,40,Zoom_in_yfactor,NIST)




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
NIST = "Dilithium-II"


beta_leaky = 638
Core_SVP_beta = 624
#
d = d_max
V_d = q**(m/d_max)

log_G_op,beta_op,log_G_beta = Log_G_beta(d,V_d,sigma,tours)


#Figure
#Figure_Two_step(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,NIST)
Figure_Two_step(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,NIST)

#Figure Zoom in
#Zoom_in_xfactor=25
Zoom_in_yfactor=20
#Figure_Two_step_Zoom_in(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST)
Figure_Two_step_Zoom_in(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,15,40,Zoom_in_yfactor,NIST)




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
NIST = "Dilithium-III"


beta_leaky = 883
Core_SVP_beta = 863
#
d = d_max
V_d = q**(m/d_max)

log_G_op,beta_op,log_G_beta = Log_G_beta(d,V_d,sigma,tours)


#Figure
#Figure_Two_step(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,NIST)
Figure_Two_step(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,NIST)


#Figure Zoom in
#Zoom_in_xfactor=20
Zoom_in_yfactor=20
#Figure_Two_step_Zoom_in(d,log_G_beta,Kyber_beta_leaky,log_G_op,beta_op,Zoom_in_xfactor,Zoom_in_yfactor,NIST)
Figure_Two_step_Zoom_in(d,log_G_beta,Core_SVP_beta,log_G_op,beta_op,25,55,Zoom_in_yfactor,NIST)



