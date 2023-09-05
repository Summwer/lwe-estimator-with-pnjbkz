import time
load("../framework/est/NIST-est.sage")

######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3: bkz-only with fixed blocksize
             4: default g6k
             "fixed-blocksize bkz-only": call fixed-blocksize bkz-only model to estimate security.
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################
#Fixed parameters

method = 2

#estimator in [DDGR20]
#ldc_param = "AGPS20"
#cal_ee = "chi" #chi-square + probabilistic + two-step


#Martin's primal usvp + our two-step mode
ldc_param = "MATZOV22"
cal_ee = "avg_sigma" #primal-martin-usvp + two-step
#------------------------------------

kyber_est(method,  ldc_param =  ldc_param, cal_ee = cal_ee)
dilithium_est(method, ldc_param =  ldc_param, cal_ee = cal_ee)








