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

    ldc_param = "MATZOV22": list decoding method in [MATZOV22]
                "AGPS20"(default): list decoding method in [AGPS20]
    
    cal_ee = "avg_sigma": martin's primal usvp + two step
             "chi" #chi-square + probabilistic + two-step
'''
#######################################
#Fixed parameters

method = 2
worst_case = False
#estimator in [DDGR20]
ldc_param = "AGPS20"
cal_ee = "chi" #chi-square + probabilistic + two-step
goal_min_cost = "gate_min" # "gate_min": find the minimal gates cost in two-step
                           # "gate_RAM_min": find the minimal (gates+RAM) cost in two-step 

#Martin's primal usvp + our two-step mode
ldc_param = "MATZOV22" #list decoding complexity proposed in [MATZOV22]
#cal_ee = "avg_sigma" #primal-martin-usvp + two-step
#------------------------------------




# WZD parameters
print("============= WZD")
q = 2**38
n = 1792
m = n
D_s = {x : 1./2 for x in range(0, 2)}
D_e =  {x : 1./2 for x in range(0, 2)}
dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s,m_range = -10)

#leaky-lwe-estimator
svp_estimate_attack( silent=False, method = 1,  dvol = dvol, dim_ = dim_,gen_GSA_gso = True,print_l = False ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = False, goal_min_cost = goal_min_cost)

#refined-two-step-lwe-estimator
svp_estimate_attack( silent=False, method = 2,  dvol = dvol, dim_ = dim_,gen_GSA_gso = True,print_l = False ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = False, goal_min_cost = goal_min_cost)



