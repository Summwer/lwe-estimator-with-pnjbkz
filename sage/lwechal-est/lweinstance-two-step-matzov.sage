load("../framework/est/lwechal-est.sage")
load("../framework/est/NIST-est.sage")
######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3: bkz-only with fixed blocksize
             4: default g6k
             
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
method = 2 # 2 - two-step trivial; 4 - default g6k 
cost_model = 1 #theoretical cost model
ldc_param = "MATZOV22"
worst_case = False
cal_ee = "chi"
goal_min_cost = "gate_min"

#------------------------------------


t = 1600
method = 2        

for sigma in [1]:
    for log2q in range(41, 40, -1):
        q = 2** log2q
        for n in range(600,601, 50):
            for success_count in [10]:
            #for success_count in range(10, 1, -1):
                m = n
                
                success_rate = 0.1 * success_count - 1e-3
               
                D_s = build_Gaussian_law(sigma, t)
                D_e = D_s
                dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
                svp_estimate_attack( silent=False, method = method, dvol = dvol, dim_ = dim_,  cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cost_model = cost_model, success_rate=success_rate, ldc_param = ldc_param)
