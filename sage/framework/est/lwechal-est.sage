load("../framework/attack.sage")
load("../framework/lwe_estimation.sage")
load("../framework/load_lwechal.sage")

######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3, parallel_ = False: EnumBS
             3, parallel_=True: EnumBS in parallel
             4: bssa
             5: bkz-only with fixed blocksize
             6: default g6k
             
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################


def lwechal_est(lwechals, method,  cost_model, worst_case, parallel_ = False,  print_l = False, ldc_param = "AGPS20", cal_ee = "chi",goal_min_cost = "gate_min", cumG = False, verbose = False, silent = False):
    for (n,alpha) in lwechals:
        (l,dim,dvol,b,sigma) = gen_lwechal_instance(n, alpha)

        svp_estimate_attack( silent=False, method = method, dvol = dvol, dim_ = dim, ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cost_model = cost_model)

