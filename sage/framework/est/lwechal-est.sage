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


lwechals = [(40,0.005),(45,0.005),(50,0.005), (55,0.005), (40,0.015), (45,0.010)]

#[(40,0.005),(45,0.005),(50,0.005), (60,0.005)]#(40,0.025), (80,0.005)] #[(40,0.025), (80,0.005)] #[(40,0.025),(45,0.020),(50,0.015),(80,0.005)]

def lwechal_est(method,  cost_model, worst_case,J=1, gap=1, J_gap=1, parallel_ = False,  print_l = False, lwechals = lwechals, default_g6k=False, ldc_param = "AGPS20", cal_ee = "chi",goal_min_cost = "gate_min", cumG = False):
    for (n,alpha) in lwechals:
        (l,dim,dvol,b,sigma) = gen_lwechal_instance(n, alpha,default_g6k)
        estimate_attack( silent=False, method = method,  parallel_ = parallel_, l = l, dim_ = dim, dvol =dvol, J=J,gap =gap,J_gap = J_gap, cost_model=cost_model,gen_GSA_gso = False, print_l = print_l, worst_case = worst_case,b = b, sigma = sigma, ldc_param = ldc_param, cal_ee = cal_ee,goal_min_cost = goal_min_cost, cumG = cumG)

