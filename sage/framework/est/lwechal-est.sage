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



def low_dim_lwe_challenge_est(method,  cost_model, worst_case, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False):
    lwechals = [(40,0.025), (45,0.020), (50,0.015), (55,0.010), (60,0.010), (70,0.005), (75,0.005)]
    for (n,alpha) in lwechals:
        dim,dvol = gen_lwechal_instance(n, alpha)

        estimate_attack( silent=False, method = method, progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)



  


def unsolved_lwe_challenge_est(method,  cost_model, worst_case, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False):
    lwechals = [(40,0.045), (45,0.035), (50,0.030), (55,0.025), (60,0.020), (65,0.015), (75,0.010), (95,0.005)]
    for (n,alpha) in lwechals:
        dim,dvol = gen_lwechal_instance(n, alpha)

        estimate_attack( silent=False, method = method, progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

