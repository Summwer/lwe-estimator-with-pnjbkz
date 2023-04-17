load("../framework/attack.sage")


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

    n = 40
    alpha = 0.025
    dvol=331.9735338747315
    dim=172
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method, progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 45
    alpha = 0.020
    dvol=373.4658972674150
    dim=185
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 50
    alpha = 0.015
    dvol=415.6552752456852
    dim=194
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 55
    alpha = 0.010
    dvol=495.0168620849964
    dim=205
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 60
    alpha = 0.010
    dvol=522.7192487542407
    dim=222
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 70
    alpha = 0.005
    dvol=641.7748006815514
    dim= 235
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)


    n = 75
    alpha = 0.005
    dvol=678.7288625607595
    dim= 251
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)





def unsolved_lwe_challenge_est(method,  cost_model, worst_case, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False):

    n = 40
    alpha = 0.045
    dvol=302.1993617
    dim=195
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method, progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 45
    alpha = 0.035
    dvol=357.0995642
    dim=211
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 50
    alpha = 0.030
    dvol=400.4076907
    dim=228
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 55
    alpha = 0.025
    dvol=439.9769224
    dim=241
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 60
    alpha = 0.020
    dvol=494.0253108
    dim=254
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 65
    alpha = 0.015
    dvol=557.6405653
    dim= 262
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)


    n = 75
    alpha = 0.010
    dvol=637.6057085
    dim= 281
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

    n = 95
    alpha = 0.005
    dvol=836.9696072
    dim= 323
    print("TU LWE Challenge, n = %3d, alpha = %3.3f" %(n,alpha))
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))

    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)


