import time
load("../framework/attack.sage")
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


def dilithium_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False, ldc_param = "AGPS20", cal_ee = "chi"):

    # Dilithium-I round-3 parameters
    print("============= Dilithium-I")
    dim_ = 2049
    dvol = 15614.2193172
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)

    
    print("============= Dilithium-II")

    # Dilithium-II round-3 parameters
    dim_ = 2654
    dvol = 19371.0238433
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)


    print("============= Dilithium-III")

    # Dilithium-III round-3 parameters
    dim_ = 3540
    dvol = 26623.1162463
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l,  ldc_param = ldc_param, cal_ee = cal_ee)



def kyber_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False, ldc_param = "AGPS20", cal_ee = "chi"):

    # Kyber-512 round-3 parameters
    print("============= Kyber-512")
    #eta1 = eta2 = 3
    #dim_ = 1025
    #dvol = 3944.9406103
    #eta1 = 3, eta2 = 2
    dim_ = 1004
    dvol = 3882.6780896
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee)


    print("============= Kyber-768")

    # Kyber-768 round-3 parameters
    dim_ = 1467
    dvol = 5661.0782118
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)


    print("============= Kyber-1024")

    # Kyber-1024 round-3 parameters
    dim_ = 1918
    dvol = 7242.6115232 
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)

   

def frodo_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False,ldc_param = "AGPS20", cal_ee = "chi"):

    # Frodo-640 round-3 parameters
    print("============= Frodo-640")
    dim_ = 1297
    dvol = 5479.4593497
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)


    print("============= Frodo-976")

    # Frodo-976 round-3 parameters
    dim_ = 1969
    dvol = 9347.2957371 
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)


    print("============= Frodo-1344")

    # Frodo-1344 round-3 parameters
    dim_ = 2634
    dvol = 13355.2889193 
    estimate_attack( silent=False, method = method,progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee)

   

