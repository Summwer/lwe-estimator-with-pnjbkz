import time

load("../framework/est/bkz_only_est.sage")
load("../framework/est/progressive_bkz_only_est.sage")
load("../framework/est/default_g6k_est.sage")
load("../framework/est/two_step_mode_est.sage")
#load("../framework/est/bssa_est.sage")
#load("../framework/est/enumbs_est.sage")



def estimate_attack(silent=False, method=1, progressive_sieve = True,parallel_ = False, l = None, dvol = None, dim_ = None, J=1, gap=1, J_gap = 1, cost_model = 1, gen_GSA_gso = True,print_l = False, b = None):
    """ Assesses the complexity of the lattice attack on the instance.
    Return value in Bikz

    method = 1: progressive bkz estimation.
             2: two step mode
             3, parallel_ = False: EnumBS
             3, parallel_=True: EnumBS in parallel
             4: bssa
             5: bkz-only with fixed blocksize
             6: default g6k 
    progressive_sieve = True:  progressive sieve
    """

    if l == None or gen_GSA_gso:
        print("Generate gs-lengths by GSA assumption.")
        delta = compute_delta(2)
        l = [log(bkzgsa_gso_len(dvol, i, dim_, delta=delta)) / log(2) for i in range(dim_)]
    else:
        print("Input real rr values.")
        dim_ = len(l)
    if(print_l):
        print(l)
    
    T0 = time.time()
    
    if(method == 1):
      
        print(" Attack Estimation via simulation + probabilistic model (cum_G) ")
        betamin, G, B = progressive_bkz_only(l, verbose=not silent, cost_model = cost_model)
        print("Cost in simulation + probabilistic model (cum_G): %.3f s" %(time.time()-T0))



    if(method == 2):
        print(" Attack Estimation via simulation + probabilistic model (simply progressive bkz+pump)")
        betamin, G1, dsvp, dsvpmax, G, B, delta = two_step_mode_estimation( l, verbose=not silent, progressive_sieve = progressive_sieve, cost_model=cost_model)
        print("Cost in simple two-mode attack: %.3f s" %(time.time()-T0))


    if(method == 3):
        
        enumBS = EnumBS()
        if(not parallel_):
            print(" Attack Estimation via simulation + probabilistic model (EnumBS)")
            betamin, G1, dsvp, dsvpmax, G, B, delta = enumBS.EnumBS_estimation(dim_, dvol, l, verbose=not silent, progressive_sieve = progressive_sieve, J=J,gap=gap,J_gap=J_gap,cost_model=cost_model)
            print("Cost in EnumBS attack: %.3f s" %(time.time()-T0))
        else:
            print(" Attack Estimation via simulation + probabilistic model (EnumBS in parallel)")
            betamin, G1, dsvp, dsvpmax, G, B, delta = enumBS.EnumBS_estimation_in_parallel(dim_, dvol, l, verbose=not silent, progressive_sieve = progressive_sieve, J=J,gap=gap,J_gap=J_gap,cost_model=cost_model)
            print("Cost in EnumBS attack: %.3f s" %(time.time()-T0))

    if(method == 4):
        print(" Attack Estimation via simulation + probabilistic model (BSSA)")
        bssa = BSSA()
        betamin, dsvp, dsvpmax, G, G1, B, delta = bssa.BSSA_est(dim_, dvol, l, verbose=not silent, progressive_sieve = progressive_sieve,J=J,gap=gap,J_gap=J_gap,cost_model=cost_model)
        print("Cost in BSSA attack: %.3f s" %(time.time()-T0))
    if(method == 5):
        print(" Attack Estimation via simulation + probabilistic model (fixed blocksize bkz-only)")

        betamin , G, B = bkz_only(dim_, dvol, l,verbose=not silent,cost_model=cost_model)
        print("Cost in fixed blocksize bkz-only: %.3f s" %(time.time()-T0))
        
    if(method==6):
        print(" Attack Estimation via simulation + probabilistic model (default g6k)")
        betamin, G, B = default_g6k_est(dim_, dvol, b, l, verbose=not silent, progressive_sieve = progressive_sieve,cost_model=cost_model)
        print("Cost in default g6k: %.3f s" %(time.time()-T0))

        

    if not silent:
        if cost_model == 1:
            time_unit = "log2(gate)"
        elif cost_model == 2:
            time_unit = "log2(sec)"
        if method == 1 or method ==5:
            logging("dim=%3d \t ln(dvol) = %4.7f \t avg_beta=%3.2f \t G=%3.2f %s \t B=%3.2f log2(bit) " % (dim_, dvol, betamin, G, time_unit, B), style="VALUE")
        elif method >= 2 and method <=4:
            logging("dim=%3d \t ln(dvol) = %4.7f \t β_strategy=[%s,...,%s] \t G1=%3.2f %s \t dim_svp=%3.2f \t dim_svp_max=%3d \t G_sieve=%3.2f %s \t G=%3.2f %s \t B=%3.2f log2(bit) " %
                             (dim_, dvol, str(betamin[0]), str(betamin[-1]),G1,time_unit, dsvp, dsvpmax, log2(2**G-2**G1),time_unit, G,time_unit, B), style="VALUE")
        elif method == 6:
            logging("dim=%3d \t ln(dvol) = %4.7f \t β_strategy=[%s,...,%s] \t G=%3.2f %s \t B=%3.2f log2(bit) " % (dim_, dvol, str(betamin[0]), str(betamin[-1]), G, time_unit, B), style="VALUE")
        logging("\n")


def div_sigma(l, sigma):
    """
    input: l -- log2(||b_i^*||)
           sigma -- Standard Deviation in gaussian distribution
    """

    return [_-log2(sigma) for _ in l]
