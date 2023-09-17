import time
load("../framework/est/bkz_only_est.sage")
load("../framework/est/default_g6k_est.sage")
load("../framework/est/progressive_bkz_est.sage")
load("../framework/est/two_step_mode_est.sage")

def estimate_attack(silent=False, method=1, worst_case = False, parallel_ = False, l = None, dvol = None, dim_ = None, J=1, gap=1, J_gap = 1, cost_model = 1, gen_GSA_gso = True,print_l = False, b = None, sigma = 0., ldc_param = "AGPS20", cal_ee = "chi", goal_min_cost = "gate_min", cumG = False):
    """ Assesses the complexity of the lattice attack on the instance.

    method = 1: progressive bkz estimation.
             2: two step mode
             3: bkz-only with fixed blocksize
             4: default g6k 
    ldc_param: classical list decoding complexity, choose  "AGPS20"(https://eprint.iacr.org/2019/1161.pdf) or "MATZOV22"(https://marketing.idquantique.com/acton/attachment/11868/f-0587a79f-5592-47fe-9bdf-a3f3e7f7d802/1/-/-/-/-/Report%20on%20the%20Security%20of%20LWE.pdf)

    :goal_min_cost: "gate_min": find the minimal gates cost in two-step
                    "gate_RAM_min": find the minimal (gates+RAM) cost in two-step

    :cumG : False: [PV21] estimate
            True: cumulated gate estimate adaptive to [PV21]
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

    
    if(worst_case):
        print("Compute cost in worst case...")
    else:
        print("Compute cost in average case...")
    
    if(method == 1):
        if(not cumG):
            print(" Attack Estimation via progressive BKZ + chi_square mode (estimator in [PV21]) ", end = "")
        else:
            print(" Attack Estimation via progressive BKZ + chi_square mode (cumulated gate estimator adaptive to [PV21]) ", end = "")
        if(cost_model == 1):
            print(", list_decoding_classical = ", ldc_param)
        print()
        betamin, G, B = pro_bkz_est(l, verbose=not silent, cost_model = cost_model, worst_case = worst_case, cumG = cumG)
        
        print("Cost: %.3f s" %(time.time()-T0))

    if(method == 2):
        print(" Attack Estimation via two-step mode (simply progressive bkz+pump)", end= "")
        if(cost_model == 1):
            print(", list_decoding_classical = ", ldc_param, end= ", ")
        if(cal_ee == "chi"):
            print("cumulated probalistic estimator form in [DDGR20],", end="")
        if(cal_ee == "avg_sigma"):
            print("estimator form in [ADPS16],", end="")
        print("goal_min_cost = ", goal_min_cost)
        betamin, G1, dsvp, dsvp_prime, G, B = two_step_mode_estimation(l, dvol, dim_, verbose=not silent, cost_model=cost_model,ldc_param = ldc_param, cal_ee = cal_ee,  worst_case = worst_case, goal_min_cost = goal_min_cost)
        print("Cost: %.3f s" %(time.time()-T0))

 
    if(method == 3):
        print(" Attack Estimation via bkz-only mode with fixed blocksize", end="")
        if(cost_model == 1):
            print(", list_decoding_classical = ", ldc_param)
        print()

        betamin , G, B = bkz_only(dim_, dvol, l,verbose=not silent,cost_model=cost_model, worst_case = worst_case,ldc_param = ldc_param)
        print("Cost: %.3f s" %(time.time()-T0))
        
    if(method==4):
        print(" Attack Estimation via defaul g6k mode  ", end="")
        if(cost_model == 1):
            print(", list_decoding_classical = ", ldc_param)
        print()
        betamin, G, B = default_g6k_est(dim_, dvol, b, l, verbose=not silent, cost_model=cost_model, worst_case = worst_case, sigma = sigma)
        print("Cost: %.3f s" %(time.time()-T0))
        

    if not silent:
        if cost_model == 1:
            time_unit = "G = %3.2f log2(gate) \t B = %3.2f log2(bit)" %(G,B)
        elif cost_model == 2:
            time_unit = "G = %3.2f sec \t B = %3.2f GB" %((2**G), 2**(B-33))
        if method == 1:
            logging("dim=%3d \t ln(dvol) = %4.7f \t avg_beta=%3.2f \t %s " % (dim_, dvol, betamin, time_unit), style="VALUE")
        if method ==3:
            logging("dim=%3d \t ln(dvol) = %4.7f \t beta=%d \t %s " % (dim_, dvol, round(betamin), time_unit), style="VALUE")
        elif method == 2:
            if(cost_model == 1):
                logging("dim=%3d \t ln(dvol) = %4.7f \t β_strategy=[%s,...,%s] \t dsvp = %d \t dsvp_prime = %d \t %s " %(dim_, dvol, str(betamin[0]), str(betamin[-1]), dsvp, dsvp_prime, time_unit), style="VALUE")
            if(cost_model == 2):
                logging("dim=%3d \t ln(dvol) = %4.7f \t β_strategy=[%s,...,%s] \t dsvp = %d \t dsvp_prime = %d \t %s " %(dim_, dvol, str(betamin[0]), str(betamin[-1]), dsvp, dsvp_prime, time_unit), style="VALUE")
        elif method == 4:
            logging("dim=%3d \t ln(dvol) = %4.7f \t β_strategy=[%s,...,%s] \t %s " % (dim_, dvol, str(betamin[0]), str(betamin[-1]), time_unit), style="VALUE")
        logging("\n")


def div_sigma(l, sigma):
    """
    input: l -- log2(||b_i^*||)
           sigma -- Standard Deviation in gaussian distribution
    """

    return [_-log2(sigma) for _ in l]
