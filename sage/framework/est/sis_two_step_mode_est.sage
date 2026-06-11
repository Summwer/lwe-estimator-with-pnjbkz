#two-step-mode-est for SIS 


load("../framework/utils.sage")


#find the minimal dimension for last pump using the upper bound of sis
def stev_d_svp_prediction(l,upper_bound, cost_model,ldc_param):
    """
    Dimension of sieve/progressive sieve chosen to find target vector.
    :l: log(||b_i^*||), i = 0,...,d-1
    :upper_bound: the square 2-norm upper bound of the SIS problem
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    :progressieve_sieve: True: progressieve sieve
                         False: normal sieve
    :ldc_param: list decoding param

    return value:
        dsvp: minimal dimension value to sieve
        G_sieve: the cumulated time cost for sieve
        B_dsvp: the maximal memory cost for sieve

    """

    d = len(l)
    l_ = [2*_*log(2.) for _ in l]

    #predict dimension of last sieve: progressive sieve
    dsvp = 50
    flag1 = False
    flag2 = False
    while(dsvp <= d):
        #2**(2 * l1[:dsvp])==2**(2 * l1dsvp)==gh
        gh = gaussian_heuristic(l_[:dsvp])
        
        #dsvp_prime = floor(dsvp - dim4free_wrapper(theo_dim4free_fun2,dsvp))
        if(cost_model == 2):
            dsvp_prime = floor(dsvp - dim4free_wrapper(default_dim4free_fun,dsvp))
        else:
            dsvp_prime = floor(dsvp - theo_dim4free_fun2(dsvp))
        
        if(gh <= upper_bound): 
            Gpump, Bsieve = pump_cost(dsvp_prime,cost_model = cost_model,ldc_param=ldc_param)

            return dsvp_prime, dsvp, Gpump, Bsieve   


def sis_two_step_mode_estimation(dvol, dim_, upper_bound, verbose=False, cost_model=1,ldc_param = "AGPS20", worst_case = False, goal_min_cost = "gate_min"):
    
    """
    SIS estimation: Simplified progressive BKZs + Pump to solve SIS problem.
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    :progressieve_sieve: True: progressieve sieve
                         False: normal sieve
    :goal_min_cost: "gate_min": find the minimal gates cost in two-step
                    "gate_RAM_min": find the minimal (gates+RAM) cost in two-step
    """
    Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = float("inf"), float("inf"), d, d, float("inf")
    betamin = []

    # Keep increasing beta to be sure to catch the second intersection
   
    betastart = 50
    out_loop = False
    for beta in range(betastart, min(d,1000)):
        #l = simulate_pnjBKZ(l, beta, 1, 1)
        delta = compute_delta(beta)
        #print(beta, delta, dvol, d)
        l = [log(bkzgsa_gso_len(dvol, i, d, delta=delta)) / log(2.) for i in range(d)]
        G1, B1 = bkz_cost(d,beta,ldc_param = ldc_param)
        G1 += log2(C)


        if(2**(2 * l[i]) <= upper_bound):
            G = G1
            B = B1
            dsvp = 0
            dsvp_prime = 0
            out_loop = True
        else:
            #d_svp prediction
            (dsvp, dsvp_prime, Gpump, Bpump) = sis_d_svp_prediction(l,cost_model,ldc_param)

            G = log2(2**G1 + 2**Gpump)
            B = max(B1, Bpump)

        if(goal_min_cost == "gate_min" and G < Gmin):
            Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = G, B, dsvp, dsvp_prime, G1
            betamin = list(range(betastart,beta+1))
        if(goal_min_cost == "gate_RAM_min" and G + B < Gmin + Bmin):
            Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = G, B, dsvp, dsvp_prime, G1
            betamin = list(range(betastart,beta+1))

       
        if verbose:
            print("β= %d, pump-{%d,%d,%d},  G=%3.2f gate,  B=%3.2f bit"%(beta, d-dsvp, dsvp, dsvp-dsvp_prime, Gmin, Bmin), end = "\r" if beta%50 != 0 else "\n")

        if(out_loop):
            break
    print()
    return betamin, G1min, dsvpmin, dsvp_prime_min, Gmin, Bmin


