load("../framework/utils.sage")
load("../framework/d_svp_prediction.sage")



def two_step_mode_estimation(l, dvol, dim_, verbose=False, cost_model=1,ldc_param = "APGS20", cal_ee = "chi", worst_case = False):
    if(cal_ee == "chi"):
        return proba_two_step_mode_estimation(l, verbose=verbose, cost_model=cost_model, ldc_param = ldc_param, worst_case = worst_case)
    if(cal_ee == "avg_sigma"):
        return stev_two_step_mode_estimation(dvol, dim_, verbose=verbose, cost_model=cost_model, ldc_param = ldc_param)

def stev_two_step_mode_estimation(dvol, d, verbose=False, cost_model=1,ldc_param = "APGS20"):
    
    """
    LWE estimation: Simplified progressive BKZs + Pump through average sigma value, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    :progressieve_sieve: True: progressieve sieve
                         False: normal sieve
    """
    Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = float("inf"), float("inf"), d, d, float("inf")
    betamin = []

    # Keep increasing beta to be sure to catch the second intersection
   
    betastart = 50
    for beta in range(betastart, min(d,1000)):
        #l = simulate_pnjBKZ(l, beta, 1, 1)
        delta = compute_delta(beta)
        #print(beta, delta, dvol, d)
        l = [log(bkzgsa_gso_len(dvol, i, d, delta=delta)) / log(2.) for i in range(d)]
        G1, B1 = pro_theo_bkz_cost(d,beta,ldc_param = ldc_param)
       
        #d_svp prediction
        (dsvp, dsvp_prime, Gpump, Bpump) = stev_d_svp_prediction(l,cost_model,ldc_param)

        G = log2(2**G1 + 2**Gpump)
        B = max(B1, Bpump)

        if(Gmin > G):
            Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = G, B, dsvp, dsvp_prime, G1
            betamin = list(range(betastart,beta+1))
       
        if verbose:
            print("β= %d, pump-{%d,%d,%d},  G=%3.2f gate,  B=%3.2f bit"%(beta, d-dsvp, dsvp, dsvp-dsvp_prime, Gmin, Bmin), end = "\r" if beta%50 != 0 else "\n")
    print()
    return betamin, G1min, dsvpmin, dsvp_prime_min, Gmin, Bmin


def proba_two_step_mode_estimation(l, verbose=False, cost_model=1,ldc_param = "APGS20", worst_case = False):
    """
    LWE estimation: Simplified progressive BKZs + Pump throught chi_square distribution, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    :progressieve_sieve: True: progressieve sieve
                         False: normal sieve
    """
    d = len(l)
    DDGR20 = False
    bbeta = None
    pprev_margin = None
    Gmin, Bmin, avgbetamin, dsvpmin = float("inf"), float("inf"), d, d

    # Keep increasing beta to be sure to catch the second intersection
   
    remaining_proba = 1.
    average_beta = 0.
    cumulated_proba = 0.
    G1cum,B1cum = 0.,0.
    GBKZ = 0.
  
    betamin = []

    betastart = 50
        
    for beta in range(betastart, d):
        l = simulate_pnjBKZ(l, beta, 1, 1)

        proba = 1.
            
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(
                2**(2 * l[i]))
        
        G1, B1 = bkz_cost(d,beta,cost_model=cost_model,ldc_param = ldc_param)

        if(not worst_case):
            #G1cum = log2(2**G1cum + ((2**G1) * remaining_proba * proba))
            GBKZ = log2(2**GBKZ + 2**G1)
            G1cum = log2(2**G1cum + ((2**GBKZ) * remaining_proba * proba))
        else:
            G1cum = log2(2**G1cum + 2**G1)
       
        cumulated_proba += remaining_proba * proba
        remaining_proba = 1. - cumulated_proba

        #d_svp prediction
        (G_sieve,B_sieve,dsvp,dsvp_prime) = d_svp_prediction(l,cumulated_proba, cost_model,ldc_param)


        #print(beta, proba, l[d-beta], cumulated_proba, G_sieve,  dsvp, dsvp_prime, dsvp-dsvp_prime, dim4free_wrapper(theo_dim4free_fun2,dsvp))

        
        G = log2(2**G1cum + 2**G_sieve)
        B = max(B1cum, B_sieve)
       
        if remaining_proba < .001:
            if(not worst_case):
                G = log2(2**G + ((2**GBKZ) * remaining_proba))
            
        if(G!= float("inf") and B!= float("inf")):
            if(G < Gmin):
                Gmin, Bmin, G1min, avgbetamin, dsvpmin, dsvp_prime_min= G, B, G1cum, average_beta, dsvp, dsvp_prime
                betamin = list(range(betastart,beta+1))
           

        if verbose:
            print("β= %d, cum-pr=%.2e, pump-{%d,%d,%d},  G=%3.2f gate,  B=%3.2f bit"%(beta, cumulated_proba, d-dsvp, dsvp, dsvp-dsvp_prime, Gmin, Bmin), end="\r" if cumulated_proba < 1e-4 else "\n")
            
        if remaining_proba < .001:
            break
        
        
    if remaining_proba > .01:
        raise ValueError("This instance may be unsolvable")
       
    print()
    return betamin,G1min, dsvpmin, dsvp_prime_min, Gmin, Bmin
