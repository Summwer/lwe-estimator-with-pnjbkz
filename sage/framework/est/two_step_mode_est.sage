
load("../framework/utils.sage")
load("../framework/d_svp_prediction.sage")

def two_step_mode_estimation(l, verbose=False, cost_model=1, progressive_sieve = True, worst_case = False):
    """
    LWE estimation: Simplified progressive BKZs + Pump, cumulate cost directly.
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

    betastart = 79
        
    for beta in range(betastart, d):
        l = simulate_pnjBKZ(l, beta, 1, 1)

        proba = 1.
            
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(
                2**(2 * l[i]))
        

        G1, B1 = bkz_cost(d,beta,1,cost_model=cost_model)

        if(not worst_case):
            #G1cum = log2(2**G1cum + ((2**G1) * remaining_proba * proba))
            GBKZ = log2(2**GBKZ + 2**G1)
            G1cum = log2(2**G1cum + ((2**GBKZ) * remaining_proba * proba))
        else:
            G1cum = log2(2**G1cum + 2**G1)
            
        B1cum = max(B1cum,B1)

        cumulated_proba += remaining_proba * proba
        remaining_proba = 1. - cumulated_proba

        #d_svp prediction
        (G_sieve,B_sieve,dsvp,dsvp_max) = d_svp_prediction(l,cumulated_proba, cost_model,progressive_sieve, worst_case)

    
        
        G = log2(2**G1cum + 2**G_sieve)
        B = max(B1cum, B_sieve)
       
        if remaining_proba < .001:
            if(not worst_case):
                G = log2(2**G + ((2**GBKZ) * remaining_proba))
            
        if(G!= float("inf") and B!= float("inf")):
            if(G < Gmin):
                Gmin, Bmin, G1min, avgbetamin, dsvpmin, dsvpmaxmin = G, B, G1cum, average_beta, dsvp, dsvp_max
                betamin = list(range(betastart,beta+1))
           

        if verbose:
            print("β= %d, cum-pr=%.2e,  dsvp_min=%3d,  G=%3.2f gate,  B=%3.2f bit"%(beta, cumulated_proba, dsvpmin, Gmin, Bmin), end="\r")

        if remaining_proba < .001:
            break
        
        
    if remaining_proba > .01:
        raise ValueError("This instance may be unsolvable")
       
    print()
    return betamin,G1min, dsvpmin, dsvpmaxmin, Gmin, Bmin, None
