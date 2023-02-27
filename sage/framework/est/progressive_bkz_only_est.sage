load("../framework/utils.sage")

def progressive_bkz_only(l, verbose=False, cost_model = 1):
    """
    LWE estimation: Simplified progressive BKZs + Pump, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    """
    d = len(l)
    G1cum,B1cum = 0.,0.

    # Keep increasing beta to be sure to catch the second intersection
    
    remaining_proba = 1.
    average_beta = 0.
    cumulated_proba = 0.
    
    for beta in [x  for x in range(50, d)]:

        #l = simBKZ(l, beta, 1)   
        l = simulate_pnjBKZ(l, beta, 1, 1)#simulate_pnjBKZ(log_GS_lengths, beta, loop, jump)
            
      
        proba = 1.
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(2**(2 * l[i]))

    
        average_beta += beta * remaining_proba * proba
                
        G1, B1 = bkz_cost(d,beta,cost_model=cost_model)
        G1cum = log2(2**G1cum + ((2**G1) * remaining_proba * proba))
        
        B1cum = max(B1cum,B1)

        cumulated_proba += remaining_proba * proba
        #remaining_proba = 1. - cumulated_proba 
        remaining_proba *= 1. - proba

        if verbose:
            print("β= %d,\t pr=%.4e, \t cum-pr=%.4e \t rem-pr=%.4e"%(beta, proba, cumulated_proba, remaining_proba), end="\r" if cumulated_proba < 1e-4 else "\n")

        if remaining_proba < .001:
            average_beta += beta * remaining_proba
            G1, B1 = bkz_cost(d,beta,cost_model=cost_model)
            G1cum = log2(2**G1cum + ((2**G1) * remaining_proba))
            B1cum = max(B1cum,B1)
            break

    if remaining_proba > .01:
        raise ValueError("This instance may be unsolvable")

    return average_beta, G1cum, B1cum

