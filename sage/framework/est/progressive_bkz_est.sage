load("../framework/utils.sage")

def pro_bkz_est(l, verbose=False, cost_model = 1, worst_case = False,  ldc_param = "AGPS20" , cumG = False):
    if(not cumG):
        avgbeta = leaky_lwe_est(l, verbose=verbose)
        G, B = bkz_cost(len(l),avgbeta,ldc_param = ldc_param)
        G += log2(C)
        return avgbeta, G, B
    if(cumG):
        return pro_bkz_cumG(l, verbose=verbose, cost_model = cost_model, worst_case = worst_case,ldc_param = ldc_param)

def pro_bkz_cumG(l, verbose=False, cost_model = 1, worst_case = False,ldc_param = "AGPS20"):
    """
    LWE estimation: Simplified progressive BKZs + Pump, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    """
    d = len(l)
    G1cum,B1cum = 0.,0.
    GBKZ = 0.

    # Keep increasing beta to be sure to catch the second intersection
    
    remaining_proba = 1.
    average_beta = 0.
    cumulated_proba = 0.
    
    prob_scumGs = []
    probs = []
    for beta in [x  for x in range(50, d)]:

        #l = simBKZ(l, beta, 1)   
        l = simulate_pnjBKZ(l, beta, 1, 1)#simulate_pnjBKZ(log_GS_lengths, beta, loop, jump)
            
        proba = 1.
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(2**(2 * l[i]))

    
        average_beta += beta * remaining_proba * proba
                
        G1, B1 = bkz_cost(d,beta,cost_model=cost_model, ldc_param = ldc_param)

        if(not worst_case):
            GBKZ = log2(2**GBKZ+ 2**G1)
            G1cum = log2(2**G1cum + ((2**GBKZ) * remaining_proba * proba))
            B1cum = log2(2**B1cum + ((2**B1) * remaining_proba * proba))
            #G1cum = log2(2**G1cum * ((2**(GBKZ * remaining_proba * proba))))
            #B1cum = log2(2**B1cum * ((2**(B1 * remaining_proba * proba))))
        else:
            G1cum = log2(2**G1cum + 2**G1)
            B1cum = max(B1cum,B1)
        

        cumulated_proba += remaining_proba * proba
        remaining_proba = 1. - cumulated_proba 


        if(cumulated_proba > 1e-4):
            prob_scumGs.append(G1cum)
            probs.append(cumulated_proba)
      

        if verbose:
            print("β= %d, G = %3.2f log2(gate), B =%3.2f log2(bit), cum-pr=%.4e, pr=%.4e"%
                        (beta, G1cum, B1cum, cumulated_proba,proba), 
                        end="\r" if cumulated_proba < 1e-4 else "\n")

        if remaining_proba < .001:
            average_beta += beta * remaining_proba
            #G1, B1 = bkz_cost(d,beta,cost_model=cost_model)
            if(not worst_case):
                G1cum = log2(2**G1cum + ((2**GBKZ)) * remaining_proba)
                B1cum = log2(2**B1cum + ((2**B1 * remaining_proba)))
                #G1cum = log2(2**G1cum * ((2**(GBKZ * remaining_proba))))
                #B1cum = log2(2**B1cum * ((2**(B1 * remaining_proba))))
            break

    if remaining_proba > .01:
        raise ValueError("This instance may be unsolvable")


    print("probG1cum: ", prob_scumGs)
    print("cum_probs: ", probs)
    return average_beta, G1cum, B1cum


def leaky_lwe_est(l, verbose=False):
    """
    LWE estimation: Simplified progressive BKZs + Pump, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    """
    d = len(l)
    # Keep increasing beta to be sure to catch the second intersection
    
    remaining_proba = 1.
    average_beta = 0.
    cumulated_proba = 0.
    
    
    for beta in [x  for x in range(50, d)]:
        #l = simBKZ(l, beta, 1)   
        l = simulate_pnjBKZ(l, beta, 1, 1)
        proba = 1.
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(2**(2 * l[i]))

        average_beta += beta * remaining_proba * proba
        cumulated_proba += remaining_proba * proba
        remaining_proba = 1. - cumulated_proba 
      

        if verbose:
            print("β= %d, G = %3.2f log2(gate), B =%3.2f log2(bit), cum-pr=%.4e, pr=%.4e"%
                        (beta, bkz_cost(d, average_beta)[0], bkz_cost(d, average_beta)[1], cumulated_proba, proba), 
                        end="\r" if cumulated_proba < 1e-4 else "\n")

        if remaining_proba < .001:
            average_beta += beta * remaining_proba
            break

    if remaining_proba > .01:
        raise ValueError("This instance may be unsolvable")

    return average_beta

