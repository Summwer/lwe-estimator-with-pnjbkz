

def d_svp_prediction(l, cumulated_proba,cost_model, progressive_sieve, worst_case):
    """
    Dimension of sieve/progressive sieve chosen to find target vector.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log(||b_i^*||), i = 0,...,d-1
    :cumulated_proba: current scuccess cumulated probability of gs-lengths in reduction 
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    :progressieve_sieve: True: progressieve sieve
                         False: normal sieve

    return value:
        dsvp/avgdsvp: average dimension value to sieve
        dsvp_r: the largest dimension value to sieve
        G_sieve: the cumulated time cost for sieve
        B_dsvp: the maximal memory cost for sieve

    """

    d = len(l)
    l_ = [2*_*log(2.) for _ in l]
    remaining_proba = 1. - cumulated_proba
    if(cumulated_proba>= 1.):
        return (0,0,0,0)
    if not progressive_sieve:
        G_sieve, B_sieve = float("inf"), float("inf")
        for dsvp in range(50, d):
            psvp = 1.
            #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
            gh = gaussian_heuristic(l_[d-dsvp:])
            psvp *= chisquared_table[dsvp].cum_distribution_function(gh)
            
            p = cumulated_proba + remaining_proba * psvp
            rp = 1. - p

            if rp < 0.001:
                dsvp += dsvp * rp #Avoid too small of dsvp
                break
        G_sieve, B_sieve = pump_cost(d,dsvp,cost_model=cost_model)
        
    else:           
        #predict dimension of last sieve: progressive sieve
        p = deepcopy(cumulated_proba)
        rp = 1. - p
        avgdsvp = 0.
        avgG_sieve,avgB_sieve = 0.,0.
        Gpump = 0.
        pre_psvp = 0.
        for dsvp in range(50, d):
            
            psvp = 1.
            #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
            gh = gaussian_heuristic(l_[d-dsvp:])
            psvp *= chisquared_table[dsvp].cum_distribution_function(gh)

            G_sieve, B_sieve = sieve_cost(d,dsvp,cost_model=cost_model)
            
            avgdsvp += dsvp * rp * psvp
            Gpump = log2(2**Gpump+2**G_sieve)
            if(not worst_case):
                #avgG_sieve = log2(2**avgG_sieve+(2**G_sieve) * rp * psvp)
                avgG_sieve = log2(2**avgG_sieve+(2**Gpump) * rp * psvp)
            else:
                avgG_sieve = log2(2**avgG_sieve+(2**G_sieve) * (1-pre_psvp))
            avgB_sieve = max(B_sieve,avgB_sieve)

            p += rp * psvp
            rp = 1. - p
            #print(dsvp, gh)
        
            if(not worst_case):
                if rp < 0.001:
                    avgdsvp += dsvp * rp #Avoid too small of dsvp
                    
                    avgG_sieve = log2(2**avgG_sieve + ((2**Gpump) * rp))
                    
                    return (avgG_sieve,avgB_sieve,avgdsvp,dsvp)
            else:
                if(1-psvp < 0.001):
                    return (avgG_sieve,avgB_sieve,avgdsvp,dsvp)
            pre_psvp = psvp
            
    return (G_sieve,B_sieve,dsvp,dsvp)
