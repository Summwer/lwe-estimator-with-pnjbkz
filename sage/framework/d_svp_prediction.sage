def d_svp_prediction(l, cumulated_proba,cost_model,ldc_param,worst_case,dsvp_model = 1):
    if(dsvp_model == 1):
        return d_svp_prediction_model1(l, cumulated_proba,cost_model,ldc_param,worst_case)
    else:
        return d_svp_prediction_model2(l, cumulated_proba,cost_model,ldc_param)

#compute pump cost in succ-fail probability
def d_svp_prediction_model2(l, cumulated_proba,cost_model,ldc_param):
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
           
    #predict dimension of last sieve: progressive sieve
    p = deepcopy(cumulated_proba)
    rp = 1. - p
    avgdsvp = 0.
    G_cum, B_cum = 0.,0.
    Gpump = 0.
    pre_psvp2 = 0.

    dsvp = 50
    flag1 = false
    flag2 = false
    while(dsvp <= d):
        #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic(l_[d-dsvp:])
        psvp1 = 1. * chisquared_table[dsvp].cum_distribution_function(gh)
        psvp2 = 1. * chisquared_table[dsvp].cum_distribution_function(4/3.*gh)
        
        if(not flag2):
            Gpump, Bpump = pump_cost(d,dsvp,cost_model=cost_model, ldc_param = ldc_param)
            G_cum = log2(pow(2,G_cum)+pow(2,Gpump) * (psvp2-pre_psvp2))
            B_cum = log2(pow(2,B_cum)+pow(2,Bpump) * (psvp2-pre_psvp2))
            #B_cum = max(B_cum,Bpump)
       
        if(psvp1 > 0.999 and not flag1):
            dsvp1 = dsvp
            flag1 = True
    
        if(psvp2 > 0.999 and not flag2):
            dsvp2 = dsvp
            flag2 = True
        
        if(flag1 and flag2):
            return  G_cum,B_cum,dsvp1, dsvp2
    
        pre_psvp2 = psvp2
        dsvp+=1 







#compute pump cost in cumulated probability
def d_svp_prediction_model1(l, cumulated_proba,cost_model,ldc_param,worst_case):
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
    #predict dimension of last sieve: progressive sieve
    p = deepcopy(cumulated_proba)
    rp = 1. - p
    avgdsvp = 0.
    avgG_sieve,avgB_sieve = 0.,0.
    pre_psvp = 0. 
    for dsvp in range(50, d):
        psvp = 1.
        #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic(l_[d-dsvp:])
        psvp *= chisquared_table[dsvp].cum_distribution_function(gh)

        dsvp_prime = floor(dsvp - dim4free_wrapper(theo_dim4free_fun2,dsvp))
        Gpump, Bpump = pump_cost(d,dsvp_prime,cost_model=cost_model,ldc_param=ldc_param)
            
        avgdsvp += dsvp * rp * psvp
        if(not worst_case):
            avgG_sieve = log2(2**avgG_sieve+(2**Gpump) * rp * psvp)
            #avgB_sieve = log2(2**avgB_sieve + (2**Bpump) * rp * psvp)
            avgB_sieve = max(Bpump,avgB_sieve)
        else:
            G_cum = log2(pow(2,G_cum)+pow(2,Gpump) * (psvp-pre_psvp))
            avgB_sieve = max(Bpump,avgB_sieve)

        p += rp * psvp
        rp = 1. - p
        #print(dsvp, gh)
        
        if rp < 0.001:
            avgdsvp += dsvp * rp #Avoid too small of dsvp
            avgG_sieve = log2(2**avgG_sieve + ((2**Gpump) * rp))        
            return (avgG_sieve,avgB_sieve,dsvp,avgdsvp)
        else:
            if(1-psvp < 0.001):
                return (avgG_sieve,avgB_sieve,dsvp,avgdsvp)
        pre_psvp = psvp
            
    return (G_sieve,B_sieve,dsvp,dsvp)



#find the minimal dimension for last pump using the average sigma formula
def stev_d_svp_prediction(l,cost_model,ldc_param):
    """
    Dimension of sieve/progressive sieve chosen to find target vector.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log(||b_i^*||), i = 0,...,d-1
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
        #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic(l_[d-dsvp:])
        
        # dsvp_prime = dsvp - d4f(dsvp) <= sqrt(4/3.) ||pi_(d-dsvp)(b_(d-dsvp))||
        if(not flag2 and dsvp <= 4/3. * gh): 
            dsvp_prime = dsvp
            flag2 = True
        if(not flag1 and dsvp <= gh):
            dsvp_full = dsvp
            flag1 = True
        if(flag1 and flag2):
            Gpump, Bpump = pump_cost(d,dsvp_prime,cost_model=cost_model, ldc_param = ldc_param)
            return dsvp, dsvp_prime, Gpump, Bpump
        dsvp+=1