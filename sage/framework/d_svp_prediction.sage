def d_svp_prediction(l, cumulated_proba,cost_model,ldc_param,worst_case,GBKZ = 0.,Gcums = [], cumulated_probas= []):
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
    if(cumulated_proba >= 0.999):
        return (-1000.,-1000.,0,0,Gcums,cumulated_probas)
           
    #predict dimension of last sieve: progressive sieve
    p = deepcopy(cumulated_proba)
    rp = 1. - p
    avgdsvp = 0.
    G_cum, B_cum = -1000.,-1000.
    pre_psvp = 0.

    for dsvp in range(30,d):
        #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic(l_[d-dsvp:])
        psvp = 1. * chisquared_table[dsvp].cum_distribution_function(gh)
        
        #dsvp_prime = floor(dsvp - dim4free_wrapper(theo_dim4free_fun2,dsvp))
        if(cost_model == 2):
            dsvp_prime = floor(dsvp - dim4free_wrapper(default_dim4free_fun,dsvp))
        else:
            dsvp_prime = floor(dsvp - theo_dim4free_fun2(dsvp))
        
        #if(cost_model == 1):
        #    Gsieve, Bsieve = sieve_cost(dsvp_prime,cost_model = cost_model, ldc_param=ldc_param)
        #    #print(cost_model, Gpump,2**Gpump,Gsieve,2**Gsieve,dsvp,dsvp_prime)
        #    Gpump = log2(2**Gpump + 2** Gsieve)
        #else:
        Gpump, Bsieve = pump_cost(dsvp_prime,cost_model = cost_model,ldc_param=ldc_param)


        if(pre_psvp >= psvp):
            continue
        
        #print(dsvp, G_cum, rp, psvp, pre_psvp)
        
        G_cum = log2(pow(2,G_cum)+ (pow(2,Gpump)+pow(2,GBKZ)) * rp * (psvp-pre_psvp))
        B_cum = log2(pow(2,B_cum)+pow(2,Bsieve) * rp * (psvp-pre_psvp))
       
        if(p + rp * psvp > 1e-4 ):
            Gcums.append(G_cum)
            cumulated_probas.append(p + rp * psvp)
            #cumulated_probas.append(psvp)

            #if(abs(Gcums[0] -131.4856784733793)<0.01):
            #    print("GBKZ = ", GBKZ, end=",")
            #    print("psvp = ", psvp, end=",")
            #    print("p + rp * psvp = ", p + rp * psvp)
            #    print("dsvp = ", dsvp)

        if(p + rp * psvp >= 0.999):
            break

        pre_psvp = psvp
        
    return  G_cum,B_cum,dsvp, dsvp_prime, Gcums,cumulated_probas




'''
def d_svp_prediction(l, cumulated_proba,cost_model,ldc_param,worst_case,GBKZ = 0.,Gcums = [], cumulated_probas= []):
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
    if(cumulated_proba >= 0.999):
        return (0,0,0,0,Gcums,cumulated_probas)
           
    #predict dimension of last sieve: progressive sieve
    p = deepcopy(cumulated_proba)
    rp = 1. - p
    avgdsvp = 0.
    G_cum, B_cum = 0.,0.
    Gpump = deepcopy(GBKZ)
    pre_psvp = 0.

    for dsvp in range(30,d):
        #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic(l_[d-dsvp:])
        psvp = 1. * chisquared_table[dsvp].cum_distribution_function(gh)
        
        #dsvp_prime = floor(dsvp - dim4free_wrapper(theo_dim4free_fun2,dsvp))
        dsvp_prime = floor(dsvp - dim4free_wrapper(default_dim4free_fun,dsvp))
        Gsieve, Bsieve = sieve_cost(dsvp_prime,cost_model = cost_model, ldc_param=ldc_param)
        #print(cost_model, Gpump,2**Gpump,Gsieve,2**Gsieve,dsvp,dsvp_prime)
        Gpump = log2(2**Gpump + 2** Gsieve)

        if(pre_psvp >= psvp):
            Gcums.append(G_cum)
            cumulated_probas.append(p)
            continue
        
        #print(dsvp, G_cum, rp, psvp, pre_psvp)
        
        G_cum = log2(pow(2,G_cum)+pow(2,Gpump) * rp * (psvp-pre_psvp))
        B_cum = log2(pow(2,B_cum)+pow(2,Bsieve) * rp * (psvp-pre_psvp))
       
        if(p + rp * psvp > 1e-4):
            Gcums.append(G_cum)
            #cumulated_probas.append(p + rp * psvp)

            p += rp * psvp
            rp = 1. - p

            cumulated_probas.append(p)

        if(p + rp * psvp >= 0.999):
            break

        pre_psvp = psvp
        
    return  G_cum,B_cum,dsvp, dsvp_prime, Gcums,cumulated_probas
'''

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
            Gpump, Bpump = pump_cost(dsvp_prime,cost_model=cost_model, ldc_param = ldc_param)
            return dsvp, dsvp_prime, Gpump, Bpump
        dsvp+=1