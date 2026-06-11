load("../framework/utils.sage")
load("../framework/simulator/pnjbkz_simulator.sage")
load("../framework/simulator/pump_simulator.sage")
load("../framework/cost.sage")


def simulate_pump(l,up_dsvp, cumulated_proba,cost_model=1, worst_case = False, sigma = 1., ldc_param = "AGPS20", GBKZ = -1000., success_rate=0.999):
    l0 = deepcopy(l)
    d = len(l)
    remaining_proba = 1. - cumulated_proba
    if(cumulated_proba>= success_rate or up_dsvp < 50):
        return (0.,0,0.,0,0,l)
               
    #predict dimension of last sieve: progressive sieve
    p = deepcopy(cumulated_proba)
    rp = 1. - p
    avg_d_svp = 0.
    avgG2,avgB2 = 0.,0.
    Gpump = 0.
    pre_psvp = 0.
    G_cum, B_cum = -1000.,-1000.
    for dsvp in range(30, up_dsvp+1):   
        #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        l_ = pump_simulator(l0,d,d-dsvp)
        #psvp *= chisquared_table[dsvp].cum_distribution_function(2**(2 * (l_[d-dsvp]- log(sigma)/log(2))))
        psvp = 1. * chisquared_table[dsvp].cum_distribution_function(2**(2 * (l_[d-dsvp])))
        
        if(cost_model == 2):
            dsvp_prime = floor(dsvp - dim4free_wrapper(default_dim4free_fun,dsvp))
        else:
            dsvp_prime = floor(dsvp - theo_dim4free_fun2(dsvp))

        Gpump, Bsieve = pump_cost(dsvp_prime,cost_model=cost_model,ldc_param = ldc_param)
        
        
        if(pre_psvp >= psvp ):
            continue

        G_cum = log2(pow(2,G_cum)+ (pow(2,Gpump)+pow(2,GBKZ)) * rp * (psvp-pre_psvp))
        B_cum = log2(pow(2,B_cum)+pow(2,Bsieve) * rp * (psvp-pre_psvp))
        cum_p = p + rp * psvp
        if(cum_p >= success_rate):
            break
            
        pre_psvp = psvp
                
    return (G_cum,Gpump, B_cum,dsvp,dsvp_prime,cum_p, l_)

#LWE estimation: Simplified progressive BKZs + Pump
def default_g6k_est( d, logvol, b, l, verbose=False,  cost_model=1, worst_case = False,threads = 32, sigma = 1,ldc_param = "AGPS20", success_rate=0.999):
    """
    Computes the beta value for given dimension and volumes
    It is assumed that the instance has been normalized and sphericized, 
    i.e. that the covariance matrices of the secret is the identity
    :d: integer
    :vol: float
    """
    


    if(b == None):
        #select b

        bbeta = None
        pprev_margin = None
        interpolate = True
        
        for beta in range(2, d):
            lhs = RR(sqrt(beta))
            rhs = bkzgsa_gso_len(logvol, d - beta, d, beta=beta)

            if lhs < rhs and bbeta is None:
                margin = rhs / lhs
                prev_margin = pprev_margin
                bbeta = beta

            if lhs > rhs:
                bbeta = None
            pprev_margin = rhs / lhs

        if bbeta is None:
            return 9999, 0, 0, 0, 0
            
        ddelta = compute_delta(bbeta) * margin**(1. / d)
        if prev_margin is not None and interpolate:
            beta_low = log(margin) / (log(margin) - log(prev_margin))
        else:
            beta_low = 0
        assert beta_low >= 0
        assert beta_low <= 1

        bbeta = bbeta - beta_low

        b = ceil(bbeta)

    goal_margin = 1.5
    target_norm = goal_margin * (sigma**2) * (d-1) + 1
    
    bbeta = None
    pprev_margin = None
    Gmin, Bmin, avgbetamin, dsvpmin = float("inf"), float("inf"), d, d

    # Keep increasing beta to be sure to catch the second intersection
    if(True):
        remaining_proba = 1.
        average_beta = 0.
        cumulated_proba = 0.
        G,B = 0.,0.
        
        betamin = []

        #blocksizes = list(range(10, 50)) + [b-20, b-17] + list(range(b - 14, b + 25, 2))
        blocksizes = list(range(10, 50)) + [b-20, b-17] + list(range(b - 14, d, 2))

        blocksizes = [ b for b in blocksizes if b >= 50]
   
        dsvp, dsvp_prime = 0, 0
        G1cum,B1cum = -1000., -1000.
        GBKZ = -1000.
        for beta in blocksizes:

            l = simulate_pnjBKZ(l, beta, 1, 1)

            proba = 1.
            
            i = d - beta

            
            #proba *= chisquared_table[beta].cum_distribution_function(
            #    2**(2 * (l[i]- log(sigma)/log(2))))

            proba *= chisquared_table[beta].cum_distribution_function(
                2**(2 * (l[i])))
            
            average_beta += beta * remaining_proba * proba

            G1, B1 = bkz_cost(d,beta,cost_model=cost_model,ldc_param = ldc_param)


            if(not worst_case):   
                #G1cum = log2(2**G1cum + ((2**G1) * remaining_proba * proba))
                GBKZ = log2(2**GBKZ + 2**G1)
                G1cum = log2(2**G1cum + ((2**GBKZ) * remaining_proba * proba))
                #B1cum = log2(2**B1cum + ((2**B1) * remaining_proba * proba))
                B1cum = max(B1cum, B1)
        
            else:
                G1cum = log2(2**G1cum + 2**G1)
                B1cum = max(B1cum, B1)

            cumulated_proba += remaining_proba * proba
            remaining_proba = 1. - cumulated_proba

            if remaining_proba < 1 - success_rate:
                G = G1cum
                B = B1cum
                if(not worst_case):
                    G = log2(2**G + ((2**GBKZ) * remaining_proba))

                break
     

            if verbose:
                if(cost_model == 1):
                    print("β= %d, slope=%f,  G=%3.2f log2(gate),  B=%3.2f log2(bit)"%
                        (beta, get_current_slope(l,0,d),  G, B))
                if(cost_model == 2):
                    print("β= %d, slope=%f,  G=%3.2f log2(sec), walltime = %3.2f sec, B=%3.2f log2(bit)"%
                        (beta, get_current_slope(l,0,d),  G, 2**G, B))
                

            

            #d_svp prediction
            
            Gpump,Bpump = 0., 0.
            
            n_max = int(58 + 2.85 * (GBKZ + log(threads,2))) #32: threads number
            
            l_ = [2.*_*log(2.) for _ in l]
        
            for n_expected in range(50, d-2):
                x = (target_norm / goal_margin) * n_expected/(1.*d)
                #if(beta == 112):
                #    print(n_expected,x,4./3 * gaussian_heuristic(l_[d-n_expected:]))
                if 4./3 * gaussian_heuristic(l_[d-n_expected:]) > x:
                    break
 

            print("Without otf, would expect solution at pump-%d. n_max=%d in the given time." % (n_expected, n_max)) # noqa
            flag = true
            if n_expected >= n_max - 1:
                continue

            n_max += 1
        
            llb = d - beta
            while gaussian_heuristic(l_[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
                llb -= 1
                if llb <= 0:
                    break
            

            if(flag):
                f = max(0,d - llb - n_max)
                if(cost_model == 1):
                    print("Starting svp pump_{%d, %d, %d}, n_max = %d, Tmax= %.2f gate" % (llb, d-llb, f, n_max, 2**G1))
                elif(cost_model == 2):
                    print("Starting svp pump_{%d, %d, %d}, n_max = %d, Tmax= %.2f sec" % (llb, d-llb, f, n_max, 2**G1))
                dsvp, dsvp_prime = get_beta_from_sieve_dim(n_max,d,default_dim4free_fun)
                
                #dsvp = get_beta_from_sieve_dim(n_max,d,theo_dim4free_fun1)
                #dsvp = get_beta_from_sieve_dim(n_max,d,theo_dim4free_fun2)
                (G_sieve, Gpump, B_sieve,dsvp,dsvp_prime,cum_p, l) = simulate_pump(l,dsvp, cumulated_proba,cost_model=cost_model, sigma = sigma,ldc_param = ldc_param, GBKZ=GBKZ, success_rate=success_rate)
                remaining_proba = 1. - cum_p


                if(cumulated_proba >= success_rate):
                    G = G1cum
                    B = B1cum
                    break
                else:
                    G = log2(2**G1cum + 2**G_sieve)
                    GBKZ = log2(2**GBKZ + 2**Gpump)

                    if(worst_case):
                        B = max(B1cum, B_sieve)
                    else:
                        B = log2(2**B1cum + 2**B_sieve)

                    if(cum_p >= success_rate):
                        break

                
                
            
                if verbose:
                    if(cost_model == 1):
                        print("slope=%f,  G=%3.2f log2(gate),  B=%3.2f bit"%
                        (get_current_slope(l,0,d),  G, B))
                    if(cost_model == 2):
                        print("slope=%f,  G=%3.2f log2(sec), walltime = %3.2f sec, B=%3.2f bit"%
                        (get_current_slope(l,0,d),  G, 2**G, B))
    
      
            
            if remaining_proba < 1 - success_rate:
                #average_beta += beta * remaining_proba 
                #if(not worst_case):
                #    G = log2(2**G + ((2**GBKZ) * remaining_proba))
                #    B = log2(2**B1cum + 2**B_sieve)
                #else:
                #    G = log2(2**G1cum + 2**G_sieve)
                #    B = max(B1cum, B_sieve)
                break
        
        if remaining_proba > .001:
            print("This instance may be unsolvable.")
    
        betamin = [ b for b in blocksizes if b >= 50 and b <= beta]
        
        return betamin, G, B, dsvp, dsvp_prime
