load("../framework/utils.sage")
load("../framework/simulator/pnjbkz_simulator.sage")
load("../framework/simulator/pump_simulator.sage")
load("../framework/cost.sage")


def simulate_pump(l,up_dsvp, cumulated_proba,progressive_sieve = False,cost_model=1, worst_case = False):
    l0 = deepcopy(l)
    d = len(l)
    remaining_proba = 1. - cumulated_proba
    if(cumulated_proba>= 1. or up_dsvp < 50):
        return (0.,0,0.,0,0,l)
    if not progressive_sieve:
        G_sieve, B_sieve = float("inf"), float("inf")
        for dsvp in range(50, up_dsvp+1):
            psvp = 1.
            #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh

            l_ = pump_simulator(l0,d,d-dsvp)
            psvp *= chisquared_table[dsvp].cum_distribution_function(2**(2 * l_[d-dsvp]))
            
            p = cumulated_proba + remaining_proba * psvp
            rp = 1 - p

            if rp < 0.001:
                G_sieve, B_sieve = pump_cost(d,dsvp,cost_model=cost_model)
                return (G_sieve,B_sieve,dsvp+dsvp*rp,dsvp,1.,l_)
        G_sieve, B_sieve = pump_cost(d,dsvp,cost_model=cost_model)
        return (G_sieve,B_sieve,dvsp,dsvp,p,l_)    
       
        
    else:           
        #predict dimension of last sieve: progressive sieve
        p = deepcopy(cumulated_proba)
        rp = 1. - p
        avg_d_svp = 0.
        avgG2,avgB2 = 0.,0.
        Gpump = 0.
        pre_psvp = 0.
        for dsvp in range(50, up_dsvp+1):
            
            psvp = 1.
            #2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
            l_ = pump_simulator(l0,d,d-dsvp)
            psvp *= chisquared_table[dsvp].cum_distribution_function(2**(2 * l_[d-dsvp]))

            G_sieve, B_sieve = pump_cost(d,dsvp,cost_model=cost_model)
        
            avg_d_svp += dsvp * rp * psvp
            Gpump = log2(2**Gpump+2**G_sieve)
            if(not worst_case):
                avgG2 = log2(2**avgG2+(2**Gpump) * rp * psvp)
            else:
                avgG2 = log2(2**avgG2+(2**Gpump) * (1 - pre_psvp))
            
            avgB2 = max(B_sieve,avgB2)

            p += rp * psvp
            rp = 1. - p
            #print(dsvp, gh)
        
            if(not worst_case):
                if rp < 0.001:
                    #raise ""
                    #print(rp,avg_d_svp,dsvp * rp)
                    avg_d_svp += dsvp * rp #Avoid too small of dsvp
                    avgG2 = log2(2**avgG2 + ((2**Gpump) * rp))
                    return (avgG2,avgB2,avg_d_svp,dsvp,1.,l_)
            else:
                if(1-psvp < 0.001):
                    return (avgG_sieve,avgB_sieve,avgdsvp,dsvp,1.,l_)
            
            pre_psvp = psvp
                
    return (G_sieve,B_sieve,avg_d_svp,dsvp,p,l_)

#LWE estimation: Simplified progressive BKZs + Pump
def default_g6k_est( d, logvol, b, l, verbose=False, progressive_sieve = True, cost_model=1, worst_case = False):
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
            return 9999, 0
            
        ddelta = compute_delta(bbeta) * margin**(1. / d)
        if prev_margin is not None and interpolate:
            beta_low = log(margin) / (log(margin) - log(prev_margin))
        else:
            beta_low = 0
        assert beta_low >= 0
        assert beta_low <= 1

        bbeta = bbeta - beta_low

        b = ceil(bbeta)



    target_norm = d / 1.5


    
    DDGR20 = False
    bbeta = None
    pprev_margin = None
    Gmin, Bmin, avgbetamin, dsvpmin = float("inf"), float("inf"), d, d

    # Keep increasing beta to be sure to catch the second intersection
    if(True):
        remaining_proba = 1.
        average_beta = 0.
        cumulated_proba = 0.
        Gcum,Bcum = 0.,0.
        G = 0.
        

        betamin = []

        blocksizes = list(range(10, 50)) + [b-20, b-17] + list(range(b - 14, b + 25, 2))
        blocksizes = [_ for _ in blocksizes if _ <= d and _ >=50]
        for beta in blocksizes:

            l = simulate_pnjBKZ(l, beta, 1, 1)

            proba = 1.
            
            i = d - beta
            proba *= chisquared_table[beta].cum_distribution_function(
                2**(2 * l[i]))
            
            average_beta += beta * remaining_proba * proba

            G1, B1 = bkz_cost(d,beta,1,cost_model=cost_model)

            if(not worst_case):
                G = log2(2**G +2**G1)
                Gcum = log2(2**Gcum + ((2**G) * remaining_proba * proba))
            else:
                Gcum = log2(2**Gcum + 2**G1)
            
            Bcum = max(Bcum,B1)

            cumulated_proba += remaining_proba * proba
            remaining_proba = 1. - cumulated_proba

            #d_svp prediction
            
            G_sieve,B_sieve = 0., 0.
            
            n_max = int(58 + 2.85 * G) #G_BKZ
            #n_max = int(53 + 2.85 * Gcum)

            for n_expected in range(2, d-2):
                x = target_norm * n_expected/(1.*d)
                l_ = [2*_*log(2.) for _ in l]
                if 4./3 * gaussian_heuristic(l_[d-n_expected:]) > x:
                    break
                   

            print("Without otf, would expect solution at pump-%d. n_max=%d in the given time." % (n_expected, n_max)) # noqa
            flag = true
            if n_expected >= n_max - 1:
                flag = false


            llb = d - beta
            l_ = [2*_*log(2.) for _ in l]
            while gaussian_heuristic(l_[llb:]) < target_norm * (d - llb)/(1.*d): # noqa
                llb -= 1
                if llb < 0:
                    break
            
            if(flag):
                f = dims4free(d-llb)
                print("Starting svp pump_{%d, %d, %d}, n_max = %d, Tmax= %.2f sec" % (llb, d-llb, f, n_max, G1))
                (G_sieve,B_sieve,avg_d_svp,dsvp,cumulated_proba,l) = simulate_pump(l,d-llb, cumulated_proba,progressive_sieve = progressive_sieve ,cost_model=cost_model)
                remaining_proba = 1. - cumulated_proba
      

            G = log2(2**Gcum + 2**G_sieve)
            B = max(Bcum, B_sieve)

            if verbose:
                print("β= %d, cum-pr=%.2e,  G=%3.2f gate,  B=%3.2f bit"%
                        (beta, cumulated_proba,  G, B), 
                        end="\n")
    
            if remaining_proba < .001:
                average_beta += beta * remaining_proba 
                if(not worst_case):
                    G = log2(2**Gcum + ((2**G) * remaining_proba))
                break
        
        if remaining_proba > .01:
            raise ValueError("This instance may be unsolvable")
    
        betamin = list(range(50,beta+1))

        return betamin, G, B
