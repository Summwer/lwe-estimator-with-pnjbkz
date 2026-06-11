load("../framework/utils.sage")
load("../framework/d_svp_prediction.sage")



def two_step_mode_estimation(l, dvol, dim_, verbose=False, cost_model=1,ldc_param = "AGPS20", cal_ee = "chi", worst_case = False, goal_min_cost = "gate_min", dsvp_model = 2, success_rate=0.999):
    if(cal_ee == "chi"):
        return proba_two_step_mode_estimation(l, dvol, verbose=verbose, cost_model=cost_model, ldc_param = ldc_param, worst_case = worst_case, goal_min_cost = goal_min_cost, success_rate=success_rate)
    if(cal_ee == "avg_sigma"):
        return stev_two_step_mode_estimation(dvol, dim_, verbose=verbose, cost_model=cost_model, ldc_param = ldc_param, goal_min_cost = goal_min_cost, success_rate=success_rate)

def stev_two_step_mode_estimation(dvol, d, verbose=False, cost_model=1,ldc_param = "AGPS20", goal_min_cost = "gate_min", success_rate=0.999):
    
    """
    LWE estimation: Simplified progressive BKZs + Pump through average sigma value, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
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
    for beta in range(betastart, min(d,1000)):
        #l = simulate_pnjBKZ(l, beta, 1, 1)
        delta = compute_delta(beta)
        #print(beta, delta, dvol, d)
        l = [log(bkzgsa_gso_len(dvol, i, d, delta=delta)) / log(2.) for i in range(d)]
        G1, B1 = bkz_cost(d,beta,ldc_param = ldc_param)
        G1 += log2(C)

        #d_svp prediction
        (dsvp, dsvp_prime, Gpump, Bpump) = stev_d_svp_prediction(l,cost_model,ldc_param)

        G = log2_add(G1 , Gpump)
        B = max(B1, Bpump)

        if(goal_min_cost == "gate_min" and G < Gmin):
            Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = G, B, dsvp, dsvp_prime, G1
            betamin = list(range(betastart,beta+1))
        if(goal_min_cost == "gate_RAM_min" and G + B < Gmin + Bmin):
            Gmin, Bmin, dsvpmin, dsvp_prime_min, G1min = G, B, dsvp, dsvp_prime, G1
            betamin = list(range(betastart,beta+1))
        
       
        if verbose:
            print("β= %d, pump-{%d,%d,%d},  G=%3.2f gate,  B=%3.2f bit"%(beta, d-dsvp, dsvp, dsvp-dsvp_prime, Gmin, Bmin), end = "\r" if beta%50 != 0 else "\n")
    print()
    return betamin, G1min, dsvpmin, dsvp_prime_min, Gmin, Bmin


def proba_two_step_mode_estimation(l, logvol, betastart = 50, verbose=False, cost_model=1,ldc_param = "AGPS20", worst_case = False, goal_min_cost = "gate_min", b = None, success_rate=0.999):
    """
    LWE estimation: Simplified progressive BKZs + Pump throught chi_square distribution, cumulate cost directly.
    Computes the probabilistic cumulated cost value for given gs-lengths.
    :l: log2(||b_i^*||), i = 0,...,d-1
    :cost_model: 1: gate model
                 2: sec model with threads=32, gpus = 2 
    :progressieve_sieve: True: progressieve sieve
                         False: normal sieve
    :goal_min_cost: "gate_min": find the minimal gates cost in two-step
                    "gate_RAM_min": find the minimal (gates+RAM) cost in two-step
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
    G1cum,B1cum = -1000., -1000.
    GBKZ = -1000.
  
    betamin = []
    #Gcums = [0.]
    #cumulated_probas = [0.]
    #betas = [0]
    Gs = []
    if(cost_model == 2):
        betastart = 50
    print_Gcums = []

    betastart  = 50
    dsvp = 0
    dsvp_prime = 0

    # default g6k strategy 
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

    #blocksizes = list(range(10, 50)) + [b-20, b-17] + list(range(b - 14, d, 2))
    #blocksizes = [ b for b in blocksizes if b >= 50]
    #trivial strategy 
    blocksizes = list(range(betastart,d))


    for beta in blocksizes:

        l = simulate_pnjBKZ(l, beta, 1, 1)
        #print("beta = %d, current slope = %.4f\n"  %(beta,get_current_slope(l,0,d)))

        proba = 1.
            
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(
                2**(2 * l[i]))
        
        G1, B1 = bkz_cost(d,beta,cost_model=cost_model,ldc_param = ldc_param)
        

        if(not worst_case):   
            #G1cum = log2(2**G1cum + ((2**G1) * remaining_proba * proba))
            GBKZ = log2(2**GBKZ + 2**G1)
            G1cum = log2(2**G1cum + ((2**GBKZ) * remaining_proba * proba))
            B1cum = log2(2**B1cum + ((2**B1) * remaining_proba * proba))
            #B1cum = max(B1cum, B1)
    
        else:
            G1cum = log2(2**G1cum + 2**G1)
            B1cum = max(B1cum, B1)

        
        cumulated_proba += remaining_proba * proba
        remaining_proba = 1. - cumulated_proba

        if(cumulated_proba >= success_rate):
            if(not worst_case):
                G = log2(2**G + ((2**GBKZ) * remaining_proba )) - log2(cumulated_proba)
                B = log2(2**B + ((2**B1) * remaining_proba )) - log2(cumulated_proba)
              
        else:
            (G_sieve,B_sieve,dsvp,dsvp_prime, Gcums1, cumulated_probas1) = d_svp_prediction(l, cumulated_proba, cost_model,ldc_param,worst_case,GBKZ=GBKZ, Gcums = [], cumulated_probas= [], success_rate=success_rate)
        
      
            
            G = log2_add(G1cum , G_sieve) - log2(success_rate)
           

            if(worst_case):
                B = max(B1cum, B_sieve)
            else:
                B = log2_add(B1cum , B_sieve) - log2(success_rate)
        
        Gs.append(G)         
        print_Gcums.append(G)

        
        if(G!= float("inf") and B!= float("inf")):
            if(goal_min_cost == "gate_min" and G < Gmin):
                Gmin, Bmin, G1min, avgbetamin, dsvpmin, dsvp_prime_min, Gcumsmin,cumulated_probasmin = G, B, G1cum, average_beta, dsvp, dsvp_prime,Gcums1,cumulated_probas1
                
                betamin = list(range(betastart,beta+1))
                if verbose:
                    if(cost_model == 1):
                        print("β= %d, G = %3.2f log2(gate), G1 = %3.2f log2(gate) B =%3.2f log2(bit), GBKZ = %3.2f log2(gate), cum-pr=%.4e"%(beta, G, G1cum, B, GBKZ, cumulated_proba), 
                            end="\r" if cumulated_proba < 1e-5 else "\n")
                    if(cost_model == 2): #G1 = %3.2f log2(sec), Gpump = %3.2f log2(sec), # G1cum, G_sieve
                        print("β= %d, G = %3.2f log2(sec),  B =%3.2f log2(bit), cum-pr=%.4e"%(beta, G, B, cumulated_proba), 
                                    end="\r" if cumulated_proba < 1e-5 else "\n")
                    
            
            if(goal_min_cost == "gate_RAM_min" and G + B < Gmin + Bmin):
                Gmin, Bmin, G1min, avgbetamin, dsvpmin, dsvp_prime_min, Gcumsmin,cumulated_probasmin= G, B, G1cum, average_beta, dsvp, dsvp_prime,Gcums1,cumulated_probas1
                betamin = [ b for b in blocksizes if b >= 50 and b <= beta]
                

        if(cumulated_proba >= success_rate):
            #prob_diff = 1 - cumulated_proba
            #log_prob_weight = log2(rp) + log2(prob_diff)
            #log_step_total_cost = log2_add(Gpump, GBKZ) + log_prob_weight
            #log_mem_step_cost = Bsieve + log_prob_weight
            #G_cum = log2_add(G_cum, log_step_total_cost)
            #B_cum = log2_add(B_cum, log_mem_step_cost)
            break

       

    #if(verbose):
        #print()
        #print("Gcumsmin: ", Gcumsmin)
        #print("cumulated_probasmin: ", cumulated_probasmin)
        #print("betas: ", list(range(betastart, beta + 1)))
        #print("Gs: ", Gs )
    print("slope after PnJBKZs: ", get_current_slope(l, start_row=0, stop_row=len(l)))

    return betamin,G1min, dsvpmin, dsvp_prime_min, Gmin, Bmin
