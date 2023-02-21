load("../framework/utils.sage")




def bkz_only(d, logvol, l, verbose=False, cost_model = 1):
    """
    Computes the beta value for given dimension and volumes
    It is assumed that the instance has been normalized and sphericized, 
    i.e. that the covariance matrices of the secret is the identity
    :d: integer
    :vol: float
    """

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


    G1cum,B1cum = 0.,0.

    # Keep increasing beta to be sure to catch the second intersection
    
    remaining_proba = 1.
    average_beta = 0.
    cumulated_proba = 0.

    max_tours = 20

        
    for beta in [ceil(bbeta)  for _ in range(max_tours)]:
            
        l = simulate_pnjBKZ(l, beta, 1, 1)#simulate_pnjBKZ(log_GS_lengths, beta, loop, jump)
            
      
        proba = 1.
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(
            2**(2 * l[i]))

        average_beta += beta * remaining_proba * proba
        G1, B1 = bkz_cost(d,beta,cost_model=cost_model)
        G1cum = log2(2**G1cum + ((2**G1) * remaining_proba * proba))
        B1cum = max(B1cum,B1)

        cumulated_proba += remaining_proba * proba
        remaining_proba = 1. - cumulated_proba

        if verbose:
            print("β= %d,\t pr=%.4e, \t cum-pr=%.4e \t rem-pr=%.4e"%
                        (beta, proba, cumulated_proba, remaining_proba), 
                        end="\r" if cumulated_proba < 1e-4 else "\n")
        

        if remaining_proba < .001:
            average_beta += beta * remaining_proba * proba
            G1, B1 = bkz_cost(d,beta,cost_model=cost_model)
            G1cum = log2(2**G1cum + ((2**G1) * remaining_proba))
            B1cum = max(B1cum,B1)
            break

    if remaining_proba > .01:
        print()
        raise ValueError("This instance may be unsolvable")

    
    
    return average_beta , G1cum, B1cum

