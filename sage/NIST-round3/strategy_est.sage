load("../framework/utils.sage")
load("../framework/d_svp_prediction.sage")
load("../framework/simulator/pnjbkz_simulator.sage")
load("../framework/cost.sage")
load("../framework/lwe_estimation.sage")

import matplotlib.pyplot as plt
import os
from random import randint



chisquared_table = {i: None for i in range(1000)}
for i in range(5000):
    chisquared_table[i] = RealDistribution('chisquared', i)

def strategy_est(l,strategy,cost_model,progressive_sieve, worst_case):
    costs = []
    d = len(l)
    cum_pr = 0.
    rem_pr = 1.
    G1cum, B1cum = 0., 0.
    slope = get_current_slope(l, 0, len(l))
    print(slope)
    for bs in strategy:
        (beta, jump, N) = bs
        if(beta>d):
            break
        for tours in range(N):
            l = simulate_pnjBKZ(l, beta, jump, 1)
            slope = get_current_slope(l, 0, len(l))
            print("slope: ", slope)
            i = d - beta
            proba = 1.*chisquared_table[beta].cum_distribution_function(2**(2 * l[i]))
            G1, B1 = bkz_cost(d,beta,J=jump,cost_model=cost_model)
            
            if(not worst_case):
                G1cum = log2(2**G1cum + ((2**G1) * rem_pr * proba))
            else:
                G1cum = log2(2**G1cum + 2**G1)
            B1cum = max(B1cum,B1)

            cum_pr += rem_pr * proba
            rem_pr *= 1. - proba

            (G_sieve,B_sieve,dsvp,dsvp_max) = d_svp_prediction(l,cum_pr,cost_model, progressive_sieve, worst_case)
            costs.append(log2(2**G1cum+2**G_sieve))
            print(beta,dsvp_max)
            if(cum_pr > 0.999):
                return costs
    if(strategy == []):
        (G_sieve,B_sieve,dsvp,dsvp_max) = d_svp_prediction(l,cum_pr,cost_model, progressive_sieve, worst_case)
        print(dsvp_max, G_sieve, B_sieve)
        return costs
   


worst_case = True
cost_model = 2
progressive_sieve = True


print("==========")
n = 40
alpha = 0.025
dim,dvol = gen_lwechal_instance(n, alpha)
dim = 172
dvol = 331.9735339
delta = compute_delta(2)
l = [log(bkzgsa_gso_len(dvol, i, dim, delta=delta)) / log(2) for i in range(dim)]
blocksize_strategy = [( 79,  8,  1),( 91,  8,  1),(112,  8,  1)]
costs = strategy_est(l,blocksize_strategy,cost_model,progressive_sieve, worst_case)



print("==========")
n = 45
alpha = 0.020
#dim,dvol = gen_lwechal_instance(n, alpha)
dim = 172
dvol = 331.9735339
delta = compute_delta(2)
l = [log(bkzgsa_gso_len(dvol, i, dim, delta=delta)) / log(2) for i in range(dim)]
blocksize_strategy = [( 79,  8,  1),( 91,  8,  1),(112,  8,  1)]
costs = strategy_est(l,blocksize_strategy,cost_model,progressive_sieve, worst_case)




"""


#n = 95
#alpha = 0.005
#dim,dvol = gen_lwechal_instance(n, alpha)
dim = 323
dvol = 836.9696072
delta = compute_delta(2)
l = [log(bkzgsa_gso_len(dvol, i, dim, delta=delta)) / log(2) for i in range(dim)]
blocksize_strategy = [(130,2,3),(130,2,1),(133,2,1),(134,2,2),(136,2,1),(138,2,1),(139,2,1),(141,2,2),(144,  4,  1),(147,  4,  1),(153, 4, 1)]
costs = strategy_est(l,blocksize_strategy,cost_model,progressive_sieve, worst_case)
    
   
print("==========")
n = 40
alpha = 0.035
#dim,dvol = gen_lwechal_instance(n, alpha)
dim = 188
dvol = 327.7388246557668
delta = compute_delta(2)
l = [log(bkzgsa_gso_len(dvol, i, dim, delta=delta)) / log(2) for i in range(dim)]
blocksize_strategy = []
costs = strategy_est(l,blocksize_strategy,cost_model,progressive_sieve, worst_case)


"""