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

def cost_convexity_test(l,strategy,cost_model,progressive_sieve, worst_case):
    costs = []
    d = len(l)
    cum_pr = 0.
    rem_pr = 1.
    G1cum, B1cum = 0., 0.
    for bs in strategy:
        (beta, jump, N) = bs
        if(beta>d):
            break
        for tours in range(N):
            l = simulate_pnjBKZ(l, beta, jump, 1)
            i = d - beta
            if(cost_model == 2 )
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

            #print(beta,jump,tours+1,G1cum,G_sieve)
            costs.append(log2(2**G1cum+2**G_sieve))
            #print(beta)
            if(cum_pr > 0.999):
                return costs
    



def random_gen_blocksize_strategy(betamin_, J, loop, dim):
    blocksize_strategy = []
    while(len(blocksize_strategy)< 5):
        blocksize_strategy = []
        betamin = betamin_
        while(betamin <= dim):
            betamin = randint(betamin, min(betamin+20,dim))
            j = max(1,min(dim4free_wrapper(default_dim4free_fun, betamin), randint(1,J)))
            t = randint(1,loop)
            blocksize_strategy.append((betamin,j,t))
            betamin = betamin+1

    #print(blocksize_strategy)
    return blocksize_strategy



worst_case = True
cost_model = 1
progressive_sieve = True
path = "cost_convexity_test/randomly/"

try:
    os.mkdir( path)
except FileExistsError:
    pass



lwechals = [#(40,0.005), (40,0.010), (40,0.015), (40,0.020), (40,0.025), (40,0.030), (40,0.035),(40,0.040),
            #(45,0.005), (45,0.010), (45,0.015), (45,0.020), (45,0.025), (45,0.030), (45,0.035), 
            #(50,0.005), (50,0.010), (50,0.015), (50,0.020), (50,0.025), (50,0.030),
            #(55,0.005), (55,0.010), (55,0.015), (55,0.020), (55,0.025), 
            #(60,0.005), (60,0.010), (60,0.015), (60,0.020), 
            #(65,0.005), (65,0.010), (65,0.015), 
            #(70,0.005), (70,0.010), (70,0.015), 
            (75,0.005), (75,0.010),
            #(80,0.005), (80,0.010), 
            #(85,0.005), (85,0.010), 
            #(90,0.005), (90,0.010), 
            #(95,0.005)
            ]

betamin = 50
J = 8
loop = 5

mincost = (149.81,0)

while(mincost[0] >= 149.81):
    plt.cla()
    plt.close("all")
    #dim,dvol = gen_lwechal_instance(n, alpha)
    #dilithium1
    print("============= Dilithium-I")
    dim = 2049
    dvol = 15614.2193172
    delta = compute_delta(2)
    l = [log(bkzgsa_gso_len(dvol, i, dim, delta=delta)) / log(2) for i in range(dim)]
    blocksize_strategy = random_gen_blocksize_strategy(betamin, J, loop, dim)
    
    xs = []
    for Tuple in blocksize_strategy:
        for t in range(Tuple[2]):
            xs.append(str((Tuple[0],Tuple[1],t+1)))

    costs = cost_convexity_test(l,blocksize_strategy,cost_model,progressive_sieve, worst_case)
    mincost = (costs[0],0)
    convexity = []
    for i in range(1, len(costs)-1):
        if(costs[i-1] > costs[i] and costs[i] < costs[i+1]):
            convexity.append(((costs[i-1],xs[i-1]),(costs[i],xs[i]),(costs[i+1],xs[i+1])))
        if(mincost[0] > costs[i]):
            mincost = (costs[i],i)
    
    print("convexity = ", convexity, ", number of convexity = ", len(convexity)," (should <= 1)") 
    print("mincost = %f, index = %d" %(mincost[0],mincost[1]))
    print()
    #alpha = int(round(alpha * 1000))
    #
    plt.figure(figsize=(20,5))
    
    #print(xs)
    #print(len(xs))
    #print(len(costs))
    #xs[:len(costs)]
    print("[", ",".join(xs[:len(costs)]),"]")

    fig = plt.gcf()
    plt.scatter(xs[:len(costs)],costs, marker=".") #, label="dilithium1"
    #plt.legend()
    plt.xticks(rotation=45, fontsize=14)
    plt.subplots_adjust(bottom=0.3)
    Title = "dilithium1"
    #Title += ": ["+ ",".join(xs[:len(costs)])+"]"
    plt.title(Title)
    plt.savefig(path+"dilithium1.png")
    plt.close(fig)

    #assert(len(convexity)<=1)
   

