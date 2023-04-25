load("../framework/attack.sage")
load("../framework/lwe_estimation.sage")
load("../framework/load_lwechal.sage")

######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3, parallel_ = False: EnumBS
             3, parallel_=True: EnumBS in parallel
             4: bssa
             5: bkz-only with fixed blocksize
             6: default g6k
             
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################

def gen_lwechal_instance(n=40, alpha=0.005):
    A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    
    print("-------------------------")
    print("Primal attack, TU LWE challenge n=%d, alpha=%.4f" % (n, alpha))

    try:
        min_cost_param = gsa_params(n=A.ncols(), alpha=alpha, q=q, decouple=True)
        (b, s, m) = min_cost_param
    except TypeError:
        raise TypeError("No winning parameters.")
   
    print("Chose %d samples. Predict solution at bkz-%d + svp-%d" % (m, b, s))
    print()

    d = m + 1

    B = primal_lattice_basis(A, c, q, m=m)

    sigma = alpha * q
    G, M = B.gram_schmidt()
    G = matrix(RDF, G)
    rr = [sum(G[i,j]**2 for j in range(d)) for i in range(d)]
    log_rr = [log(rr[i])/2 - log(sigma) for i in range(d)]
    dvol = sum(log_rr)

    dim = m + 1
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))


    
    return (dim, dvol)




def low_dim_lwe_challenge_est(method,  cost_model, worst_case, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False):
    lwechals = [(40,0.025), (45,0.020), (50,0.015), (55,0.010), (60,0.010), (70,0.005), (75,0.005)]
    for (n,alpha) in lwechals:
        dim,dvol = gen_lwechal_instance(n, alpha)

        estimate_attack( silent=False, method = method, progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)



  


def unsolved_lwe_challenge_est(method,  cost_model, worst_case, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False, progressive_sieve = True, print_l = False):
    lwechals = [(40,0.045), (45,0.035), (50,0.030), (55,0.025), (60,0.020), (65,0.015), (75,0.010), (95,0.005)]
    for (n,alpha) in lwechals:
        dim,dvol = gen_lwechal_instance(n, alpha)

        estimate_attack( silent=False, method = method, progressive_sieve = progressive_sieve, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim,J=J,gap =gap,J_gap = J_gap,cost_model=cost_model,gen_GSA_gso=gen_GSA_gso,print_l = print_l, worst_case = worst_case)

