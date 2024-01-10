import time
from math import factorial as fac
load("../framework/attack.sage")

######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3: bkz-only with fixed blocksize
             4: default g6k
             "fixed-blocksize bkz-only": call fixed-blocksize bkz-only model to estimate security.
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################

#binomial distribution
def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2 * k, x + k) / 2.**(2 * k)


def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k + 1):
        D[i] = centered_binomial_pdf(k, i)
    return D

def build_Gaussian_law(sigma, t):
    D = {}
    for i in range(0, t + 1):
        D[i] = exp(-i ** 2 / (2 * sigma ** 2))
        D[-i] = D[i]
    normalization = sum([D[i] for i in D])
    for i in D:
        D[i] = D[i] / normalization
    assert abs(sum([D[i] for i in range(-t, t + 1)]) - 1.) <= 10 ** -10
    return D


def build_uniform_law(p):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(p):
        D[i - p / 2] = 1. / p
    return D




#ROUNDING_FACTOR = 2**64

#def round_to_rational(x):
#    A = ZZ(round(x * ROUNDING_FACTOR))
#    return QQ(A) / QQ(ROUNDING_FACTOR)

def average_variance(D):
    mu = 0.
    s = 0.
    
    for (v, p) in D.items():
        mu += v * p
        s += v * v * p

    s -= mu * mu
    return mu, s
    #return round_to_rational(mu), round_to_rational(s)


def compute_beta(d, logvol):
    """
    Computes the beta value for given dimension and volumes
    It is assumed that the instance has been normalized and sphericized, 
    i.e. that the covariance matrices of the secret is the identity
    :d: integer
    :vol: float
    """
    bbeta = None
    pprev_margin = None

    # Keep increasing beta to be sure to catch the second intersection
    for beta in range(2, d):
        lhs = RR(sqrt(beta))
        rhs = bkzgsa_gso_len(logvol, d - beta, d, beta=beta)
        
        if lhs < rhs and bbeta is None:
            margin = rhs / lhs
            prev_margin = pprev_margin
            bbeta = beta
            #print(rhs,lhs)

        if lhs > rhs:
            bbeta = None
        pprev_margin = rhs / lhs

    if bbeta is None:
        return None

    ddelta = compute_delta(bbeta) * margin**(1. / d)
    if prev_margin is not None:
        beta_low = log(margin) / (log(margin) - log(prev_margin))
    else:
        beta_low = 0
    assert(beta_low >= 0)
    assert(beta_low <= 1)
    return bbeta - beta_low


def initialize_from_LWE_instance(n, q, m, D_e, D_s, verbosity = True):
    """
    constructor that builds a DBDD instance from a LWE instance
    :n: (integer) size of the secret s
    :q: (integer) modulus
    :m: (integer) size of the error e
    :D_e: distribution of the error e (dictionnary form)
    :D_s: distribution of the secret s (dictionnary form)
    """
    if verbosity:
        logging("     Build LWE  instance  ", style="HEADER")
        logging("n=%3d \t m=%3d \t q=%d" % (n, m, q), style="VALUE")
    # define the mean and sigma of the instance
    mu_e, s_e = average_variance(D_e)
    mu_s, s_s = average_variance(D_s)
    #Find the best m.
    min_beta = -1
    for test_m in range(m,0,-1):
        dvol = test_m * log(q) - (test_m*log(abs(s_e))+n*log(abs(s_s))) / 2.
        d = test_m+n+1
        est_beta = compute_beta(d, dvol)
        if(dvol <= 0 or est_beta is None):
            break
        if(min_beta == - 1 or min_beta > est_beta ):
            min_beta = est_beta
            best_m = test_m
        
    m = best_m
    dvol = m * log(q) - (m*log(abs(s_e))+n*log(abs(s_s))) / 2.
    dim = m + n + 1
    print("dim = %d, m = %d, dvol = %3.11f, β = %3.4f" %(dim, m, dvol, min_beta))


    core_SVP_est(dim, dvol)
    return dim, dvol


def core_SVP_est(d, dvol):
    for beta in range(50,d):
        delta0 = compute_delta(beta)
        #[ADPS16] uses 2beta-d-1, it is 2beta-d actually.
        if(1/2.*log(beta) <= (2*beta-d-1)*log(delta0) + 1./d * dvol ):
            f = dim4free_wrapper(theo_dim4free_fun1,beta)
            print("core-SVP estimate, β = ", beta, "beta-f = ", beta-f, "cost =", 0.292*(beta-f))
            return



def dilithium_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False,  print_l = False, ldc_param = "AGPS20", cal_ee = "chi", worst_case = False, goal_min_cost = "gate_min", cumG = False):
    
    # Dilithium-II round-3 parameters
    print("============= Dilithium-II")
    n = 4*256
    m = 4*256
    q = 8380417
    eta = 2
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)}
    D_e = D_s
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    
    print("============= Dilithium-III")

    # Dilithium-III round-3 parameters
    dim_ = 2654
    dvol = 19371.0238433
    n = 5*256
    m = 6*256
    q = 8380417
    eta = 4
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)}
    D_e = D_s
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)
    

    print("============= Dilithium-V")

    # Dilithium-V round-3 parameters
    dim_ = 3540
    dvol = 26623.1162463
    n = 7*256
    m = 8*256
    q = 8380417
    eta = 2
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)}
    D_e = D_s
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l,  ldc_param = ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

def kyber_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False,  print_l = False, ldc_param = "AGPS20", cal_ee = "chi", worst_case = False,goal_min_cost = "gate_min", cumG = False):

    # Kyber-512 round-3 parameters
    print("============= Kyber-512 with σs = 1.5 and σe = 1.5")
    #eta1 = eta2 = 3
    #dim_ = 1025
    #dvol = 3944.9406103
    #eta1 = 3, eta2 = 2
    #eta1 =3, eta2 =2
    n = 512
    m = 512
    q = 3329
    D_s = build_centered_binomial_law(3)
    D_e = D_s
    #D_e = build_centered_binomial_law(2)
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    
    print("============= Kyber-512 with σs = 1.5 and σe = 1")
    #eta1 = eta2 = 3
    #dim_ = 1025
    #dvol = 3944.9406103
    #eta1 = 3, eta2 = 2
    #eta1 =3, eta2 =2
    n = 512
    m = 512
    q = 3329
    D_s = build_centered_binomial_law(3)
    D_e = build_centered_binomial_law(2)
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)




    print("============= Kyber-768")

    # Kyber-768 round-3 parameters
    #dim_ = 1467
    #dvol = 5661.0782118
    n = 768
    m = 768
    q = 3329
    D_s = build_centered_binomial_law(2)
    D_e = D_s
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)


    print("============= Kyber-1024")

    # Kyber-1024 round-3 parameters
    #dim_ = 1918
    #dvol = 7242.6115232 
    n = 1024
    m = 1024
    q = 3329
    D_s = build_centered_binomial_law(2)
    D_e = D_s
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l, ldc_param = ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)




def hufu_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False,  print_l = False, ldc_param = "AGPS20", cal_ee = "chi", worst_case = False,goal_min_cost = "gate_min", cumG = False):

    # HuFu-I parameters
    print("============= HuFu-I")
    n = 848
    m = 736
    q = int(pow(2,4))
    p = int(pow(2,12))
    Q = p * q
    D_s = build_centered_binomial_law(1) #binomial distribution on B1
    D_e = build_centered_binomial_law(1) #binomial distribution on B1
    dim_, dvol = initialize_from_LWE_instance(n, Q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)




    # HuFu-III parameters
    print("============= HuFu-III")
    n = 1232
    m = 1024
    q = int(pow(2,4))
    p = int(pow(2,13))
    Q = p * q
    D_s = build_centered_binomial_law(1) #binomial distribution on B1
    D_e = build_centered_binomial_law(1) #binomial distribution on B1
    dim_, dvol = initialize_from_LWE_instance(n, Q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)



    # HuFu-V parameters
    print("============= HuFu-V")
    n = 1552
    m = 1312
    q = int(pow(2,4))
    p = int(pow(2,13))
    Q = p * q
    D_s = build_centered_binomial_law(1) #binomial distribution on B1
    D_e = build_centered_binomial_law(1) #binomial distribution on B1
    dim_, dvol = initialize_from_LWE_instance(n, Q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)



def eaglesign_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False,  print_l = False, ldc_param = "AGPS20", cal_ee = "chi", worst_case = False,goal_min_cost = "gate_min", cumG = False):
    
    q = 12289

    # EagleSign2 parameters
    p = 512
    k = 2
    l = 1
    n = p*l
    m = p*k
    print("============= EagleSign2LPk: Longterm secret recovery")
    eta = 6 #etag=etad
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over [-etag,etag]
    D_e = D_s #uniformly random over [-etad,etad]
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    print("============= EagleSign2EPk: Ephemeral secret recovery")
    t = 140
    eta = 64 #etay2
    D_s = {1: t/(2.*p), -1: t/(2.*p), 0: (2.*p -t)/p} #SparseTernary: t/2 entries of 1,  t/2 entries of -1, rest is 0.
    D_e = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over [-eta,eta], eta = etay2
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)
    
    
    # EagleSign3 parameters
    p = 1024
    k = 1
    l = 1
    n = p*l
    m = p*k
    #print("============= EagleSign3LPk: Longterm secret recovery")
    eta = 1
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over {-1,1}
    D_e = D_s #uniformly random over {-1,1}
    #dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    #estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)
    
    print("============= EagleSign3EPk: Ephemeral secret recovery")
    t = 140
    eta = 64 #etay2
    D_s = {1: t/(2.*p), -1: t/(2.*p), 0: (2.*p -t)/p} #SparseTernary: t/2 entries of 1,  t/2 entries of -1, rest is 0.
    D_e = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over [-eta,eta], eta = etay2
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)


    # EagleSign5 parameters
    p = 1024
    q = 12289
    k = 1
    l = 2
    n = p*l
    m = p*k
    print("============= EagleSign5LPk: Longterm secret recovery")
    eta = 1
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over {-1,1}
    D_e = D_s #uniformly random over {-1,1}
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    print("============= EagleSign5EPk: Ephemeral secret recovery")
    t = 86
    eta = 32 #etay2
    D_s = {1: t/(2.*p), -1: t/(2.*p), 0: (2.*p -t)/p} #SparseTernary: t/2 entries of 1,  t/2 entries of -1, rest is 0.
    D_e = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over [-eta,eta], eta = etay2
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    




def haetae_est(method, J=1, gap=1, J_gap=1, l = None, gen_GSA_gso = True, parallel_ = False,  print_l = False, ldc_param = "AGPS20", cal_ee = "chi", worst_case = False,goal_min_cost = "gate_min", cumG = False):

    # HAETAE-II parameters
    print("============= HAETAE-II")
    q = 64513
    p = 256
    k = 2
    l = 4
    n = p*(l-1)
    m = p*k
    eta = 1
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over {-1,1}
    D_e = D_s #uniformly random over {-1,1}
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    # HAETAE-III parameters
    print("============= HAETAE-III")
    q = 64513
    p = 256
    k = 3
    l = 6
    n = p*(l-1)
    m = p*k
    eta = 1
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over {-1,1}
    D_e = D_s #uniformly random over {-1,1}
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

    # HAETAE-V parameters
    print("============= HAETAE-V")
    q = 64513
    p = 256
    k = 4
    l = 7
    n = p*(l-1)
    m = p*k
    eta = 1
    D_s = {x : 1./(2*eta+1) for x in range(-eta, eta+1)} #uniformly random over {-1,1}
    D_e = D_s #uniformly random over {-1,1}
    dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
    estimate_attack( silent=False, method = method, parallel_ = parallel_, l = l, dvol = dvol, dim_ = dim_,J=J,gap =gap,J_gap = J_gap,gen_GSA_gso=gen_GSA_gso,print_l = print_l  ,ldc_param =  ldc_param, cal_ee = cal_ee, worst_case = worst_case, goal_min_cost = goal_min_cost, cumG = cumG)

