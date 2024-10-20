#!/usr/bin/env python
# -*- coding: utf-8 -*-


from math import e, lgamma, log, pi
from fpylll import *
load("../framework/load_lwechal.sage")
load("../framework/simulator/pnjbkz_simulator.sage")


def dim4free_wrapper(dim4free_fun, blocksize):
    """
    Deals with correct dim4free choices for edge cases when non default
    function is chosen.

    :param dim4free_fun: the function for choosing the amount of dim4free
    :param blocksize: the BKZ blocksize

    """
    if blocksize < 40:
        return 0
    dim4free = dim4free_fun(blocksize)
    return int(min((blocksize - 40)/2, dim4free))


def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)


def delta_0f(k):
    """
    Auxiliary function giving root Hermite factors. Small values
    experimentally determined, otherwise from [Chen13]

    :param k: BKZ blocksize for which the root Hermite factor is required

    """
    small = (( 2, 1.02190),  # noqa
             ( 5, 1.01862),  # noqa
             (10, 1.01616),
             (15, 1.01485),
             (20, 1.01420),
             (25, 1.01342),
             (28, 1.01331),
             (40, 1.01295))

    k = float(k)
    if k <= 2:
        return (1.0219)
    elif k < 40:
        for i in range(1, len(small)):
            if small[i][0] > k:
                return (small[i-1][1])
    elif k == 40:
        return (small[-1][1])
    else:
        return (k/(2*pi*e) * (pi*k)**(1./k))**(1/(2*(k-1.)))


def log_gh_svp(d, delta_bkz, svp_dim, n, q):
    """
    Calculates the log of the Gaussian heuristic of the context in which
    SVP will be ran to try and discover the projected embedded error.

    The volume component of the Gaussian heuristic (in particular the lengths
    of the appropriate Gram--Schmidt vectors) is estimated using the GSA
    [Schnorr03] with the multiplicative factor = delta_bkz ** -2.

    NB, here we use the exact volume of an n dimensional sphere to calculate
    the ``ball_part`` rather than the usual approximation in the Gaussian
    heuristic.

    :param d: the dimension of the embedding lattice = n + m + 1
    :param delta_bkz: the root Hermite factor given by the BKZ reduction
    :param svp_dim: the dimension of the SVP call in context [d-svp_dim:d]
    :param n: the dimension of the LWE secret
    :param q: the modulus of the LWE instance

    """
    d = float(d)
    svp_dim = float(svp_dim)
    ball_part = ((1./svp_dim)*lgamma((svp_dim/2.)+1))-(.5*log(pi))
    vol_part = ((1./d)*(d-n-1)*log(q))+((svp_dim-d)*log(delta_bkz))
    return ball_part + vol_part


def gsa_params(n, alpha, q=None, samples=None, d=None, decouple=False):
    """
    Finds winning parameters (a BKZ reduction dimension and a final SVP call
    dimension) for a given Darmstadt LWE instance (n, alpha).

    :param n: the dimension of the LWE secret
    :param alpha: the noise rate of the LWE instance
    :param q: the modulus of the LWE instance. ``None`` means determine by
        reloading the challenge
    :param samples: maximum number of LWE samples to use for the embedding
        lattice. ``None`` means ``5*n``
    :param d: find best parameters for a dimension ``d`` embedding lattice
    :param decouple: if True the BKZ dimension and SVP dimension may differ

    """
    if q is None or samples is None:
        A, _, q = load_lwe_challenge(n, alpha)
        samples = A.nrows()

    stddev = alpha*q
    
    params = decoupler(decouple, n, samples, q, stddev, d)
    min_cost_param = find_min_complexity(params)
    if min_cost_param is not None:
        return min_cost_param


def decoupler(decouple, n, samples, q, stddev, d):
    """
    Creates valid (bkz_dim, svp_dim, d) triples, as determined by
    ``primal_parameters`` and determines which succeed in the recovery of the
    embedded error.

    :param decouple: if True the BKZ dimension and SVP dimension may differ
    :param n: the dimension of the LWE secret
    :param samples: maximum number of LWE samples to use for the embedding
        lattice. ``None`` means ``5*n``
    :param q: the modulus of the LWE instance
    :param stddev: the standard deviation of the distribution from which the
        error vector components were uniformly and indepedently drawn
    :param d: find best parameters for dimension ``d`` embedding lattice

    """
    params = []

    ms = range(n, min(5*n+1, samples+1))
    
    for m in ms:
        beta_bound = min(m+1, 400+default_dim4free_fun(400))
        svp_bound = min(m+1, 400)
        for bkz_block_size in range(40, beta_bound):
            delta_0 = delta_0f(bkz_block_size)
            if decouple:
                svp_dims = range(40, svp_bound)
            else:
                svp_dims = [min(bkz_block_size, svp_bound)]

            for svp_dim in svp_dims:
                d = float(m + 1)
                rhs = log_gh_svp(d, delta_0, svp_dim, n, q)
                if rhs - log(stddev) - log(svp_dim)/2. >= 0:
                    params.append([bkz_block_size, svp_dim, m+1])

    return params


def find_min_complexity(params):
    """
    For each valid and solving triple (bkz_dim, svp_dim, d) determines an
    approximate (!) cost and minimises.

    :param params: a list of all solving (bkz_dim, svp_dim, d) triples

    """
    min_cost = None
    min_cost_param = None
    expo = .349

    for param in params:

        bkz_block_size = param[0] - default_dim4free_fun(param[0])
        svp_dim = param[1] - default_dim4free_fun(param[1])
        d = param[2]

        bkz_cost = 2 * d * (2 ** (expo * bkz_block_size))
        finisher_svp_cost = 2 ** ((expo * svp_dim))
        new_cost = bkz_cost + finisher_svp_cost

        if min_cost is None or new_cost < min_cost:
            min_cost = new_cost
            min_cost_param = param

    return min_cost_param


def sim_params(n, alpha):
    A, c, q = load_lwe_challenge(n, alpha)
    stddev = alpha*q
    winning_params = []
    for m in range(60, min(2*n+1, A.nrows+1)):
        B = primal_lattice_basis(A, c, q, m=m)
        M = GSO.Mat(B)
        M.update_gso()
        beta_bound = min(m+1, 500+default_dim4free_fun(500)+1)
        svp_bound = min(m+1, 500)
        rs = [M.get_r(i, i) for i in range(M.B.nrows)]
        for beta in range(40, beta_bound):
            rs, _ = simulate(rs, fplll_bkz.EasyParam(beta, max_loops=1))
            for svp_dim in range(40, svp_bound):
                gh = gaussian_heuristic(rs[M.B.nrows-svp_dim:])
                if svp_dim*(stddev**2) < gh:
                    winning_params.append([beta, svp_dim, m+1])
                    break
    min_param = find_min_complexity(winning_params)
    return min_param


def primal_lattice_basis(A, c, q, m=None):
    """
    Construct primal lattice basis for LWE challenge
    ``(A,c)`` defined modulo ``q``.

    :param A: LWE matrix
    :param c: LWE vector
    :param q: integer modulus
    :param m: number of samples to use (``None`` means all)

    """
    if m is None:
        m = A.nrows()
    elif m > A.nrows():
        raise ValueError("Only m=%d samples available." % A.nrows())
    n = A.ncols()

    B = IntegerMatrix(m+n+1, m+1)
    for i in range(m):
        for j in range(n):
            B[j, i] = A[i, j]
        B[i+n, i] = q
        B[-1, i] = c[i]
    B[-1, -1] = 1

    B = LLL.reduction(B)
    assert(B[:n] == IntegerMatrix(n, m+1))
    B = B[n:]

    return B




def gen_lwechal_instance(n=40, alpha=0.005, default_g6k = False):
    A, c, q = load_lwe_challenge(n=n, alpha=alpha)
    
    print("-------------------------")
    print("Primal attack, TU LWE challenge n=%d, alpha=%.4f, q = %d. " % (n, alpha, q))

    try:
        min_cost_param = gsa_params(n=A.ncols(), alpha=alpha, q=q, decouple=True)
        (b, s, m) = min_cost_param
    except TypeError:
        raise TypeError("No winning parameters.")
   
    print("Chose %d samples. Predict solution at bkz-%d + svp-%d." % (m, b, s))
    
    d = m + 1

    B = primal_lattice_basis(A, c, q, m=m)

    sigma = alpha * q
    
    M = GSO.Mat(B)
    M.update_gso()
    rr = [M.get_r(i,i) for i in range(d)]
    if(default_g6k):
        log_rr = [log(rr[i],2)/2. for i in range(d)]
    else:
        log_rr = [log(rr[i],2)/2. - log(sigma,2) for i in range(d)]
        sigma = 0.
    print("Initial slope: ", get_current_slope(log_rr,0,d))
    
    dvol = sum(log_rr) * log(2)  #ln(vol)

    dim = m + 1
    print("dim = %3d, dvol = %3.7f" %(dim, dvol))
    print()

    # return (dim, dvol)
    return (log_rr,dim,dvol,b,sigma)
