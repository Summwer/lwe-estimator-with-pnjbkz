from numpy.random import choice as np_random_choice
from numpy.linalg import slogdet, qr, cholesky
from numpy import array
from math import pi, exp, lgamma
import importlib
importlib.reload(sys)

print_style = {'SUCCESS': '\x1b[1;37;42m',
               'FAILURE': '\x1b[1;37;41m',
               'REJECT': '\x1b[3;31m',
               'ACCEPT': '\x1b[3;32m',
               'VALUE': '\x1b[1;33m',
               'DATA': '\x1b[3;34m',
               'WARNING': '\x1b[6;30;43m',
               'ACTION': '\x1b[1;37m',
               'HEADER': '\x1b[4;37m',
               'NORMAL': '\x1b[0m',
               }


def logging(message, style='NORMAL', newline=True):
    if style is not None:
        print(print_style[style], end=' ')
    print(message, end=' ')
    if newline:
        print('\x1b[0m')
    else:
        print('\x1b[0m', end=' ')
    sys.stdout.flush()

def GH_sv_factor_squared(k):
    return ((pi * k)**(1. / k) * k / (2. * pi * e))

def compute_delta(k):
    """Computes delta from the block size k. Interpolation from the following
    data table:
    Source : https://bitbucket.org/malb/lwe-estimator/
    src/9302d4204b4f4f8ceec521231c4ca62027596337/estima
    tor.py?at=master&fileviewer=file-view-default
    :k: integer
    estimator.py table:
    """

    small = {0: 1e20, 1: 1e20, 2: 1.021900, 3: 1.020807, 4: 1.019713, 5: 1.018620,
             6: 1.018128, 7: 1.017636, 8: 1.017144, 9: 1.016652, 10: 1.016160,
             11: 1.015898, 12: 1.015636, 13: 1.015374, 14: 1.015112, 15: 1.014850,
             16: 1.014720, 17: 1.014590, 18: 1.014460, 19: 1.014330, 20: 1.014200,
             21: 1.014044, 22: 1.013888, 23: 1.013732, 24: 1.013576, 25: 1.013420,
             26: 1.013383, 27: 1.013347, 28: 1.013310, 29: 1.013253, 30: 1.013197,
             31: 1.013140, 32: 1.013084, 33: 1.013027, 34: 1.012970, 35: 1.012914,
             36: 1.012857, 37: 1.012801, 38: 1.012744, 39: 1.012687, 40: 1.012631,
             41: 1.012574, 42: 1.012518, 43: 1.012461, 44: 1.012404, 45: 1.012348,
             46: 1.012291, 47: 1.012235, 48: 1.012178, 49: 1.012121, 50: 1.012065}

    if k != round(k):
        x = k - floor(k)
        d1 = compute_delta(floor(k))
        d2 = compute_delta(floor(k) + 1)
        return x * d2 + (1 - x) * d1

    k = int(k)
    if k < 50:
        return small[k]
    else:
        delta = GH_sv_factor_squared(k)**(1. / (2. * k - 2.))
        return delta.n()


def bkzgsa_gso_len(logvol, i, d, beta=None, delta=None):
    if delta is None:
        delta = compute_delta(beta)

    return RR(delta**(d - 1 - 2 * i) * exp(logvol / d))


rk = (
    0.789527997160000,
    0.780003183804613,
    0.750872218594458,
    0.706520454592593,
    0.696345241018901,
    0.660533841808400,
    0.626274718790505,
    0.581480717333169,
    0.553171463433503,
    0.520811087419712,
    0.487994338534253,
    0.459541470573431,
    0.414638319529319,
    0.392811729940846,
    0.339090376264829,
    0.306561491936042,
    0.276041187709516,
    0.236698863270441,
    0.196186341673080,
    0.161214212092249,
    0.110895134828114,
    0.0678261623920553,
    0.0272807162335610,
    -0.0234609979600137,
    -0.0320527224746912,
    -0.0940331032784437,
    -0.129109087817554,
    -0.176965384290173,
    -0.209405754915959,
    -0.265867993276493,
    -0.299031324494802,
    -0.349338597048432,
    -0.380428160303508,
    -0.427399405474537,
    -0.474944677694975,
    -0.530140672818150,
    -0.561625221138784,
    -0.612008793872032,
    -0.669011014635905,
    -0.713766731570930,
    -0.754041787011810,
    -0.808609696192079,
    -0.859933249032210,
    -0.884479963601658,
    -0.886666930030433,
)



cd = [rk[-i] - sum(rk[-i:]) / i for i in range(1, 46)]
cd += [(lgamma(beta_ / 2.0 + 1) * (1.0 / beta_) - log(sqrt(pi))) / log(2.0) for beta_ in range(46, 2000)]


chisquared_table = {i: None for i in range(2000)}


for i in range(5000):
    chisquared_table[i] = RealDistribution('chisquared', i)


rk = [0.789527997160000, 0.780003183804613, 0.750872218594458, 0.706520454592593, 0.696345241018901, 0.660533841808400, 0.626274718790505, 0.581480717333169, 0.553171463433503, 0.520811087419712, 0.487994338534253, 0.459541470573431, 0.414638319529319, 0.392811729940846, 0.339090376264829, 0.306561491936042, 0.276041187709516, 0.236698863270441, 0.196186341673080, 0.161214212092249, 0.110895134828114, 0.0678261623920553, 0.0272807162335610, -
      0.0234609979600137, -0.0320527224746912, -0.0940331032784437, -0.129109087817554, -0.176965384290173, -0.209405754915959, -0.265867993276493, -0.299031324494802, -0.349338597048432, -0.380428160303508, -0.427399405474537, -0.474944677694975, -0.530140672818150, -0.561625221138784, -0.612008793872032, -0.669011014635905, -0.713766731570930, -0.754041787011810, -0.808609696192079, -0.859933249032210, -0.884479963601658, -0.886666930030433]
simBKZ_c = [None] + [rk[-i] - sum(rk[-i:]) / i for i in range(1, 46)]

pruning_proba = .5
simBKZ_c += [RR(log(GH_sv_factor_squared(d)) / 2. -
                log(pruning_proba) / d) / log(2.) for d in range(46, 1000)]
#simBKZ_c += [(lgamma(beta_ / 2.0 + 1) * (1.0 / beta_) - log(sqrt(pi))) / log(2.0) for #beta_ in range(46, 1000)]

def simBKZ(l, beta, tours=1, c=simBKZ_c):

    n = len(l)
    l2 = copy(vector(RR, l))

    for k in range(n - 1):
        d = min(beta, n - k)
        f = k + d
        logV = sum(l2[k:f])
        lma = logV / d + c[d]
    
        if lma >= l2[k]:
            continue

        diff = l2[k] - lma
        l2[k] -= diff
        for a in range(k + 1, f):
            l2[a] += diff / (f - k - 1)

    return l2


def gaussian_heuristic(log_rr):
    """
    Return squared norm of shortest vector as predicted by the Gaussian heuristic.

    :param log_rr: vector of ln(squared) Gram-Schmidt norms

    """
    n = len(list(log_rr))
    log_vol = sum(log_rr)
    log_gh =  1./n * (log_vol - 2 * ball_log_vol(n))
    return exp(log_gh)



def ball_log_vol(n):
    """
    Return volume of `n`-dimensional unit ball

    :param n: dimension

    """
    return (n/2.) * log(pi) - lgamma(n/2. + 1)



