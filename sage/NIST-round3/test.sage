import time


import math
import matplotlib.pyplot as plt
import os
from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt
from math import pi, exp, lgamma
from numpy.random import choice as np_random_choice
from numpy.linalg import slogdet, qr, cholesky
from numpy import array


'''BKZ simulator with pump and jump.'''

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



cd = [(rk[-i] - sum(rk[-i:])) / i for i in range(1, 46)]
cd += [(lgamma(beta_ / 2.0 + 1) * (1.0 / beta_) - log(sqrt(pi))) / log(2.0) for beta_ in range(46, 1000)]


def simulate_pnjBKZ_above_45(l,beta,jump=1,N=1,cd=cd):
    """
    BKZ simulator with pump and jump: block size larger than or equal to 45

    Args:
        l (list):list of log_2(||b_i^*||)
        beta (int): block size in BKZ prediction
        N (int): loop of BKZ-beta
        jump (int): gap in reduction

    Returns:
        list: list of log_2(||b_i^*||) after N times of BKZ-beta prediction
    """

    d = len(l)
    if N == 0 or beta == 0:
        return l
  
    l_ = [0 for _ in range(d)]
    J = jump
    k_ = min (len(rk), beta)

    for tours in range(N):
        flag = True
        sumf=0.
        sumk=0.
        ff = 0
    
        for k in range(d - beta):
            beta_ = min(beta, d - k)
            f = min(k - (k % J) + beta, d);
            sumf += sum(l[ff:f])
            if(k!=0):
                sumk += l_[k-1]
            logV = sumf - sumk
            l_k = logV / (beta_-(k % J)) + cd[(beta_-(k % J))-1]
            ff = f
            if flag:
                if l_k < l[k]:
                    l_[k] = l_k
                    flag = False
                else:
                    l_[k] = l[k]
            else:
                l_[k] = l_k

        sumf += sum(l[ff:d])
        for k in range(d - beta,d - 45):
            beta_ = d - k
            sumk += l_[k-1]
            logV = sumf - sumk
            l_k = logV / beta_ + cd[beta_-1]
            if flag == True:
                if l_k < l[k]:
                    l_[k] = l_k
                    flag = False
                else:
                    l_[k] = l[k]
            else:
                l_[k] = l_k
                
            #print(l_)
        
        #last 45 norms
        logV -= l_[d-46]
        for k in range(d - k_, d):
            l_[k] = logV/k_ + rk[k + 45 - d]
           
         #Set l[i] as l_[i] 
        l = deepcopy(l_)
        l_ = [0 for _ in range(d)]
    return l



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




chisquared_table = {i: None for i in range(1000)}


for i in range(5000):
    chisquared_table[i] = RealDistribution('chisquared', i)

dim_ = 1004
dvol = 3882.6780896

delta = compute_delta(2)

beta = 200

l = [log(bkzgsa_gso_len(dvol, i, dim_, delta=delta)) / log(2.) for i in range(dim_)]

l = simulate_pnjBKZ_above_45(l,beta)

proba =1.* chisquared_table[beta].cum_distribution_function(2**(2 * l[i]))

print(proba)