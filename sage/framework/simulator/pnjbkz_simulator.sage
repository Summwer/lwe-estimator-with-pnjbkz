import math
import matplotlib.pyplot as plt
import os
from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt

load("../framework/utils.sage")
#load("../framework/simulator/CN11.sage")

'''BKZ simulator with pump and jump.'''


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
    return int(min((blocksize - 40)/2., dim4free))


def default_dim4free_fun(blocksize):
    """
    Return expected number of dimensions for free, from exact-SVP experiments.

    :param blocksize: the BKZ blocksize

    """
    return int(11.5 + 0.075*blocksize)


def theo_dim4free_fun1(blocksize):
    """
    Theoretical Dimension-for-free function 1 without e in [Duc18]
    """

    return max(0, ceil(blocksize*log(4/3.)/log(blocksize/2./pi)) )


def theo_dim4free_fun2(blocksize):
    """
    Theoretical Dimension-for-free function 2 with e in [Duc18]
    """

    return max(0, ceil(blocksize*log(4/3.)/log(blocksize/2./pi/e)) )


#d4f according to specific lattice: ||pi_f(s)||<= sqrt(4/3) GH(L_f)
#def accurate_d4f(l_, dsvp):
#    for f in range(0, d - dsvp):
#        gh = gaussian_heuristic(l_[d-dsvp - f:])
#        psvp = 1. * chisquared_table[dsvp].cum_distribution_function(4/3. * gh)
#        if()





def get_beta_from_sieve_dim(sieve_dim,d,dim4free_fun):
    for beta in range(sieve_dim,d):
        f = dim4free_wrapper(dim4free_fun,beta)
        # print(beta,f,beta-f,sieve_dim)
        if beta - f >= sieve_dim:
            return beta

    

def simulate_pnjBKZ(l0, beta, jump, loop):
    """
    BKZ simulator with pump and jump

    Args:
        l0 (list): list of log_2(||b_i^*||)
        beta (int): block size in BKZ prediction
        jump (int): gap in reduction
        loop (int): loop of BKZ-beta

    Returns:
        list: list of log_2(||b_i^*||) after N times of BKZ-beta prediction
    """
    
    l = deepcopy(l0)

    #extra_dim4free = 12
    #f = dim4free_wrapper(default_dim4free_fun, beta)
    #if jump <=2:
    #    beta_ = beta
    #elif jump>=3 and jump <=4:
    #    beta_ = get_beta_from_sieve_dim(beta-f,d,theo_dim4free_fun2)
    #elif jump>=5:
    #    beta_ = get_beta_from_sieve_dim(beta-f,d,theo_dim4free_fun1)


    #if beta <=50 and jump ==1:
    #    return simBKZ(l, beta, tours=loop) #bkz with enumeration oracle

    #bkz with sieve oracle
    if beta <= 45:
        return simulate_pnjBKZ_below_45(l, beta, loop) 
    else:
        return simulate_pnjBKZ_above_45(l, beta, jump, loop)
        
        

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


#Implement the simulator in beta below 45 as CN11
def simulate_pnjBKZ_below_45(l,beta,N=1,cd=cd):
    """
    BKZ simulator with pump and jump: block size less than 45

    Args:
        l0 (list):list of log_2(||b_i^*||)
        beta (int): block size in BKZ prediction
        N (int): loop of BKZ-beta

    Returns:
        list: list of log_2(||b_i^*||) after N times of BKZ-beta prediction
    """

    d = len(l)
    
    if N == 0 or beta == 0:
        return l

    l_ = deepcopy(l)
    k_ = min (len(rk), beta)
    
    for j in range(N):
        flag = True
        sumf=0.
        sumk=0.
        ff = 0
        for k in range(d-beta):
            beta_ = min(beta, d - k)
            f = k+beta_
            sumf += sum(l[ff:f])
            if(k!=0):
                sumk += l_[k-1]
            logV = sumf - sumk
            l_k = logV / beta_ + cd[beta_-1]
            ff = f
            if flag == True:
                if l_k < l[k]:
                    l_[k] = l_k
                    flag = False
            else:
                l_[k] = l_k
        # early termination
        if flag or l_ == l:
            break
        else:
            #last beta elements
            sumf += l[d-1]
            sumk += l_[d-beta-1]
            logV = sumf - sumk
            tmp = sum(rk[-beta:]) / beta
            rk1 = [r_ - tmp for r_ in rk[-beta :]]
            for k in range(d - beta, d):
                l_[k] = logV/beta + rk1[k + beta - d]

            #Set l[i] as l_[i]         
            l = deepcopy(l_)
            l_ = [0 for _ in range(d)]
    return l

    




def show_gs_slope_figure(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error,beta,N):
    plt.figure(figsize=(15, 10), dpi=100)
    # plt.ylim(4,17) #set range of y_ticks
    # plt.xlim(-5,210)
    t = 0
    plt.scatter([_+1 for _ in range(dimension)],log_gs_length,marker="*")#,c = color)
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="*")#,c = color)
    plt.title("n=%d, alpha = %s, dimension = %d, Blocksize = %d, Current Tours= %d,square error = %f" %(n,alpha_,dimension,beta,N,square_error))
    
    try:
        os.mkdir(dir+"gs-lengths-gh simulator/")
    except FileExistsError:
        pass
    #plt.savefig(dir+"n=%d,alpha=%s,d=%d.png" %(n,alpha_,dimension))
    plt.show()

def compute_square_error(list1,list2,flag = 1):
    square_error = 0
    l = min(len(list1),len(list2))
    L = max(len(list1),len(list2))
    for i in range(L):

        if i == 0 and list2[1]-list2[0] > 1 and flag == 1:
            #remove the error point.
            continue
        elif i < l :
            square_error += (list1[i]-list2[i])**2
        elif i>=l:
            try:
                square_error += list1[i]**2
            except IndexError:
                square_error += list2[i]**2
            
    return square_error/sum(_**2 for _ in list1)
             


def get_current_slope(r, start_row=0, stop_row=-1):
    """
    A Python re-implementation of ``MatGSO.get_current_slope``.

        >>> from fpylll import IntegerMatrix, GSO, LLL, FPLLL
        >>> FPLLL.set_random_seed(1337)
        >>> A = IntegerMatrix.random(100, "qary", bits=30, k=50)
        >>> _ = LLL.reduction(A)
        >>> M = GSO.Mat(A); _ = M.update_gso()
        >>> from fpylll.tools.quality import get_current_slope
        >>> M.get_current_slope(0, 100)  # doctest: +ELLIPSIS
        -0.085500625...
        >>> get_current_slope(M.r(), 0, 100) # doctest: +ELLIPSIS
        -0.085500625...

    """
    x = [2.*r[i]*log(2.) for i in range(start_row, stop_row)]
    n = stop_row - start_row
    i_mean = (n - 1) * 0.5 + start_row
    x_mean = sum(x)/n
    v1, v2 = 0.0, 0.0
    for i in range(stop_row - start_row):
        v1 += (i - i_mean) * (x[i] - x_mean)
        v2 += (i - i_mean) * (i - i_mean)
    return v1 / v2



