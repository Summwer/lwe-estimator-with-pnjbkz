import math
import matplotlib.pyplot as plt
import os
import shutil
from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma
from math import sqrt

'''v2: Add the jump value into the simulator'''

load("../framework/utils.sage")



#beta>45
def pump_simulator(log_rr0,d,llb,cd=cd):
    
    # print(db_size_base)
    N=1
    r = d
    beta = r - llb
    #dim4f = dim4free_wrapper(theo_dim4free_fun1, beta)

    db_size_base = 0 #log(sqrt(4/3. * (beta / (beta-dim4f))))
    # print(db_size_base)
    if beta == 0:
        return log_rr0
    l = [_ for _ in log_rr0]
    #print(len(l))
    l_ = []
    
    
    for tours in range(N):
        flag = True
        for k in range(d - beta):
            l_.append(l[k])

        for k in range(d - beta ,d - 45):
            beta_ = d - k
            f = d
            sumf = 0
            sumk = 0
            for i in range(f):
                sumf = sumf + l[i]
            for i in range(k):
                sumk = sumk + l_[i]
            logV = sumf - sumk
            if flag == True:
                if logV / beta_ + cd[beta_-1] < l[k]:
                    l_k = logV / beta_ + cd[beta_-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
                    flag = False
                else:
                    l_k = l[k]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
            else:
                l_k = db_size_base + logV / beta_ + cd[beta_-1]
                if tours == 0:
                    l_.append(l_k)
                else:
                    l_[k] = l_k
                
            #print(l_)
        k_ = min (len(rk), beta)
        #last 45 norms
        sumf=0
        sumk=0
        for i in range(d):
            sumf = sumf + l[i]
        for i in range(d - 45):
            sumk = sumk + l_[i]
        logV = sumf - sumk
        for k in range(d - k_, d):
            l_k =  db_size_base + logV/k_ + rk[k + 45 - d]
            if tours == 0:
                l_.append(l_k)
            else:
                l_[k] = l_k
                #Set l[i] as l_[i] 
        for i in range(d):
            l[i] = l_[i]
    return [_ for _ in l_]


def sim_pump_d_beta(log_rr0,d,llb,cd=cd):
    
    # print(db_size_base)
    N=1
    r = d
    beta = r - llb
    #dim4f = dim4free_wrapper(theo_dim4free_fun1, beta)

    db_size_base = 0 #log(sqrt(4/3. * (beta / (beta-dim4f))))
    # print(db_size_base)
    if beta == 0:
        return log_rr0
    l = [_ for _ in log_rr0]
    #print(len(l))
    l_ = []
    
    
    for tours in range(N):
        flag = True
        for k in range(d - beta):
            l_.append(l[k])

        for k in range(d - beta ,d - 45):
            beta_ = d - k
            f = d
            sumf = 0
            sumk = 0
            for i in range(f):
                sumf = sumf + l[i]
            for i in range(k):
                sumk = sumk + l_[i]
            logV = sumf - sumk
            if flag == True:
                if logV / beta_ + cd[beta_-1] < l[k]:
                    l_k = logV / beta_ + cd[beta_-1]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
                    flag = False
                else:
                    l_k = l[k]
                    if tours == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
            else:
                l_k = db_size_base + logV / beta_ + cd[beta_-1]
                if tours == 0:
                    l_.append(l_k)
                else:
                    l_[k] = l_k
            if k == d - beta:
                return l_[d-beta]
                
            #print(l_)
        k_ = min (len(rk), beta)
        #last 45 norms
        sumf=0
        sumk=0
        for i in range(d):
            sumf = sumf + l[i]
        for i in range(d - 45):
            sumk = sumk + l_[i]
        logV = sumf - sumk
        for k in range(d - k_, d):
            l_k =  db_size_base + logV/k_ + rk[k + 45 - d]
            if tours == 0:
                l_.append(l_k)
            else:
                l_[k] = l_k
                #Set l[i] as l_[i] 
            if k == d - beta:
                return l_[d-beta]
     




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
             



def gen_pump_figures(dir):
    file_name = "../../PnjBKZcost_test_result_test6_pumpData_theo_down_sieve.txt"
    log_GS_lengths,GS_lengths,_ = read_file(file_name) #read files
    d = len(log_GS_lengths[0])
    

    # try:
    #     shutil.rmtree(dir, ignore_errors=False, onerror=None)
    # except FileNotFoundError:
    #     pass
    
    i = 0
    tours = 2
    log_GS_lengths_dic = {}
    log_GS_lengths_dic[(0,d,0,0)] = log_GS_lengths[0]
    for llb in range(d-50,0,-5):
        T_tmp = []
        RAM_tmp = []
        beta = d - llb
        f = dim4free_wrapper(theo_dim4free_fun1, beta)
        if beta - f < 131:
            for t in range(tours):
                if 2*tours*(i+1) >= len(log_GS_lengths):
                    break
                # show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t],0,d,0,t)
                show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t],log_GS_lengths[2*tours*i + 2*t + 1],llb,beta,f,t)

                log_GS_lengths_dic[(llb,beta,f,t)] = log_GS_lengths[2*tours*i + 2*t + 1]

                sim_gs = pump_simulator(log_GS_lengths[2*tours*i + 2*t],d,llb,d)
                # sim_gs = pump_simulator_test(log_GS_lengths[2*tours*i + 2*t],d,llb,d)
    
                show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t+1],sim_gs,llb,beta,f,2)
    
            i += 1

        

    return log_GS_lengths_dic




def gen_pump_figures_more_rounds(dir):

    file_name = "../../PnjBKZcost_test_result_test6_pumpData_g6k_down_sieve.txt"
    log_GS_lengths,GS_lengths,_ = read_file(file_name) #read files
    d = len(log_GS_lengths[0])
    

    try:
        shutil.rmtree(dir, ignore_errors=False, onerror=None)
    except FileNotFoundError:
        pass
    
    i = 0
    tours = 2
    log_GS_lengths_dic = {}
    log_GS_lengths_dic[(0,d,0,0)] = log_GS_lengths[0]
    for llb in range(d-50,0,-5):
        T_tmp = []
        RAM_tmp = []
        beta = d - llb
        f = dim4free_wrapper(theo_dim4free_fun1, beta)
        if beta - f < 131:
            for t in range(tours):
                if (tours+1)*(i+1) >= len(log_GS_lengths):
                    break
                # show_pump_gs_figure(dir,log_GS_lengths[2*tours*i + 2*t],0,d,0,t)
                show_pump_gs_figure(dir,log_GS_lengths[(tours+1)*i],log_GS_lengths[(tours+1)*i + t +1],llb,beta,f,t+1)

                # log_GS_lengths_dic[(llb,beta,f,t+1)] = log_GS_lengths[tours*i + t]

                sim_gs = pump_simulator(log_GS_lengths[(tours+1)*i],d,llb,d)
    
                show_pump_gs_figure(dir,log_GS_lengths[(tours+1)*i+t+1],sim_gs,llb,beta,f,11+t)

                # print(11+t)
    
            i += 1

        

    # return log_GS_lengths_dic


