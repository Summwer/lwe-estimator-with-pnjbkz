# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time
import math
import matplotlib.pyplot as plt
import os
import numpy as np

from collections import OrderedDict # noqa
from math import pi,exp,log,sqrt,lgamma


#read GS length 
def read_gs_lengths(filename):
    f = open(filename,'r')
    GS_lengths = f.read()
    GS_lengths = GS_lengths.split('\n')
    GS_lengths.remove('')
    log_GS_lengths = []
    for i in range(len(GS_lengths)):
        gs = GS_lengths[i]
        gs = gs.split(' ')
        gs.remove('')
        GS_lengths[i] = [float(_) for _ in gs]
        log_GS_lengths.append([log(_) for _ in GS_lengths[i]])
    return log_GS_lengths,GS_lengths

def read_blocksizes(file_name):
    f = open(file_name,'r')
    blocksizes = f.read()
    f.close()
    blocksizes = blocksizes.split(' ')
    blocksizes.remove('')
    blocksizes = [int(_) for _ in blocksizes]
    return blocksizes

def read_file(n,alpha_,jump,Pumpdown,Blocksize,Tours,test_number,d):
    dir = "gs-lengths-simulator/"+"n=%d,alpha=%s,jump=%d,Pumpdown=%s,Blocksize=%d,Tours=%d,test_number=%dd=%d,/" %(n,alpha_,jump,Pumpdown,Blocksize,Tours,test_number,d)

    log_GS_lengths,GS_lengths = read_gs_lengths(dir+"rr_set.txt")
    
    return log_GS_lengths,GS_lengths,dir

#When jump=1, pnj-bkz simulator degenerates to CN11 simulator
def pnjBKZ_simulator(log_rr0,beta,N,d,jump):
    l = copy.deepcopy(log_rr0)
    d = len(l)
    if N == 0 or beta == 0:
        return deepcopy(l)
    l_ = [0 for _ in range(d)]
    J = jump
    randomSquaredAverages = [2.95747208 , 2.905419231 , 2.808209731 , 2.704122291 , 2.589140986 , 2.473617892 , 2.363752978 , 2.254288677 , 2.160395606 , 2.056417068 , 1.962571394 , 1.86329581 , 1.781471788 , 1.691206302 , 1.607855519 , 1.526258286 , 1.45048576 , 1.377169955 , 1.304800173 , 1.238714637 , 1.172164783 , 1.11055388 , 1.04872671 , 0.9933346386 , 0.9385826564 , 0.8858378414 , 0.8346947761 , 0.7869433841 , 0.743196743 , 0.6999110693 , 0.6612589966 , 0.620890347 , 0.5840212199 , 0.5494461008 , 0.5153120689 , 0.4851862344 , 0.4542255962 , 0.426615261 , 0.4011004324 , 0.3761683243 , 0.3521434119 , 0.3305054658 , 0.3109962405 , 0.2976205498 , 0.2897254705]
    #print(len(randomSquaredAverages))
    rk = []
    for i in range(len(randomSquaredAverages)):
        RK = log(math.sqrt(randomSquaredAverages[i]))
        rk.append(RK)
    
    k_ = min (len(rk), beta)
    
    cd = []
    for i in range(beta):
        CD = lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2
        cd.append(CD)
    #print(cd)
    print(N)
    for tours in range(N):
        flag = True
        sumf = 0.
        sumk = 0.
        ff = 0
        for k in range(d - beta):
            beta_ = min(beta, d - k)
            f = min(k - (k % J) + beta, d)
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
                
        #last 45 norms
        logV -= l_[d-46]
        for k in range(d - k_, d):
            l_[k] = logV/k_ + rk[k + 45 - d]
            
         #Set l[i] as l_[i] 
        l = copy.deepcopy(l_)
        l_ = [0 for _ in range(d)]
    return l


def CN11_simulator(log_rr0,beta,N,d):
    l = []
    for i in range(d):
        temp = log_rr0[i]
        l.append(temp)
    #print(len(l))
    l_ = []
    
    randomSquaredAverages = [2.95747208 , 2.905419231 , 2.808209731 , 2.704122291 , 2.589140986 , 2.473617892 , 2.363752978 , 2.254288677 , 2.160395606 , 2.056417068 , 1.962571394 , 1.86329581 , 1.781471788 , 1.691206302 , 1.607855519 , 1.526258286 , 1.45048576 , 1.377169955 , 1.304800173 , 1.238714637 , 1.172164783 , 1.11055388 , 1.04872671 , 0.9933346386 , 0.9385826564 , 0.8858378414 , 0.8346947761 , 0.7869433841 , 0.743196743 , 0.6999110693 , 0.6612589966 , 0.620890347 , 0.5840212199 , 0.5494461008 , 0.5153120689 , 0.4851862344 , 0.4542255962 , 0.426615261 , 0.4011004324 , 0.3761683243 , 0.3521434119 , 0.3305054658 , 0.3109962405 , 0.2976205498 , 0.2897254705]
    #print(len(randomSquaredAverages))
    rk = []
    for i in range(len(randomSquaredAverages)):
        RK = log(math.sqrt(randomSquaredAverages[i]))
        rk.append(RK)
    #print(rk)
    k_ = min (len(randomSquaredAverages), beta)
    cd = []
    for i in range(beta):
        CD = lgamma(((i+1)/2.0 + 1))*(1/(i+1))-log(pi)/2
        cd.append(CD)
    #print(cd)
    print(N)
    for j in range(N):
        flag = True
        for k in range(d):
            beta_ = min(beta, d - k)
            f = min(k + beta, d)
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
                    if j == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
                    flag = False
                else:
                    l_k = l[k]
                    if j == 0:
                        l_.append(l_k)
                    else:
                        l_[k] = l_k
            else:
                l_k = logV / beta_ + cd[beta_-1]
                if j == 0:
                    l_.append(l_k)
                else:
                    l_[k] = l_k
            #print(l_)
        #Set l[i] as l_[i]         
        for i in range(d):
            l[i] = l_[i]
    return l_

def show_gs_slope_figure(dir,log_gs_length,sim_log_gs_lengths,n,dimension,alpha_,square_error,Blocksize,Jump,N,test_number):
    plt.figure(figsize=(15, 10), dpi=100)

    t = 0
    plt.scatter([_+1 for _ in range(t,dimension)],[log_gs_length[_] for _ in range(t,dimension) ],marker="o",c='none',edgecolors='b')
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="x",c='r')
    #plt.title("square error = %f" %(square_error),x=0.9,y=0.875)
    plt.rcParams.update({'font.size':20})
    plt.legend(['Experimental','jump simulator'])
   

    plt.savefig(dir+"n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    plt.close()
    
def show_gs_slope_figure2(dir,log_gs_length,sim_log_gs_lengths,CN_11_mid_log_gs,n,dimension,alpha_,square_error1,square_error2,Blocksize,Jump,N,test_number):
    plt.figure(figsize=(15, 10), dpi=100)

    t = 0
    plt.scatter([_+1 for _ in range(t,dimension)],[log_gs_length[_] for _ in range(t,dimension) ],marker="^",c='none',edgecolors='k')
    plt.scatter([_+1 for _ in range(t,len(sim_log_gs_lengths)) ],[sim_log_gs_lengths[_] for _ in range(t,len(sim_log_gs_lengths)) ],marker="x",c='r')
    plt.scatter([_+1 for _ in range(t,len(CN_11_sim)) ],[CN_11_sim[_] for _ in range(t,len(CN_11_sim)) ],marker="o",c='none',edgecolors='b')
    plt.title("pnj-BKZ simulator prediction error = %f,\n BKZ2.0 simulator prediction error = %f" %(square_error1,square_error2))
    plt.xlabel("Index")
    plt.ylabel("log GS norm")
    plt.legend(['Experimental','pnj-BKZ simulator','BKZ2.0 simulator'])
  

    plt.savefig(dir+"Comparing figure, n=%d, alpha = %s, dimension = %d, Blocksize = %d, Jump = %d, #Current Tours= %d.png" %(n,alpha_,dimension,Blocksize,Jump,N+1))
    plt.close()    

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

    return square_error
             
        
#=========================================input==============================================

if __name__ == "__main__":
    n, d, alpha_, jump, Pumpdown, Blocksize, Tours, test_number  = 75,252,"005",9,"True",95,12,20
    Blocksizes=[95,95,95,95,95,95,95,95,95,95,95,95]

    # log_GS_lengths,GS_lengths,dir = read_file(n,alpha_,jump,Pumpdown,Blocksize,Tours,test_number,d) #read
    
    #Total_test_number represents the number of pnj-BKZ reduction experiments
    Total_test_number = 20
    GS_sum = []
    for i in range(Total_test_number):
        log_GS_lengths,GS_lengths,dir = read_file(n,alpha_,jump,Pumpdown,Blocksize,Tours,i+1,d)
        if i==0:
            GS_sum = copy.deepcopy(log_GS_lengths)
        else:
            for k in range(Tours):
                for j in range(d):
                    GS_sum[k][j] = GS_sum[k][j] + log_GS_lengths[k][j]
    
    average = GS_sum.copy()
    
    #Calculate the average GS values of multiple experiments
    for i in range(Tours):
        for j in range(d):
            average[i][j] = average[i][j]/Total_test_number
        
    log_GS_lengths = average
    
    
    #log_GS_lengths we get are 2*log||bi*||, so we do following pre-processingg: 
    for i in range(len(Blocksizes)+1):
        for j in range(d):
            log_GS_lengths[i][j] = log_GS_lengths[i][j]/2    
    
    mid_log_gs=[]
    for N in range(len(Blocksizes)):
        if N == 0:
            if Blocksizes[0]>44:
                sim_log_gs_lengths = pnjBKZ_simulator(log_GS_lengths[0], Blocksizes[0], 1, d, jump).copy()
                CN_11_sim = CN11_simulator(log_GS_lengths[0], Blocksizes[0], 1, d).copy()
            else:
                sim_log_gs_lengths = CN11_simulator(log_GS_lengths[0], Blocksizes[0], 1, d).copy()
                CN_11_sim = CN11_simulator(log_GS_lengths[0], Blocksizes[0], 1, d).copy()
            mid_log_gs = sim_log_gs_lengths.copy()
            CN_11_mid_log_gs = CN_11_sim.copy()

            square_error1 = compute_square_error(sim_log_gs_lengths,log_GS_lengths[N+1],1)
            square_error2 = compute_square_error(CN_11_sim,log_GS_lengths[N+1],1)
            print("\nOur square error: %f, CN11 square error: %f" %(square_error1,square_error2))
            show_gs_slope_figure(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
            show_gs_slope_figure2(dir,log_GS_lengths[N+1],sim_log_gs_lengths,CN_11_mid_log_gs,n,d,alpha_,square_error1,square_error2,Blocksizes[N],jump,N,test_number)
        else:
            if Blocksizes[N]>44:
                sim_log_gs_lengths = pnjBKZ_simulator(mid_log_gs, Blocksizes[N], 1, d, jump).copy()
                CN_11_sim = CN11_simulator(CN_11_mid_log_gs, Blocksizes[N], 1, d).copy()
            else:
                sim_log_gs_lengths = CN11_simulator(mid_log_gs, Blocksizes[N], 1, d).copy()
                CN_11_sim = CN11_simulator(CN_11_mid_log_gs, Blocksizes[N], 1, d).copy()
            mid_log_gs = sim_log_gs_lengths.copy()
            CN_11_mid_log_gs = CN_11_sim.copy()
            
            square_error1 = compute_square_error(sim_log_gs_lengths,log_GS_lengths[N+1],1)
            square_error2 = compute_square_error(CN_11_sim,log_GS_lengths[N+1],1)
            print("\nOur square error: %f, CN11 square error: %f" %(square_error1,square_error2))
            show_gs_slope_figure(dir,log_GS_lengths[N+1],sim_log_gs_lengths,n,d,alpha_,square_error1,Blocksizes[N],jump,N,test_number)
            show_gs_slope_figure2(dir,log_GS_lengths[N+1],sim_log_gs_lengths,CN_11_mid_log_gs,n,d,alpha_,square_error1,square_error2,Blocksizes[N],jump,N,test_number)
            
        


    
