import time
load("../framework/est/NIST-est.sage")

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
#Fixed parameters

method = 1
cost_model = 1
worst_case = False

#------------------------------------


dilithium_est(method,  cost_model, worst_case)

kyber_est(method,  cost_model, worst_case)

frodo_est(method,  cost_model, worst_case)