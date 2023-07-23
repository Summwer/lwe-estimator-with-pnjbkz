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
#ldc_param = "MATZOV22" #change list decoding parameters

#------------------------------------
#The original leaky lwe estimator
kyber_est(method)
dilithium_est(method)
frodo_est(method)
