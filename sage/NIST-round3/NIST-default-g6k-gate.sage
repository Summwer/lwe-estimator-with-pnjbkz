import time
load("../framework/est/NIST-est.sage")

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
             "fixed-blocksize bkz-only": call fixed-blocksize bkz-only model to estimate security.
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################
#Fixed parameters

method = 6
progressive_sieve = True
cost_model = 1

#------------------------------------


dilithium_est(method,  cost_model, progressive_sieve = progressive_sieve)

kyber_est(method,  cost_model, progressive_sieve = progressive_sieve)

frodo_est(method,  cost_model, progressive_sieve = progressive_sieve)