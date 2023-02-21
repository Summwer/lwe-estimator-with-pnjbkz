load("../framework/est/lwechal-est.sage")

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
             
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################
#Fixed parameters
method = 5
cost_model = 2
#------------------------------------

unsolved_lwe_challenge_est(method,  cost_model)