load("../framework/est/lwechal-est.sage")

######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3: bkz-only with fixed blocksize
             4: default g6k
             
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
'''
#######################################
#Fixed parameters
method = 3
cost_model = 2
worst_case = True
#------------------------------------

lwechals = [  (40,0.025),(80,0.005)]
lwechal_est(lwechals, method,  cost_model, worst_case)
