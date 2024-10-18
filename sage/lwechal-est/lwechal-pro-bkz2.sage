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

    :cumG : False: [PV21] estimate
            True: cumulated gate estimate adaptive to [PV21]
'''

'''
#######################################
#Fixed parameters
method = 1
cost_model = 1
worst_case = True
#------------------------------------

lwechal_est(method,  cost_model, worst_case)
'''


#######################################
#Fixed parameters
method = 1
cost_model = 1
ldc_param = "MATZOV22"
worst_case = False
cumG = True
#------------------------------------

lwechals = [(40,0.020), (40,0.025), (40,0.030), (40,0.035), (40,0.040), 
            (45,0.020), (45,0.025), (45,0.030), (45,0.035), (45,0.040), 
            (50,0.020), (50,0.025), (50,0.030), (50,0.035), (50,0.040), 
            (55,0.020), (55,0.025), (55,0.030), (55,0.035), (55,0.040), 
            (60,0.020), (60,0.025), (60,0.030), (60,0.035), (60,0.040), 
            (65,0.020), (65,0.025), (65,0.030), (65,0.035), (65,0.040),
            (70,0.020), (70,0.025), (70,0.030), (70,0.035), (70,0.040), 
            (75,0.020), (75,0.025), (75,0.030), (75,0.035), (75,0.040), 
            (80,0.020), (80,0.025), (80,0.030), (80,0.035), (80,0.040)] 
lwechal_est(lwechals, method,  cost_model, worst_case, cumG = cumG)