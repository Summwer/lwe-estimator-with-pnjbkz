

# This file was *autogenerated* from the file lwechal-pro-bkz-matzov.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_40 = Integer(40); _sage_const_0p025 = RealNumber('0.025'); _sage_const_80 = Integer(80); _sage_const_0p005 = RealNumber('0.005')
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
cost_model = 2
worst_case = True
#------------------------------------

lwechal_est(method,  cost_model, worst_case)
'''


#######################################
#Fixed parameters
method = _sage_const_1 
cost_model = _sage_const_1 
worst_case = False
cumG = True
ldc_param = "MATZOV22"
#------------------------------------

lwechals = [  (_sage_const_40 ,_sage_const_0p025 ),(_sage_const_80 ,_sage_const_0p005 )]
lwechal_est(lwechals, method,  cost_model, worst_case, cumG = cumG, ldc_param = ldc_param)

