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

    ldc_param = "MATZOV22": list decoding method in [MATZOV22]
                "AGPS20"(default): list decoding method in [AGPS20]

    cal_ee = "avg_sigma": martin's primal usvp + two step
             "chi" #chi-square + probabilistic + two-step
'''

'''
#######################################
#Fixed parameters
method = 2
cost_model = 2
worst_case = True
#------------------------------------

lwechal_est(method,  cost_model, worst_case)
'''


#######################################
#Fixed parameters
method = 2
cost_model = 1
ldc_param = "MATZOV22"
worst_case = False
cal_ee = "chi"
goal_min_cost = "gate_min"
#------------------------------------

lwechal_est(method,  cost_model, worst_case, ldc_param =  ldc_param, cal_ee = cal_ee, goal_min_cost = goal_min_cost)