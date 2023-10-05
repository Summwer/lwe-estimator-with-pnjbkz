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
#cost_model = 1
worst_case = False
cost_model = 2
cal_ee = "chi"
goal_min_cost = "gate_min"
#------------------------------------


lwechals = [ (40,0.005),(45,0.005),(50,0.005), (55,0.005), (40,0.015), (45,0.010)] 
lwechal_est(lwechals, method,  cost_model, worst_case, cal_ee = cal_ee, goal_min_cost = goal_min_cost)


cost_model = 1
ldc_param = "AGPS20"
lwechals = [  (40,0.025),(80,0.005)]
lwechal_est(lwechals, method,  cost_model, worst_case, ldc_param =  ldc_param, cal_ee = cal_ee, goal_min_cost = goal_min_cost)

ldc_param = "MATZOV22"
lwechal_est(lwechals, method,  cost_model, worst_case, ldc_param =  ldc_param, cal_ee = cal_ee, goal_min_cost = goal_min_cost)