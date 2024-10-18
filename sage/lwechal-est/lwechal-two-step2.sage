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

lwechals = [(40,0.015), (40,0.020), (40,0.025), (40,0.030), (40,0.035), (40,0.040), 
            (45,0.015), (45,0.020), (45,0.025), (45,0.030), (45,0.035), (45,0.040), 
            (50,0.010), (50,0.015), (50,0.020), (50,0.025), (50,0.030), (50,0.035), (50,0.040), 
            (55,0.010), (55,0.015), (55,0.020), (55,0.025), (55,0.030), (55,0.035), (55,0.040), 
            (60,0.010), (60,0.015), (60,0.020), (60,0.025), (60,0.030), (60,0.035), (60,0.040), 
            (65,0.010), (65,0.015), (65,0.020), (65,0.025), (65,0.030), (65,0.035), (65,0.040),
            (70,0.010), (70,0.015), (70,0.020), (70,0.025), (70,0.030), (70,0.035), (70,0.040), 
            (75,0.005), (75,0.010), (75,0.015), (75,0.020), (75,0.025), (75,0.030), (75,0.035), (75,0.040), 
            (80,0.005), (80,0.010), (80,0.015), (80,0.020), (80,0.025), (80,0.030), (80,0.035), (80,0.040),
            (85,0.005),
            (90,0.005), 
            (100,0.005)]

             
            
            
             
            
            
            
lwechals = [(60,0.005), (70,0.005), (85,0.005), (95,0.005)]; #[(60,0.010)]          

lwechals = [(90,0.010), (95,0.010)]; #[(60,0.010)]       

cost_model = 2

lwechal_est(lwechals, method,  cost_model, worst_case, cal_ee = cal_ee, goal_min_cost = goal_min_cost)