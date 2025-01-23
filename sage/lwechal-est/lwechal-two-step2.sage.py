

# This file was *autogenerated* from the file lwechal-two-step2.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_40 = Integer(40); _sage_const_0p015 = RealNumber('0.015'); _sage_const_0p020 = RealNumber('0.020'); _sage_const_0p025 = RealNumber('0.025'); _sage_const_0p030 = RealNumber('0.030'); _sage_const_0p035 = RealNumber('0.035'); _sage_const_0p040 = RealNumber('0.040'); _sage_const_45 = Integer(45); _sage_const_50 = Integer(50); _sage_const_0p010 = RealNumber('0.010'); _sage_const_55 = Integer(55); _sage_const_60 = Integer(60); _sage_const_65 = Integer(65); _sage_const_70 = Integer(70); _sage_const_75 = Integer(75); _sage_const_0p005 = RealNumber('0.005'); _sage_const_80 = Integer(80); _sage_const_85 = Integer(85); _sage_const_90 = Integer(90); _sage_const_100 = Integer(100); _sage_const_95 = Integer(95)
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
method = _sage_const_2 
#cost_model = 1
worst_case = False
cost_model = _sage_const_2 
cal_ee = "chi"
goal_min_cost = "gate_min"
#------------------------------------

lwechals = [(_sage_const_40 ,_sage_const_0p015 ), (_sage_const_40 ,_sage_const_0p020 ), (_sage_const_40 ,_sage_const_0p025 ), (_sage_const_40 ,_sage_const_0p030 ), (_sage_const_40 ,_sage_const_0p035 ), (_sage_const_40 ,_sage_const_0p040 ), 
            (_sage_const_45 ,_sage_const_0p015 ), (_sage_const_45 ,_sage_const_0p020 ), (_sage_const_45 ,_sage_const_0p025 ), (_sage_const_45 ,_sage_const_0p030 ), (_sage_const_45 ,_sage_const_0p035 ), (_sage_const_45 ,_sage_const_0p040 ), 
            (_sage_const_50 ,_sage_const_0p010 ), (_sage_const_50 ,_sage_const_0p015 ), (_sage_const_50 ,_sage_const_0p020 ), (_sage_const_50 ,_sage_const_0p025 ), (_sage_const_50 ,_sage_const_0p030 ), (_sage_const_50 ,_sage_const_0p035 ), (_sage_const_50 ,_sage_const_0p040 ), 
            (_sage_const_55 ,_sage_const_0p010 ), (_sage_const_55 ,_sage_const_0p015 ), (_sage_const_55 ,_sage_const_0p020 ), (_sage_const_55 ,_sage_const_0p025 ), (_sage_const_55 ,_sage_const_0p030 ), (_sage_const_55 ,_sage_const_0p035 ), (_sage_const_55 ,_sage_const_0p040 ), 
            (_sage_const_60 ,_sage_const_0p010 ), (_sage_const_60 ,_sage_const_0p015 ), (_sage_const_60 ,_sage_const_0p020 ), (_sage_const_60 ,_sage_const_0p025 ), (_sage_const_60 ,_sage_const_0p030 ), (_sage_const_60 ,_sage_const_0p035 ), (_sage_const_60 ,_sage_const_0p040 ), 
            (_sage_const_65 ,_sage_const_0p010 ), (_sage_const_65 ,_sage_const_0p015 ), (_sage_const_65 ,_sage_const_0p020 ), (_sage_const_65 ,_sage_const_0p025 ), (_sage_const_65 ,_sage_const_0p030 ), (_sage_const_65 ,_sage_const_0p035 ), (_sage_const_65 ,_sage_const_0p040 ),
            (_sage_const_70 ,_sage_const_0p010 ), (_sage_const_70 ,_sage_const_0p015 ), (_sage_const_70 ,_sage_const_0p020 ), (_sage_const_70 ,_sage_const_0p025 ), (_sage_const_70 ,_sage_const_0p030 ), (_sage_const_70 ,_sage_const_0p035 ), (_sage_const_70 ,_sage_const_0p040 ), 
            (_sage_const_75 ,_sage_const_0p005 ), (_sage_const_75 ,_sage_const_0p010 ), (_sage_const_75 ,_sage_const_0p015 ), (_sage_const_75 ,_sage_const_0p020 ), (_sage_const_75 ,_sage_const_0p025 ), (_sage_const_75 ,_sage_const_0p030 ), (_sage_const_75 ,_sage_const_0p035 ), (_sage_const_75 ,_sage_const_0p040 ), 
            (_sage_const_80 ,_sage_const_0p005 ), (_sage_const_80 ,_sage_const_0p010 ), (_sage_const_80 ,_sage_const_0p015 ), (_sage_const_80 ,_sage_const_0p020 ), (_sage_const_80 ,_sage_const_0p025 ), (_sage_const_80 ,_sage_const_0p030 ), (_sage_const_80 ,_sage_const_0p035 ), (_sage_const_80 ,_sage_const_0p040 ),
            (_sage_const_85 ,_sage_const_0p005 ),
            (_sage_const_90 ,_sage_const_0p005 ), 
            (_sage_const_100 ,_sage_const_0p005 )]

             
            
            
             
            
            
            
lwechals = [(_sage_const_60 ,_sage_const_0p005 ), (_sage_const_70 ,_sage_const_0p005 ), (_sage_const_85 ,_sage_const_0p005 ), (_sage_const_95 ,_sage_const_0p005 )]; #[(60,0.010)]          

lwechals = [(_sage_const_90 ,_sage_const_0p010 ), (_sage_const_95 ,_sage_const_0p010 )]; #[(60,0.010)]       

cost_model = _sage_const_2 

lwechal_est(lwechals, method,  cost_model, worst_case, cal_ee = cal_ee, goal_min_cost = goal_min_cost)
