from timu import timu

load("../framework/est/NIST-est.sage")


method = 2
worst_case = False
ldc_param = "MATZOV22"
goal_min_cost = "gate_min"
cal_ee = "chi"






saitinumber = 6
n, q, m, _, _, target_norm = timu(saitinumber)
eta1 = 2
eta2 = 2
D_s =build_centered_binomial_law(eta1) 
D_e = build_centered_binomial_law(eta2) 

dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
svp_estimate_attack(method=method, ldc_param =  ldc_param, cal_ee = cal_ee, goal_min_cost = goal_min_cost, worst_case = worst_case, dvol = dvol, dim_ = dim_)



saitinumber = 7
n, q, m, _, _, target_norm = timu(saitinumber)
eta2 = 2
t = 120
p = 256
D_s = {1: t/(2.*p), -1: t/(2.*p), 0: (p -t)/p} #SparseTernary: t/2 entries of 1,  t/2 entries of -1, rest is 0.
D_e = build_centered_binomial_law(eta2) 

dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
svp_estimate_attack(method=method, ldc_param =  ldc_param, cal_ee = cal_ee, goal_min_cost = goal_min_cost, worst_case = worst_case, dvol = dvol, dim_ = dim_)



saitinumber = 8
n, q, m, _, _, target_norm = timu(saitinumber)
n = 2*n
m = 2*m
eta1 = 3
eta2 = 2    
D_s =build_centered_binomial_law(eta1) 
D_e = build_centered_binomial_law(eta2) 

dim_, dvol = initialize_from_LWE_instance(n, q, m, D_e, D_s)
svp_estimate_attack(method=method, ldc_param =  ldc_param, cal_ee = cal_ee, goal_min_cost = goal_min_cost, worst_case = worst_case, dvol = dvol, dim_ = dim_)

