
#include "../framework/est.h"


//params input in main function
//argv[0]: implemented file name
//argv[1]: jump value J
//argv[2]: max_loop value
//argv[3]: cost model: 1.gate model 2.practical sec model
//argv[4]: maximal dimension in enumeration
//argv[5]: threads number
//argv[6]: practical_pump_d4f
//argv[7]: start beta value.
//argv[8]: est model in dsvp_prediction for last pump
int main(int argc,char **argv){
   Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->threads = 60;
    params->cost_model = 1; //sec model;
    // params->progressive_sieve = true;
    params->verbose = true;
    params->debug = false;
    params->worst_case = false;
    params->method = atoi(argv[1]); //1:enumbs;2:bssa
    params->gap = 1;
    
    params->J_gap = 1;
    params->enumbs_G_prec = 0.001;
    params->enumbs_slope_prec = 1e-6;
    params->max_loop = atoi(argv[3]); 

    params->J = min(atoi(argv[4]),100); 
    params->max_dim = 1500; 
    params->min_G_prec = 0.001;
    if(atoi(argv[2]) == 1)
        params->list_decoding = "agps20"; //"matzov22"
    if(atoi(argv[2]) == 2)
        params->list_decoding = "matzov22"; //"matzov22"
    // params->bssa_tradion = true;  
    

    int n, m, q,eta;
    map<int,double> D_s,D_e;
    LWEchal* lwechal;




    // timu 6 parameters
    printf("============= timu 6\n");
    n = 256, m = 256, q = 3329;
    D_s = build_centered_binomial_law(2);
    D_e = D_s;
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params->verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);
    
    // timu 7 parameters
    printf("============= timu 7\n");
    n = 256, m = 256, q = 3329;
    D_s={};
    double t = 120, p = 256;
    D_s[-1] = t/(2*p);
    D_s[1] = t/(2*p);
    D_s[0] = (p -t)/p;
    
    D_e = build_centered_binomial_law(2);
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params->verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);



    // timu8 parameters
    printf("============= timu8\n");
    n = 512, m = 512, q = 3329;
    D_s = build_centered_binomial_law(3);
    D_e = build_centered_binomial_law(2);
    lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params->verbose);
    gsa_est(lwechal->dim, lwechal->dvol, params);

    

    return 1;
}
