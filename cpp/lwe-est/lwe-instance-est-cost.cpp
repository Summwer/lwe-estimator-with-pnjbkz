
#include "../framework/est.h"





// params input in main function
// argv[0]: implemented file name
// argv[1]: method -- 1:enumba; 2: bssa
int main(int argc,char **argv){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->threads = 60;
    params->cost_model = 1; //sec model;
    // params->progressive_sieve = true;
    params->verbose = false;
    params->debug = false;
    params->worst_case = false;
    params->method = atoi(argv[1]); //1:enumbs;2:bssa
    params->gap = 1;
    
    params->J_gap = 1;
    params->enumbs_G_prec = 0.001;
    params->enumbs_slope_prec = 1e-6;
    params->max_loop = atoi(argv[3]); 
    params->practical_pump_d4f  = 3; //pump_d4f = default g6k.
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

    
    //input distribution: discrete Gaussian distribution
    

    double sigma = 1;
    int t = 1600;
    for( int log2q = 30 ; log2q <= 60; log2q += 1){
        // double log_q = (double) log2q / log2(M_E); //pow(2,log2q);
        for(n = 450; n <= 600; n += 50){
            for(int success_count = 10; success_count <= 10; success_count ++ ){
                params->succ_prob = 0.1 * success_count;
                m = n;
                D_s=build_Gaussian_law(sigma, t);
                D_e=build_Gaussian_law(sigma, t);
                printf("============= lwe-instance-est\n");
                lwechal = gen_LWE_instance_parameter( n, log2q, m, D_e, D_s, params->verbose);
                gsa_est(lwechal->dim, lwechal->dvol, params);

            }
    }}

    return 1;
}


