#include "est.h"

void gsa_est(int dim, FP_NR<FT> dvol, Params* params){
    printf_input(dim,dvol);
    printf("Generate gs-lengths by GSA assumption.");
    vector<double> l = gen_simulated_gso(dim, dvol);
    if(params->method == 1)
        call_enumbs(l,params);
}



void lwechal_est(LWEchal lwechal, Params* params){
    int n = lwechal.n, dim = lwechal.dim;
    double alpha = lwechal.alpha;
    FP_NR<FT> dvol = lwechal.dvol;
    
    printf("TU LWE Challenge, n = %3d, alpha = %3.4f\n", n,alpha);
    gsa_est(dim, dvol, params);
}

// solved lwe challenges: (n,alpha,dim,dvol)
vector<LWEchal> load_solved_lwechallenges(){
    return  {
        {45, 0.030, 202, 358.0663365218434},
        {50, 0.025, 220, 412.4659700420454},
        {55, 0.020, 231, 454.6342882113398},
        {60, 0.015, 242, 516.7000963118265},
        {85, 0.005, 287, 756.4334688012609},
        {90, 0.005, 312, 834.0984244884119}
    };
}


   

// unsolved lwe challenges: (n,alpha,dim,dvol)
vector<LWEchal> load_unsolved_lwechallenges(){
    return {
        {40, 0.045,	195,	302.1993617},
        // {45, 0.035,	211,	390.6236363},
        // {50, 0.030,	228,	400.4076907},
        // {55, 0.025,	241,	439.9769224},
        // {60, 0.020,	254,	494.0253108},
        // {65, 0.015,	262,	557.6405653},
        // {75, 0.010,	281,	637.6057085},
        {95, 0.005,	323,	836.9696072}

    };
}

