#include "est.h"


//dvol = sum(log2(||b_i^*||/sigma))
void gsa_est(int dim, FP_NR<FT> dvol, Params* params){
    printf_input(dim,dvol);
    printf("Generate gs-lengths by GSA assumption.");
    vector<double> l = gen_simulated_gso(dim, dvol);
    switch(params->method){
        case 1:
            call_enumbs(l,params);
            break;
        case 2:
            call_bssa(l,params,50, dim);
            // call_bssa(l,params,50, int(0.9*(double)dim));
            break;
        default:
            cout<<"Tere's no method named: "<<params->method<<endl;
    }
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


// low-dim lwe challenges: (n,alpha,dim,dvol)
vector<LWEchal> load_low_dim_lwechallenges(){
    return  {
        // {40, 0.035,	188,	327.7388246557668}
        {40, 0.025,	172,	331.9735338747315},
        {45, 0.020,	185,	373.4658972674150},
        {50, 0.015,	194,	415.6552752456852},
        {55, 0.010,	205,	495.0168620849964},
        {60, 0.010,	222,	522.7192487542407},
        {70, 0.005,	235,	641.7748006815514},
        // {75, 0.005,	252,	678.7288625607595}

    };
}



   

// unsolved lwe challenges: (n,alpha,dim,dvol)
vector<LWEchal> load_unsolved_lwechallenges(){
    return {
        {40, 0.045,	195,	302.1993617},
        {45, 0.035,	211,	357.0995642},
        {50, 0.030,	228,	400.4076907},
        {55, 0.025,	241,	439.9769224},
        {60, 0.020,	254,	494.0253108},
        {65, 0.015,	262,	557.6405653},
        {75, 0.010,	281,	637.6057085},
        {95, 0.005,	323,	836.9696072}
    };
}

