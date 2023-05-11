#include "../framework/est.h"

//./lwe-chal-est-enumbs-sec.sh | tee enumbs_result/lwe-chal-est-enumbs-sec.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 2; //sec model;
    params->progressive_sieve = true;
    params->threads = 10;
    params->verbose = true;
    params->debug = false;
    params->worst_case = true;
    params->J = 8; 

    vector<pair<int,double>> lwes;

    //low_dim_lwechallenge_est. 
    // lwes = {{40,0.035}, 
    //{40, 0.025}, {45, 0.020},
    // // lwes= { {50, 0.015}, {55, 0.010}, {60, 0.010}, {70, 0.005}, {75, 0.005}};
    // for(int i = 0; i < int(lwes.size());i++){
    //     int n = lwes[i].first;
    //     double alpha  = lwes[i].second;
    //     lwechal_est(n, alpha, params);
    // }
    
    // //solve_lwechal_est
    // lwes =  {{45,0.030}, {50, 0.025}, {55, 0.020}, {60, 0.015}, {85, 0.005}, {90, 0.005}};
    // for(int i = 0; i < int(lwes.size());i++){
    //     int n = lwes[i].first;
    //     double alpha  = lwes[i].second;
    //     lwechal_est(n, alpha, params);
    // }

    // lwes = {{40,0.025}};
    lwes = {{40,0.045}, {45, 0.035}, {50, 0.030}, {55, 0.025}, {60, 0.020}, {65, 0.015}, {75, 0.010}, {95, 0.005}};
    // lwes = {{95, 0.005}};
    for(int i = 0; i < int(lwes.size());i++){
        int n = lwes[i].first;
        double alpha  = lwes[i].second;
        lwechal_est(n, alpha, params);
    }
    
    return 1;
}
