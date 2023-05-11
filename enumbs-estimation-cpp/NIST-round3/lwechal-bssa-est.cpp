#include "../framework/est.h"


// ./lwe-chal-est-enumbs-gate.sh | tee lwe-chal-est-enumbs-gate.log


int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 2; //sec model;
    // params->debug = true;
    params->method = 2; 
    params->J = 8;
    params->verbose = false;
    params->worst_case = true;
    params->progressive_sieve = true;
    // params->max_RAM = 43.58; //1.5T = 43.58

    vector<pair<int,double>> lwes;
    lwes = {{40, 0.025}, {45, 0.020}, {50, 0.015}, {55, 0.010}, {60, 0.010}, {70, 0.005}, {75, 0.005}, {40,0.035}};
    for(int i = 0; i < int(lwes.size());i++){
        int n = lwes[i].first;
        double alpha  = lwes[i].second;
        lwechal_est(n, alpha, params);
    }

    lwes = {{40,0.045}, {45, 0.035}, {50, 0.030}, {55, 0.025}, {60, 0.020}, {65, 0.015}, {75, 0.010}, {95, 0.005}};
    for(int i = 0; i < int(lwes.size());i++){
        int n = lwes[i].first;
        double alpha  = lwes[i].second;
        lwechal_est(n, alpha, params);
    }

    return 1;
}