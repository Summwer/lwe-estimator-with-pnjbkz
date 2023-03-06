#include "../framework/est.h"


// ./lwe-chal-est-enumbs-gate.sh | tee lwe-chal-est-enumbs-gate.log


int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 2; //sec model;
    // params->debug = true;
    params->method = 2; 
    params->max_loop = 10;
    params->verbose = false;

    vector<LWEchal> unsolved_lwechallenges = load_unsolved_lwechallenges();
    for(int i = 0; i < int(unsolved_lwechallenges.size()); i++){
        lwechal_est(unsolved_lwechallenges[i], params);
    }
    return 1;
}