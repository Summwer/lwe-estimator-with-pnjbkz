#include "../framework/est.h"


// ./lwe-chal-est-enumbs-gate.sh | tee lwe-chal-est-enumbs-gate.log


int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 1; //sec model;
    params->debug = false;
    params->progressive_sieve =  true; 
  
    vector<LWEchal> unsolved_lwechallenges = load_unsolved_lwechallenges();
    for(int i = 0; i < int(unsolved_lwechallenges.size()); i++){
        lwechal_est(unsolved_lwechallenges[i], params);
    }
    return 1;
}