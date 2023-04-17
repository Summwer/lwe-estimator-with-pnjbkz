#include "../framework/est.h"

//./lwe-chal-est-enumbs-sec.sh | tee enumbs_result/lwe-chal-est-enumbs-sec.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 2; //sec model;
    params->progressive_sieve = true;
    params->threads = 15;
    params->debug = false;
    params->worst_case = true;

    // vector<LWEchal> unsolved_lwechallenges = load_unsolved_lwechallenges();
    
    // for(int i = 0; i < int(unsolved_lwechallenges.size()); i++){
    //     lwechal_est(unsolved_lwechallenges[i], params);
    // }

    params->J = 8; 
    vector<LWEchal> low_dim_lwechallenges = load_low_dim_lwechallenges();
    for(int i = 0; i < int(low_dim_lwechallenges.size()); i++){
        lwechal_est(low_dim_lwechallenges[i], params);
    }
    return 1;
}