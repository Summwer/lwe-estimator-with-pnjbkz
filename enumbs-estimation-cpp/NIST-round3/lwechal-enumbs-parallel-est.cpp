#include "../framework/est.h"

//./lwe-chal-est-enumbs-sec.sh | tee enumbs_result/lwe-chal-est-enumbs-sec.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 2; //sec model;
    // params->J = 20;
    // params->max_loop = 20;
    // params->verbose = false;
    // params->threads = 50;
    params->threads = 50;
    // params->threads = 1;
    // params->enumbs_prec = 1e3;
    // params->debug = true;

    vector<LWEchal> unsolved_lwechallenges = load_unsolved_lwechallenges();
    for(int i = 0; i < int(unsolved_lwechallenges.size()); i++){
        lwechal_est(unsolved_lwechallenges[i], params);
    }
    return 1;
}