#include "../framework/est.h"

//./lwe-chal-est-enumbs-sec.sh | tee enumbs_result/lwe-chal-est-enumbs-sec.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 300;
    params->cost_model = 2; //sec model;
    params->enumbs_add_strategy = 2;
    params->J = 15;
    params->max_loop = 30;
    params->progressive_sieve = true;
    // params->verbose = false;
    // params->threads = 50;
    params->threads = 50;
    // params->threads = 1;
    // params->enumbs_slope_prec = 0; //1e-6; //1e-4;
    // params->enumbs_G_prec = 0;//1e-5; //1e-3;
    params->debug = false;

    vector<LWEchal> unsolved_lwechallenges = load_unsolved_lwechallenges();
    for(int i = 0; i < int(unsolved_lwechallenges.size()); i++){
        lwechal_est(unsolved_lwechallenges[i], params);
    }
    return 1;
}