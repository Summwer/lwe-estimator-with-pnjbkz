#include "../framework/est.h"


//./dilithium2-enumbs-est-gate.sh | tee enumbs_result/dilithium2-enumbs-est-gate.log

int main(){

    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 700; // maximal selected blocksize value
    // params->J = 15; //maximal jump value

    printf("============= Dilithium-II\n");
    int dim =  2817;
    FP_NR<FT> dvol = 21814.858106487554;
    gsa_est(dim, dvol, params);

    return 1;

}



