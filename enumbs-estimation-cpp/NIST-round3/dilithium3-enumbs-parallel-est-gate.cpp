#include "../framework/est.h"


//./dilithium3-enumbs-parallel-est-gate.sh | tee enumbs_result/dilithium3-enumbs-parallel-est-gate.log


int main(){

    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 1000; // maximal selected blocksize value
    // params->J = 15; //maximal jump value
    params->threads = 50;

    // Dilithium-III round-3 parameters
    printf("============= Dilithium-III\n");
    int dim =  3841;
    FP_NR<FT> dvol = 31317.16147360077;
    gsa_est(dim, dvol, params);

    return 1;

}






