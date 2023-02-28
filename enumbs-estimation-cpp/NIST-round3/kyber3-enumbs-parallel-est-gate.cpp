
#include "../framework/est.h"

//./kyber3-enumbs-parallel-est-gate.sh | tee enumbs_result/kyber3-enumbs-parallel-est-gate.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 1000; // maximal selected blocksize value
    // params->J = 15; //maximal jump value
    params->threads = 50;
    params->enumbs_prec = 1e2;

     // kyber-III(Kyber-1024) round-3 parameters
    printf("============= Kyber-III\n");
    int dim =  1918;
    FP_NR<FT> dvol = 7242.6115232;
    gsa_est(dim, dvol, params);

    return 1;

}


