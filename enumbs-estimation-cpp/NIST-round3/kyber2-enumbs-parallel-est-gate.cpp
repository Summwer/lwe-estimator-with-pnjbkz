
#include "../framework/est.h"

//./kyber2-enumbs-parallel-est-gate.sh | tee enumbs_result/kyber2-enumbs-parallel-est-gate.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 700; // maximal selected blocksize value
    // params->J = 15; //maximal jump value
    params->threads = 50;
    // params->threads = 1;
    params->enumbs_prec = 1e2;


    // Kyber-II(Kyber-768) round-3 parameters
    printf("============= Kyber-II\n");
    int dim =  1467;
    FP_NR<FT> dvol =  5661.0782118;
    gsa_est(dim, dvol, params);

    return 1;

}


