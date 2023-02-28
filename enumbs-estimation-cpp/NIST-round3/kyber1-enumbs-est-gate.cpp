
#include "../framework/est.h"

//./kyber1-enumbs-est-gate.sh | tee enumbs_result/kyber1-enumbs-est-gate.log

int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 500; // maximal selected blocksize value
    // params->J = 15; //maximal jump value
    params->enumbs_prec = 1e2;

    // Kyber-I(Kyber-512) round-3 parameters
    printf("============= Kyber-I\n");
    //eta1 = eta2 = 3
    // int dim = 1025;
    // FP_NR<FT> dvol = 3944.9406103;

    //eta1 = 3, eta2 = 2
    int dim = 1004;
    FP_NR<FT> dvol = 3882.6780896;
    gsa_est(dim, dvol, params);
    return 1;

}


