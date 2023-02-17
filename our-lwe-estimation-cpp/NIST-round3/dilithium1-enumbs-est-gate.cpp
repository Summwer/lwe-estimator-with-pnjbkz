#include "../framework/est.h"

//./dilithium1-enumbs-parallel-est-gate.sh | tee enumbs_result/dilithium1-enumbs-parallel-est-gate.log


int main(){

    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 500; // maximal selected blocksize value
    // params->J = 15; //maximal jump value
    

    printf("============= Dilithium-I\n");
    int dim =  2049;
    FP_NR<FT> dvol = 15614.219317244602;
    gsa_est(dim, dvol, params);


    return 1;

}