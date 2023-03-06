
#include "../framework/est.h"

//./NIST-round3-enumbs-parallel-est-gate.sh | tee enumbs_result/NIST-round3-enumbs-parallel-est-gate.log


int main(){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->max_dim = 500; // maximal selected blocksize value
    params->threads = 50;
    params->debug = true;
    params->enumbs_prec = 1e3;

    int dim;
    FP_NR<FT> dvol;

    // Kyber-I(Kyber-512) round-3 parameters
    printf("============= Kyber-I\n");
    //eta1 = eta2 = 3
    // dim = 1025;
    // dvol = 3944.9406103;

    //eta1 = 3, eta2 = 2
    dim = 1004;
    dvol = 3882.6780896;
    // gsa_est(dim, dvol, params);
    

    // Kyber-II(Kyber-768) round-3 parameters
    params->max_dim = 700; 
    printf("============= Kyber-II\n");
    dim =  1467;
    dvol =  5661.0782118;
    // gsa_est(dim, dvol, params);



    // kyber-III(Kyber-1024) round-3 parameters
    params->max_dim = 1000; 
    printf("============= Kyber-III\n");
    dim =  1918;
    dvol = 7242.6115232;
    // gsa_est(dim, dvol, params);


    // Dilithium-I round-3 parameters
    params->max_dim = 1000; 
    printf("============= Dilithium-I\n");
    dim =  2049;
    dvol = 15614.219317244602;
    gsa_est(dim, dvol, params);

    // Dilithium-II round-3 parameters
    params->max_dim = 700; 
    printf("============= Dilithium-II\n");
    dim =  2817;
    dvol = 21814.858106487554;
    gsa_est(dim, dvol, params);


    // Dilithium-III round-3 parameters
    params->max_dim = 1000; 
    printf("============= Dilithium-III\n");
    dim =  3841;
    dvol = 31317.16147360077;
    gsa_est(dim, dvol, params);


    return 1;



}


