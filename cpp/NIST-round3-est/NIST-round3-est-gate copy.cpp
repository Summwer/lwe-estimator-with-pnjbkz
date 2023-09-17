
#include "../framework/est.h"





// params input in main function
// argv[0]: implemented file name
// argv[1]: method -- 1:enumba; 2: bssa
int main(int argc,char **argv){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->threads = 4;
    params->cost_model = 1; //sec model;
    // params->progressive_sieve = true;
    params->verbose = false;
    params->debug = false;
    params->worst_case = false;
    params->method = atoi(argv[1]); //1:enumbs;2:bssa
    params->gap = 1;
    params->J = 1; 
    params->J_gap = 1;
    params->enumbs_G_prec = 1./atoi(argv[2]);
    params->max_loop = 1;
    params->max_dim = 1500; 
    // params->bssa_tradion = true;  

    //dim and dvol value is selected from leaky-lwe-estimator.
    int dim;    
    FP_NR<FT> dvol;

    
    
    // Kyber-I(Kyber-512) round-3 parameters
    printf("============= Kyber-I\n");
    //eta1 = 3, eta2 = 2
    dim = 1004;
    dvol = 3882.6780896;
    gsa_est(dim, dvol, params);

    
    // Kyber-II(Kyber-768) round-3 parameters
    printf("============= Kyber-II\n");
    dim =  1467;
    dvol =  5661.0782118;
    gsa_est(dim, dvol, params);



    // kyber-III(Kyber-1024) round-3 parameters
    printf("============= Kyber-III\n");
    dim =  1918;
    dvol = 7242.6115232;
    gsa_est(dim, dvol, params);


    // Dilithium-I round-3 parameters 
    printf("============= Dilithium-I\n");
    dim =  2049;
    dvol = 15614.219317244602;
    gsa_est(dim, dvol, params);



    // Dilithium-II round-3 parameters
    printf("============= Dilithium-II\n");
    // dim =  2817;
    // dvol = 21814.858106487554;
    dim = 2654;
    dvol = 19371.0238433;
    gsa_est(dim, dvol, params);
    

    // Dilithium-III round-3 parameters
    printf("============= Dilithium-III\n");
    //dim =  3841;
    //dvol = 31317.16147360077;

    dim = 3540;
    dvol = 26623.1162463;
    gsa_est(dim, dvol, params);


    // // Frodo-I round-3 parameters 
    // printf("============= Frodo-I\n");
    // dim =  1297;
    // dvol = 5479.4593497;
    // gsa_est(dim, dvol, params);



    // // Frodo-II round-3 parameters
    // printf("============= Frodo-II\n");
    // dim = 1969;
    // dvol = 9347.2957371;
    // gsa_est(dim, dvol, params);
    

    // // Frodo-III round-3 parameters
    // printf("============= Frodo-III\n");
    // dim = 2634;
    // dvol = 13355.2889193;
    // gsa_est(dim, dvol, params);

    return 1;
}


