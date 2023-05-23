#include "../framework/est.h"


//params input in main function
//argv[0]: implemented file name
//argv[1]: jump value J
//argv[2]: max_loop value
//argv[3]: cost model: 1.gate model 2.practical sec model
//argv[4]: maximal dimension in enumeration
//argv[5]: threads number
//argv[6]: practical_pump_d4f
int main(int argc,char **argv){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->J = atoi(argv[1]); 
    params->max_loop = atoi(argv[2]);
    params->cost_model = atoi(argv[3]); //sec model;
    params->max_dim = atoi(argv[4]);
    params->threads = atoi(argv[5]);
    params->practical_pump_d4f = atoi(argv[6]);
    // params->max_RAM = 43.58; //1.5T = 43.58
    // params->enumbs_min_G = false;

    vector<pair<int,double>> lwes;

    //low_dim_lwechallenge_est. 
    lwes= {{40, 0.025}, {45, 0.020}, {50, 0.015}, {55, 0.010}, {60, 0.010}, {70, 0.005}, {75, 0.005},{40,0.035},{40,0.045}, {45, 0.035}, {50, 0.030}, {55, 0.025}, {60, 0.020}, {65, 0.015}, {75, 0.010}, {95, 0.005}};
    for(int i = 0; i < int(lwes.size());i++){
        int n = lwes[i].first;
        double alpha  = lwes[i].second;
        lwechal_est(n, alpha, params);
    }

    
    // //solve_lwechal_est
    // lwes =  {{45,0.030}, {50, 0.025}, {55, 0.020}, {60, 0.015}, {85, 0.005}, {90, 0.005}};
    // for(int i = 0; i < int(lwes.size());i++){
    //     int n = lwes[i].first;
    //     double alpha  = lwes[i].second;
    //     lwechal_est(n, alpha, params);
    // }

 
    
    return 1;
}
