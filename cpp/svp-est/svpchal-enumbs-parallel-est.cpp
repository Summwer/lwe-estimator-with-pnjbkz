#include "../framework/est.h"


//params input in main function
//argv[0]: implemented file name
//argv[1]: jump value J
//argv[2]: max_loop value
//argv[3]: cost model: 1.gate model 2.practical sec model
//argv[4]: maximal dimension in enumeration
//argv[5]: threads number
//argv[6]: practical_pump_d4f
//argv[7]: start beta value.
//argv[8]: est model in dsvp_prediction for last pump
int main(int argc,char **argv){
    Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    params->J = atoi(argv[1]); 
    params->max_loop = atoi(argv[2]);
    params->cost_model = atoi(argv[3]); //sec model;
    params->max_dim = atoi(argv[4]);
    params->threads = atoi(argv[5]);
    params->practical_pump_d4f = atoi(argv[6]);
    
    // params->enumbs_min_G = false;
    params->beta_start = atoi(argv[7]);
    // params->est_model = atoi(argv[8]);
    params->worst_case = false;
    params->enumbs_slope_prec = 1e-6;
    params->min_G_prec = 1e-3;
    params->enumbs_G_prec = 1e-3;
    if(params->cost_model == 1){
        if(atoi(argv[8]) == 1)
            params->list_decoding = "agps20"; //"matzov22"
        if(atoi(argv[8]) == 2)
            params->list_decoding = "matzov22"; //"matzov22"
        // params->print_Gcums = true;
    }
    if(atoi(argv[8]) == 1)
        params->max_RAM = 30; //1.5T = 43.58, 500G = 38.9

    params->dsvp_predict_param = 2;

    params->debug =  0;

    vector<int> svps;
    svps = {180};

    for(int i = 0; i < int(svps.size());i++){
        int n = svps[i];
        svpchal_est(n, params);
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
