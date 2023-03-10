#include "attack.h"


void call_enumbs(vector<double> l, Params* params){
    EnumBS* enumbs = new EnumBS(params);

    
    if(params->threads == 1){
        auto start = system_clock::now();
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS)"<<endl;
        printf("gap = %d, J = %d, J_gap = %d, max_loop = %d, cost_model = %d, prec = %e, strategy_add_method = %d\n", params->gap, params->J, params->J_gap, params->max_loop, params->cost_model, params->enumbs_prec, params->enumbs_add_strategy);
        enumbs->enumbs_est(l);
        auto finish = system_clock::now();
        duration<double> diff = finish - start;
        cout<<"EnumBS cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    }else if (params->threads > 1){
        auto start = system_clock::now();
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS in parallel)"<<endl;
        printf("gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, prec = %e, strategy_add_method = %d\n", params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_prec, params->enumbs_add_strategy);
        enumbs->enumbs_est_in_parallel(l);
        auto finish = system_clock::now();
        duration<double> diff = finish - start;
        cout<<"EnumBS in parallel cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    }
    else{
        cerr<<"Set bas threads = "<<params->threads<<"."<<endl;
    }
}



void call_bssa(vector<double> l, Params* params, int sbeta, int gbeta){
    BSSA* bssa = new BSSA(params);
    
    auto start = system_clock::now();
    cout<<" Attack Estimation via simulation + probabilistic model (BSSA)"<<endl;
    printf("gap = %d, J = %d, J_gap = %d, max_loop = %d, cost_model = %d\n", params->gap, params->J, params->J_gap, params->max_loop, params->cost_model);
    bssa->bssa_est(l, sbeta, gbeta);
    auto finish = system_clock::now();
    duration<double> diff = finish - start;
    cout<<"BSSA cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
}


