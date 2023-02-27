#include "call_enumbs.h"


void call_enumbs(vector<double> l, Params* params){
    EnumBS* enumbs = new EnumBS(params);

    
    if(params->threads == 1){
        auto start = system_clock::now();
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS)"<<endl;
        printf("gap = %d, J = %d, J_gap = %d, max_loop = %d, cost_model = %d, prec = %e, enum_bound = %e\n", params->gap, params->J, params->J_gap, params->max_loop, params->cost_model, params->enumbs_prec, params->enumbs_bound);
        enumbs->enumbs_est(l);
        auto finish = system_clock::now();
        duration<double> diff = finish - start;
        cout<<"EnumBS cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    }else if (params->threads > 1){
        auto start = system_clock::now();
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS in parallel)"<<endl;
        printf("gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, prec = %e, enum_bound = %e\n", params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_prec, params -> enumbs_bound);
        enumbs->enumbs_est_in_parallel(l);
        auto finish = system_clock::now();
        duration<double> diff = finish - start;
        cout<<"EnumBS in parallel cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    }
    else{
        cerr<<"Set bas threads = "<<params->threads<<"."<<endl;
    }
}


