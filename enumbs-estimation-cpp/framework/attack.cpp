#include "attack.h"


void call_enumbs(vector<double> l, Params* params){
    EnumBS* enumbs = new EnumBS(params);

    
    if(params->threads == 1){
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS)"<<endl;
    }else if (params->threads > 1){
        cout<<" Attack Estimation via simulation + probabilistic model (EnumBS in parallel)"<<endl;
    }
    else{
        cerr<<"Set bad threads = "<<params->threads<<"."<<endl;
    }
    printf("beta_start= %d, gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, G_prec = %e,  slope_prec = %e,  progressive_sieve = %d, ", params->beta_start, params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_G_prec, params->enumbs_slope_prec, params->progressive_sieve);
    if(params->worst_case)
        printf("worst_case, ");
    else
        printf("average_case, ");
    
    if(params->enum_add_G2)
        printf("Min G2 Strategy. \n");
    else
        printf("Min G Strategy. \n");
    auto start = system_clock::now();
    if(params->threads == 1)
        enumbs->enumbs_est(l);
    else
        enumbs->enumbs_est_in_parallel(l);
    auto finish = system_clock::now();
    duration<double> diff = finish - start;
    cout<<"EnumBS cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
    
}

// void call_search_tree(vector<double> l, Params* params){
//     SearchTree* search_tree = new SearchTree(params);

//     if(params->threads == 1){
//         cout<<" Attack Estimation via simulation + probabilistic model (Blocksize Strategy Search Tree)"<<endl;
//     }else if (params->threads > 1){
//         cout<<" Attack Estimation via simulation + probabilistic model (Blocksize Strategy Search Tree in parallel)"<<endl;
//     }
//     else{
//         cerr<<"Set bad threads = "<<params->threads<<"."<<endl;
//     }
//     // printf("beta_start= %d, gap = %d, J = %d, J_gap = %d, cost_model = %d, max_loop = %d, threads = %d, G_prec = %e,  slope_prec = %e,  progressive_sieve = %d, ", params->beta_start, params->gap, params->J, params->J_gap, params->cost_model, params->max_loop, params->threads, params->enumbs_G_prec, params->enumbs_slope_prec, params->progressive_sieve);
//     if(params->worst_case)
//         printf("worst_case, ");
//     else
//         printf("average_case, ");
    
//     // if(params->enum_add_G2)
//     //     printf("Min G2 Strategy. \n");
//     // else
//     //     printf("Min G Strategy. \n");
//     auto start = system_clock::now();
//     if(params->threads == 1)
//         search_tree->search_tree_est(l);
//     else
//         search_tree->search_tree_est_in_parallel(l);
//     auto finish = system_clock::now();
//     duration<double> diff = finish - start;
//     cout<<"Blocksize Strategy Search Tree cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
// }

void call_bssa(vector<double> l, Params* params, int sbeta, int gbeta){
    BSSA* bssa = new BSSA(params);
    
    auto start = system_clock::now();
    cout<<" Attack Estimation via simulation + probabilistic model (BSSA)"<<endl;
    printf("gap = %d, J = %d, J_gap = %d, max_loop = %d, cost_model = %d, mul_node = %d, progressive_sieve = %d\n", params->gap, params->J, params->J_gap, params->max_loop, params->cost_model, params->mul_node, params->progressive_sieve);
     if(params->worst_case)
        printf("worst_case. \n");
    else
        printf("average_case. \n");
    if(params->mul_node)
        bssa->bssa_est_mul_node(l, sbeta, gbeta);
    else
        bssa->bssa_est(l, sbeta, gbeta);
    auto finish = system_clock::now();
    duration<double> diff = finish - start;
    cout<<"BSSA cost:"<< setprecision(2)<<diff.count()<<"s."<<endl;
    
}


