#include "bssa.h"


void BSSA::print_strategy(vector<BSSA::strategy> S){
    cout<<"Block size strategy generated by BSSA (beta,jump,tours):[";
    for(int i =0; i < int(S.size()); i ++){
        printf("(%3d,%3d,%3d)",S[i].beta,S[i].jump,S[i].tours);
        if(i!=int(S.size()) - 1)
            printf(",");
    }
    cout<<"]"<<endl;
}


void BSSA::print_BS(map<int,BSSA::blocksize_strategy> BS){
    cout<<"BS = ";
    cout<<"{";
    for (map<int, BSSA::blocksize_strategy> ::iterator it = BS.begin(); it != BS.end(); it++) {
        cout<<it-> first<<" : ";
        print_strategy(it->second.S);
        printf("(cum_avg_GB_BKZ = %3.2f gate, B_BKZ = %3.2f bit cum-pr = %3.2f)\n", it->second.cum_avg_GB_BKZ.first, it->second.cum_avg_GB_BKZ.second, it->second.cum_pr); 
	}
    cout<<"}"<<endl;
}

void BSSA::print_bs(blocksize_strategy bs){
    tuple<int,int,double,double,double> dsvp_t_ =  dsvp_predict(bs.l, bs.cum_pr, cost, params->cost_model,  bs.cum_GB_BKZ);
    int d = int(bs.l.size());
    
    cout<<"bs = ";

    if(params->cost_model == 1)
        printf("(slope = %3.5f, G_BKZ = %3.2f log2(gate), G_avg_BKZ = %3.2f log2(gate), B_BKZ = %3.2f log2(bit), cum-pr = %3.2f,  pump-{%d,%d,%d}, G_dsvp = %3.2f log2(gate), G_avg_dsvp = %3.2f log2(gate), B_dsvp = %3.2f log2(bit), G = %3.2f log2(gate), B = %3.2f log2(bit))\n", get_current_slope(bs.l,0,int(bs.l.size())), bs.cum_GB_BKZ.first, bs.cum_avg_GB_BKZ.first, bs.cum_avg_GB_BKZ.second, bs.cum_pr, d - get<1>(dsvp_t_), get<1>(dsvp_t_), get<1>(dsvp_t_) -  get<0>(dsvp_t_), get<4>(dsvp_t_), get<2>(dsvp_t_),get<3>(dsvp_t_), log2(pow(2,get<2>(dsvp_t_)) + pow(2,bs.cum_avg_GB_BKZ.first)), max(bs.cum_avg_GB_BKZ.second,get<3>(dsvp_t_)) );  
    if(params->cost_model == 2)
       printf("(slope = %3.5f, G_BKZ = %3.2f log2(sec), G_avg_BKZ = %3.2f log2(sec), B_BKZ = %3.2f log2(bit), cum-pr = %3.2f,  pump-{%d,%d,%d}, G_dsvp = %3.2f log2(sec), G_avg_dsvp = %3.2f log2(sec), B_dsvp = %3.2f log2(bit), G = %3.2f log2(sec), B = %3.2f log2(bit))\n", get_current_slope(bs.l,0,int(bs.l.size())), bs.cum_GB_BKZ.first, bs.cum_avg_GB_BKZ.first, bs.cum_avg_GB_BKZ.second, bs.cum_pr, d - get<1>(dsvp_t_), get<1>(dsvp_t_), get<1>(dsvp_t_) -  get<0>(dsvp_t_),get<4>(dsvp_t_), get<2>(dsvp_t_),get<3>(dsvp_t_), log2(pow(2,get<2>(dsvp_t_)) + pow(2,bs.cum_avg_GB_BKZ.first)), max(bs.cum_avg_GB_BKZ.second,get<3>(dsvp_t_)) ); 

    print_strategy(bs.S);
}

bool BSSA::pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB, pair<double,double> &cum_avg_GB,  double &cum_pr, int beta, int jump, tuple<int,int,double,double,double> &dsvp_t_, double &slope){

    //simulate pnj-bkz more precisely
    int d = l.size(),beta_ = get_beta_(params,beta,jump,d);


    double rem_pr = 1. - cum_pr;

    vector<double> l_;
    
    sim -> simulate(l_,l,beta,jump,1);

    if(l_[d-beta_] == l[d-beta_]){
        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model, cum_GB);
        slope = get_current_slope(l, 0, d);
        return false;
    }
    else{
        l = l_;
        slope = get_current_slope(l, 0, d);
        boost::math::chi_squared chisquare(beta_);
        double pr = boost::math::cdf(chisquare,pow(2,2.*l[d-beta_]));

        pair<double,double> GB = cost->bkz_cost(d,beta,jump,params->cost_model);
        cum_GB.first = log2(pow(2,cum_GB.first)+pow(2,GB.first));
        cum_GB.second = GB.second;
        // if(not params->worst_case){
        cum_avg_GB.first = log2(pow(2,cum_avg_GB.first)+(pow(2,cum_GB.first)*rem_pr*pr));
        cum_avg_GB.second = log2(pow(2,cum_avg_GB.second)+(pow(2,cum_GB.second)*rem_pr*pr));


        // }
        // else{
        //     //cum_avg_GB = cum_GB;
        //     cum_avg_GB.first = log2(pow(2,cum_avg_GB.first)+(pow(2,cum_GB.first)*(pr-pre_pr)));
        //     cum_avg_GB.second = log2(pow(2,cum_avg_GB.second)+pow(2,cum_GB.second) * (pr-pre_pr));
        // }

        cum_pr += rem_pr * pr;
        rem_pr = 1. - cum_pr;

        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model,  cum_GB);

        return true;
    }
}


pair<double,double> BSSA::max_tour_for_pnjbkz_beta(vector<double> l, int beta){
    
    tuple<int,int,double,double,double> dsvp_t1;
    double G21;
   
    // vector<strategy> S = bs.S;
    double cum_pr = 0.;
    pair<double,double> cum_GB= make_pair(0.,0.), cum_avg_GB= make_pair(0.,0.);
    int  loop = 0;
    
    double slope0  = get_current_slope(l, 0, l.size()), slope1;
    bool sim_term = pnjbkz_beta_loop(l, cum_GB, cum_avg_GB, cum_pr, beta, 1,dsvp_t1, slope1);

    G21 = get<4>(dsvp_t1);

    if(not sim_term)
        return make_pair(G21, slope1);
    
    // assert(G21 >= 0.);

    // cout<<slope1<<", "<<slope0<<endl;
    
    while( (slope1 - slope0) > 0 && loop < params->max_loop){

        loop +=1;
        slope0 = slope1;
        if(cum_pr >= params->succ_prob or not sim_term)
            break;
    
        sim_term = pnjbkz_beta_loop( l, cum_GB, cum_avg_GB, cum_pr, beta, 1, dsvp_t1, slope1);
        G21 = get<4>(dsvp_t1);
        // assert(G21 >= 0.);
    }
    // assert(loop>=1);
    // if(loop == 0)
    //     loop = 1;
    // return dsvp_t1;
    return make_pair(G21, slope1);
}



pair<double,double> BSSA::max_tour_for_pnjbkz_beta(blocksize_strategy bs, int beta){
    
    tuple<int,int,double,double,double> dsvp_t1;
    double G21;
    vector<double> l = bs.l;
   
    // vector<strategy> S = bs.S;
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB= bs.cum_GB_BKZ, cum_avg_GB= bs.cum_avg_GB_BKZ;
    int  loop = 0;
    
    double slope0  = get_current_slope(l, 0, l.size()), slope1;
    bool sim_term = pnjbkz_beta_loop(l, cum_GB, cum_avg_GB, cum_pr, beta, 1,dsvp_t1, slope1);

    G21 = get<4>(dsvp_t1);

    if(not sim_term)
        return make_pair(G21, slope1);
    
    // assert(G21 >= 0.);

    // cout<<slope1<<", "<<slope0<<endl;
    
    while( (slope1 - slope0) > 0 && loop < params->max_loop){

        loop +=1;
        slope0 = slope1;
        if(cum_pr >= params->succ_prob or not sim_term)
            break;
    
        sim_term = pnjbkz_beta_loop( l, cum_GB, cum_avg_GB, cum_pr, beta, 1, dsvp_t1, slope1);
        G21 = get<4>(dsvp_t1);
        // assert(G21 >= 0.);
    }
    // assert(loop>=1);
    // if(loop == 0)
    //     loop = 1;
    // return dsvp_t1;
    return make_pair(G21, slope1);
}






BSSA::blocksize_strategy BSSA::min_tour_to_each_goal_beta(BSSA::blocksize_strategy bs, int beta, int jump, pair<double,double> goal_quality){

    tuple<int,int,double,double,double> dsvp_t1;
    double G20 =  params->max_num, G21, goal_G2 = goal_quality.first, goal_slope = goal_quality.second;

    // int dsvp0 = params->max_num;
   
    vector<double> l = bs.l; 
    
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB= bs.cum_GB_BKZ, cum_avg_GB= bs.cum_avg_GB_BKZ;

    int  loop = 0;
    
    
    double slope0 = get_current_slope(l, 0, l.size()), slope1;
    bool sim_term = pnjbkz_beta_loop(l, cum_GB, cum_avg_GB, cum_pr, beta, jump, dsvp_t1, slope1);
    
    if(not sim_term)
        return {};
    
    G21 = get<4>(dsvp_t1);
    // if(G21<0)
    //     cout<<"G21: "<<G21<<endl;
    // assert(G21 >= 0.);
    
    while( (slope1 - slope0) > params->enumbs_slope_prec && loop < params->max_loop){
        loop +=1;
        G20 = G21;
        slope0 = slope1;
        if(G20 <= goal_G2 and slope0 >= goal_slope)
            break;
    
        sim_term = pnjbkz_beta_loop( l, cum_GB, cum_avg_GB, cum_pr, beta, jump, dsvp_t1, slope1);
        G21 = get<4>(dsvp_t1);
        // assert(G21 >= 0.);

    }

    if(loop>=1 and G20 <= goal_G2 and slope0 >= goal_slope ){
        vector<BSSA::strategy> S = bs.S;
        S.insert(S.end(),{beta, jump, loop});
        return {S, l, cum_GB, cum_avg_GB, cum_pr};
    }
    else
        return {};
}



//beta mode, the node is depend on cost of bkz-betas.
void BSSA::bssa_est(vector<double> l0, int sbeta, int gbeta){
    /*
    param: node_start, node_goal, node_mid, node
    */
    int d= l0.size(), len_S = 0;
    pair<double,double> goal_quality;

    BS.insert(pair<int,BSSA::blocksize_strategy>(sbeta,{{}, l0, make_pair(0.,0.), make_pair(0.,0.), 0.}));

    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

    for(int beta = sbeta; beta <= gbeta; beta+=params->gap ){
        double G_BKZ_min = params->max_num;

        if(params->bssa_tradion)
            goal_quality = max_tour_for_pnjbkz_beta(l0, beta);//1: G2, 2: slope

        for(int ssbeta = sbeta; ssbeta < beta; ssbeta+=params->gap){
            BSSA::blocksize_strategy bs = BS[ssbeta], bs_min = {}, bs_tmp = {};
            if(bs.cum_pr >= params->succ_prob){
                // cout<<"ssbeta = "<<ssbeta<<endl;
                beta = ssbeta;
                gbeta = ssbeta;
                break;
            }
                
            double G_tmp_min = params->max_num; 

            if(not params->bssa_tradion)
                goal_quality = max_tour_for_pnjbkz_beta(bs, beta);//1: G2, 2: slope
            
            bool flag = true;
            
            for(int beta_alg = beta+1; beta_alg < min(d,params->max_dim); beta_alg++){
                // int len_S = bs_tmp.S.size();
                // if ( len_S != 0 && beta_alg >= bs_min.S[len_S-1].beta + 20)
                //         break;
                for(int j = params->J; j > 0; j-=params->J_gap){
                    if(params->verbose)
                        printf("\r Blocksize strategy selection process: %4d --> %4d --> (%4d, %4d) --> %4d --> %4d", sbeta,ssbeta,beta_alg,j,beta,gbeta);
                    pair<double, double> GB_alg = cost->bkz_cost(d,beta_alg,j,params->cost_model);

                    // int f = get_f_for_pnjbkz(params,beta_alg);
                    int f = get_f_for_pnjbkz(params,beta_alg);

                    //((f == 0 or beta_alg < 79 )&& j > 1) ||
                    if(((f == 0 or beta_alg < 79 )&& j > 1) || ( f!=0 && j > floor((double) f/2.))){ 
                        continue;
                    }

                    if(beta_alg < 55 and j > 1)
                        continue;

                    // if(params->cost_model == 1 and j >= 0.1 * beta)
                    //     continue;
                    // else{
                    //     if(f==0 && j>1)
                    //         cout<<"f="<<f<<", j="<<j <<"beta_alg="<<beta_alg<<endl;
                    // //     if((f == 0 or beta_alg < 79 )&& j > 1)
                    // //         cout<<"f="<<f<<", j="<<j<<endl;
                    // }

                    //if(((f == 0 or beta_alg < 79 )&& j > 1) || (f!=0 && j >= min(f,ceil(0.1* beta_alg))))
                       

                    // if(GB_alg.first > G_tmp_min  and (j == params->J or j < min((double) f,ceil(0.1* beta_alg)))){
                    //If then, we don't need to check the larger beta_alg
                    if(GB_alg.first > G_tmp_min  and (j == params->J or j <= floor((double) f/2.))){
                        flag = false;
                        break;
                    }


                    // if(GB_alg.first > G_tmp_min){
                    //     cout<<"beta_alg = "<< beta_alg << ", j = "<<j<<endl;
                    //     flag = false;
                    //     break;
                    // }
                    
                    bs_tmp = min_tour_to_each_goal_beta(bs, beta_alg, j, goal_quality);
                    len_S = bs_tmp.S.size();

                    // printf("======================\n");
                    // cout<<"goal_G2 = "<<get<2>(dsvp_t0)<<", beta = "<<beta<<endl;
                    // print_bs(bs);
                    // 
                    // printf("----------------------\n");

                    if(len_S > 0 && G_tmp_min > bs_tmp.cum_GB_BKZ.first){
                        G_tmp_min = bs_tmp.cum_GB_BKZ.first;
                        bs_min = bs_tmp;
                        continue;
                    }

                    // if(len_S > 0 and bs_tmp.S[len_S-1].tours == 1 and (j == params->J or j == min((double) f,ceil(0.1* beta_alg)))){
                        
                    //If then, we don't need to check the larger beta_alg
                    if(len_S > 0 and bs_tmp.S[len_S-1].tours == 1 and (j == params->J or j == floor((double) f/2.))){
                        flag = false;
                        break;
                    }
                }
                if(not flag)
                    break;
            }
            
            if(bs_min.S.size() !=0 && G_BKZ_min > G_tmp_min){
                G_BKZ_min = G_tmp_min;
                if(BS.find(beta)==BS.end())
                    BS.insert(pair<int,BSSA::blocksize_strategy>(beta,bs_min));
                else
                    BS[beta] = bs_min;
                // cout<<endl;
                // printf("beta = %d, goal_G2 = %3.7f, goal_slope = %3.7f", beta, goal_quality.first, goal_quality.second);
                // cout<<"beta = "<<beta<< endl;
                // print_bs(bs_min);
            }
        }
    }
    if(params->verbose)
        cout<<endl;


    //Find the optimized strategy
    double Gmin = params->max_num, Bmin  = params->max_num, G1, G2, G,  B;
    BSSA::blocksize_strategy bsmin;
    for (map<int, BSSA::blocksize_strategy> ::iterator it = BS.begin(); it != BS.end(); it++) {
        tuple<int,int,double,double,double> dsvp_t0 = dsvp_predict(it->second.l, it->second.cum_pr, cost,params->cost_model, it->second.cum_GB_BKZ);
        G1 = it->second.cum_avg_GB_BKZ.first;
        G2 = get<2>(dsvp_t0);
        
        // G1 = it->second.cum_GB_BKZ.first;
        // G2 = get<4>(dsvp_t0);
        G = log2(pow(2,G1)+pow(2,G2));
        B = log2(pow(2,get<3>(dsvp_t0))+pow(2,it->second.cum_avg_GB_BKZ.second));
        if(G < Gmin and B < params->max_RAM){
            bsmin = it->second;
            Gmin = G;
            Bmin = B;
        }
    }
    printf("Find the optimized Strategy through BSSA!!\n");
    print_bs(bsmin);
    
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec) = %3.2f s, Memory Cost = %3.2f log2(bit) = %3.2f GB \n", Gmin, pow(2,Gmin), Bmin, pow(2,Bmin-33));

}


// pair<double,double> BSSA::strategy_verification(vector<double> l,vector<strategy> S){

//     int d = l.size();
//     double cum_pr = 0., rem_pr = 1., proba, G1cum=0., B1cum = 0.;
//     // BKZJSim* sim = new BKZJSim(params);
//     // COST* cost = new COST();
//     for(int i = 0; i< int(S.size()); i++){
//         BSSA::strategy bs = S[i];
//         int beta = bs.beta, jump = bs.jump, N = bs.tours;
//         for(int tour = 0; tour < N; tour++){
        
//             sim -> simulate(l,l,beta,jump,1);

//             boost::math::chi_squared chisquare(beta);
//             proba = boost::math::cdf(chisquare,pow(2,2.*l[d-beta]));
            

//             pair<double,double> GB = cost -> bkz_cost(d,beta,jump,params->cost_model);
//             if(not params->worst_case)
//                 G1cum = log2(pow(2,G1cum) + (pow(2,GB.first) * rem_pr * proba));
//             else
//                 G1cum = log2(pow(2,G1cum) + pow(2,GB.first));

//             B1cum = max(B1cum,GB.second);

//             cum_pr += rem_pr * proba;
//             rem_pr *= 1. - proba;
//         }
//     }

//     tuple<int,int,double,double,double> dsvp_t =  dsvp_predict(l, cum_pr, cost,params->cost_model,  params->worst_case);
//     double G2 = get<2>(dsvp_t);
//     double G = log2(pow(2,G1cum)+pow(2,G2));
//     printf("Verified cum_pr = %e \n ", cum_pr);
//     printf("Verified G1 = %3.2f, G2 = %3.2f, dsvp = %3.2f\n", G1cum,G2,get<0>(dsvp_t));
//     printf("G = %3.2f\n", G );

//     return make_pair(G1cum, cum_pr);

// }
