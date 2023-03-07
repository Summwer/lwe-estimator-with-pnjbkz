#include "bssa.h"


void BSSA::print_strategy(vector<BSSA::strategy> S){
    cout<<"Block size strategy generated by BSSA (beta,jump,tours):{";
    for(int i =0; i < int(S.size()); i ++){
        printf("(%3d,%3d,%3d)",S[i].beta,S[i].jump,S[i].tours);
        if(i!=int(S.size()) - 1)
            printf(",");
    }
    cout<<"}"<<endl;
}


void BSSA::print_BS(map<int,BSSA::blocksize_strategy> BS){
    cout<<"BS = ";
    cout<<"{";
    for (map<int, BSSA::blocksize_strategy> ::iterator it = BS.begin(); it != BS.end(); it++) {
        cout<<it-> first<<" : ";
        print_strategy(it->second.S);
        printf("(GB_BKZ = %3.2f gate, B_BKZ = %3.2f bit cum-pr = %3.2f)\n", it->second.GB_BKZ.first, it->second.GB_BKZ.second, it->second.cum_pr); 
	}
    cout<<"}"<<endl;
}

void BSSA::print_bs(blocksize_strategy bs){
    tuple<double,int,double,double> dsvp_t_ =  dsvp_predict(bs.l, bs.cum_pr, cost, params->cost_model, params->progressive_sieve);
    cout<<"bs = ";

    printf("(G_BKZ = %3.2f gate, B_BKZ = %3.2f bit cum-pr = %3.2f, dsvp = %3.2f, dsvp_r = %3d, G_dsvp = %3.2f gate, B_dsvp = %3.2f bit, G = %3.2f gate, B = %3.2f bit)\n",  bs.GB_BKZ.first, bs.GB_BKZ.second, bs.cum_pr, get<0>(dsvp_t_), get<1>(dsvp_t_),get<2>(dsvp_t_),get<3>(dsvp_t_), log2(pow(2,get<2>(dsvp_t_)) + pow(2,bs.GB_BKZ.first)), max(bs.GB_BKZ.second,get<3>(dsvp_t_)) );   

    print_strategy(bs.S);
}

bool BSSA::pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB, double &cum_pr, int beta, int jump, tuple<double,int,double,double> &dsvp_t_){

    //simulate pnj-bkz more precisely
    int f, beta_, d = l.size();
    f = default_dim4free_fun(beta);
    if(jump <= 2)
        beta_ = beta;
    else if(jump >=3 and jump <=4)
        beta_ = get_beta_from_sieve_dim(beta-f,d,2);
    else if(jump>=5)
        beta_ = get_beta_from_sieve_dim(beta-f,d,1);

    double rem_pr = 1. - cum_pr;

    vector<double> l_;
    
    sim -> simulate(l_,l,beta_,jump,1);

    if(l_[d-beta_] == l[d-beta_])
        return false;
    else{
        l = l_;
        boost::math::chi_squared chisquare(beta_);
        double pr = boost::math::cdf(chisquare,pow(2,2.*l[d-beta_]));

        pair<double,double> GB = cost->bkz_cost(d,beta,jump,params->cost_model);

        
        cum_GB.first = log2(pow(2,cum_GB.first)+(pow(2,GB.first)*rem_pr*pr));
        cum_GB.second = max(cum_GB.second, GB.second);


        cum_pr += rem_pr * pr;
        rem_pr = 1. - cum_pr;


        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model, params->progressive_sieve);

        return true;

    }

}


tuple<double,int,double,double> BSSA::max_tour_for_pnjbkz_beta(BSSA::blocksize_strategy bs, int beta){
    
    tuple<double,int,double,double> dsvp_t1;
    double G20 =  MAX_NUM, G21;
   
    vector<double> l = bs.l; 
    // vector<strategy> S = bs.S;
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB=bs.GB_BKZ;

    int  loop = 0;
    
    bool sim_term = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, 1,dsvp_t1);
    
    G21 = get<2>(dsvp_t1);

    assert(G21 >= 0.);
    //
    while((abs(G20 - MAX_NUM)<0.01 or G20 > G21)  and loop < params->max_loop){

        loop +=1;
        G20 = G21;

        if(cum_pr >= 0.999 or not sim_term)
            break;
    
        sim_term = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, 1, dsvp_t1);
        G21 = get<2>(dsvp_t1);
        assert(G21 >= 0.);
    }
    assert(loop>=1);
    
    return dsvp_t1;
}


BSSA::blocksize_strategy BSSA::min_tour_to_each_goal_beta(BSSA::blocksize_strategy bs, int beta, int jump, double G2_star){

    tuple<double,int,double,double> dsvp_t1;
    double G20 =  MAX_NUM, G21;
   
    vector<double> l = bs.l; 
    
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB=bs.GB_BKZ;

    int  loop = 0;
    
    
    bool sim_term = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, jump, dsvp_t1);
    
    G21 = get<2>(dsvp_t1);

    assert(G21 >= 0.);
    
    while((abs(G20 - MAX_NUM)<0.01 or G20 > G21) and loop < params->max_loop){
        loop +=1;
        G20 = G21;

        if(G20 <= G2_star or cum_pr >= 0.999 or not sim_term)
            break;
    
        sim_term = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, jump, dsvp_t1);
        G21 = get<2>(dsvp_t1);
        assert(G21 >= 0.);

    }

    if(loop>=1 and G20 <= G2_star){
        vector<BSSA::strategy> S = bs.S;
        S.insert(S.end(),{beta, jump, loop});
        return {S, l, cum_GB, cum_pr};
    }
    else
        return {};
}



// //node mode, we can set nodes more freely.
// // void BSSA::bssa_generate(vector<double> l0, int sbeta, int gbeta, int jump, int gap, int J_gap){
// //     /*
// //     param: node_start, node_goal, node_mid, node
// //     */
// //     int d= l0.size();
// //     int node, node_start, node_goal = 0, node_mid, beta;
// //     int dsvp_star = ceil(dsvp_predict(l0, 0, 0.));
// //     node_start = -dsvp_star; //quality increase while value of dsvp decreases.

// //     double G_min,G_tmp;

// //     for(node_mid = node_start; node_mid <= node_goal; node_mid+= max(gap,1)){
// //         G_min = MAX_NUM;
// //         for(node = node_start; node < node_mid; node += max(gap,1)){
// //             printf("\r Blocksize strategy selection process: %4d --> %4d  --> %4d --> %4d", node_start,node,node_mid,node_goal);
// //             if(node == node_start){
// //                 bs0 = {{},l0,0.,0.};
// //             }
// //             else if(node > node_start){
// //                 if(BS.find(node)==BS.end()){
// //                     bssa_generate(l0, node_start, node);
// //                 }else{
// //                     bs0 = BS[node];
// //                 }
// //             }

// //             G_tmp = MAX_NUM;
// //             //node_tmp = MAX_DIM;
// //             //beta = max(bs0.S);
// //             print_strategy(bs0.S);
// //             // dsvp_star = d - ssbeta
// //             beta = 50;
// //             for(int beta_alg = beta+1; beta_alg < min(MAX_DIM,d); beta_alg++){
// //                 printf("\r Blocksize strategy selection process: %4d --> %4d --> (%4d) --> %4d --> %4d", node_start,node,beta_alg,node_mid,node_goal);
// //                 if (beta_alg >= beta_tmp + 3)
// //                     break
// //                 for(int j=1; j < J+1; j+=J_gap){
// //                     if(G_tmp == 0)
// //                         break;
// //                     // l_, _, loop_, G_, cumulated_proba_= min_tour_to_each_goal_dsvp(l,cum_pr,beta,jump=j,dsvp_star = dsvp_star);
                        
// //                 //     if G_tmp > G_:
// //                 //         S_tmp, l_tmp, G_tmp, cumulated_proba_tmp,beta_tmp =[(beta_alg,j) for _ in range(loop_)], l_,  G_,cumulated_proba_, beta_alg
// //                 }
// //             }
// //         }
// //     }
// //     printf("\n");
// // }



//beta mode, the node is depend on cost of bkz-betas.
void BSSA::bssa_est(vector<double> l0, int sbeta, int gbeta){
    /*
    param: node_start, node_goal, node_mid, node
    */
    int d= l0.size(), len_S;

    BS.insert(pair<int,BSSA::blocksize_strategy>(sbeta,{{}, l0, make_pair(0.,0.), 0.}));

    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

    for(int beta = sbeta; beta <= gbeta; beta++ ){
        double G_BKZ_min = MAX_NUM;
        for(int ssbeta = sbeta; ssbeta < beta; ssbeta++){
            BSSA::blocksize_strategy bs = BS[ssbeta], bs_min = {}, bs_tmp = {};
            if(bs.cum_pr >= 0.999){
                // cout<<"ssbeta = "<<ssbeta<<endl;
                beta = ssbeta;
                gbeta = ssbeta;
                break;
            }
                
            // cout<<endl;
            // print_bs(bs);
            // len_S = bs.S.size();
            // if(len_S == 0){
            //     beta_start = sbeta;
            //     j_start = params->J;
            // }
            // else{
            //     j_start = bs.S[len_S-1].jump;
            //     beta_start = bs.S[len_S-1].beta;
            //     if(j_start==1){
            //         beta_start += 1;
            //         j_start = params->J;
            //     }
            //     else
            //         j_start -= params->J_gap;
            // }

            double G_tmp_min = MAX_NUM;

            tuple<double,int,double,double> dsvp_t0 = max_tour_for_pnjbkz_beta(bs, beta);

            
            for(int beta_alg = beta+1; beta_alg < d; beta_alg++){
                len_S = bs_tmp.S.size();
                // if ( len_S != 0 && beta_alg >= bs_min.S[len_S-1].beta + 3)
                //         break;
                for(int j = params->J; j > 0; j--){
                    if(params->verbose)
                        printf("\r Blocksize strategy selection process: %4d --> %4d --> (%4d, %2d) --> %4d --> %4d", sbeta,ssbeta,beta_alg,j,beta,gbeta);
                    pair<double, double> GB_alg = cost->bkz_cost(d,beta_alg,j,params->cost_model);

                    if(GB_alg.first > G_tmp_min)
                        break;

                    int f = dims4free(beta_alg);
                    if((f == 0 && j > 1) || (f!=0 && j >= f))
                        continue;
                    
                    bs_tmp = min_tour_to_each_goal_beta(bs, beta_alg, j, get<2>(dsvp_t0));
                    
                    // printf("======================\n");
                    // cout<<"G2_star = "<<get<2>(dsvp_t0)<<", beta = "<<beta<<endl;
                    // print_bs(bs);
                    // 
                    // printf("----------------------\n");

                    if(bs_tmp.S.size() != 0 && G_tmp_min > bs_tmp.GB_BKZ.first){
                        G_tmp_min = bs_tmp.GB_BKZ.first;
                        bs_min = bs_tmp;
                    }

                    
                }
            }
            
            if(bs_min.S.size() !=0 && G_BKZ_min > G_tmp_min){
                G_BKZ_min = G_tmp_min;
                if(BS.find(beta)==BS.end())
                    BS.insert(pair<int,BSSA::blocksize_strategy>(beta,bs_min));
                else
                    BS[beta] = bs_min;
                // cout<<endl;
                // cout<<"beta = "<<beta<< endl;
                // print_bs(bs_min);
            }
        }
    }
    if(params->verbose)
        cout<<endl;


    //Find the optimized strategy
    double Gmin = MAX_NUM, Bmin  = MAX_NUM, G1, G2, G,  B;
    BSSA::blocksize_strategy bsmin;
    for (map<int, BSSA::blocksize_strategy> ::iterator it = BS.begin(); it != BS.end(); it++) {
        tuple<double,int,double,double> dsvp_t0 = dsvp_predict(it->second.l, it->second.cum_pr, cost,params->cost_model, params->progressive_sieve);
        G1 = it->second.GB_BKZ.first;
        G2 = get<2>(dsvp_t0);
        G = log2(pow(2,G1)+pow(2,G2));
        B = max(get<3>(dsvp_t0),it->second.GB_BKZ.second);
        if(G<Gmin){
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
        printf("Min Cost = %3.2f log2(sec), Memory Cost = %3.2f log(bit)\n", Gmin, Bmin);

}


pair<double,double> BSSA::strategy_verification(vector<double> l,vector<strategy> S){

    int d = l.size();
    double cum_pr = 0., rem_pr = 1., proba, G1cum=0., B1cum = 0.;
    BKZJSim* sim = new BKZJSim();
    COST* cost = new COST();
    for(int i = 0; i< int(S.size()); i++){
        BSSA::strategy bs = S[i];
        int beta = bs.beta, jump = bs.jump, N = bs.tours;
        for(int tour = 0; tour < N; tour++){
        
            sim -> simulate(l,l,beta,jump,1);

            boost::math::chi_squared chisquare(beta);
            proba = boost::math::cdf(chisquare,pow(2,2.*l[d-beta]));
            

            pair<double,double> GB = cost -> bkz_cost(d,beta,jump,params->cost_model);
            G1cum = log2(pow(2,G1cum) + (pow(2,GB.first) * rem_pr * proba));
            B1cum = max(B1cum,GB.second);

            cum_pr += rem_pr * proba;
            rem_pr *= 1. - proba;
        }
    }

    tuple<double,int,double,double> dsvp_t =  dsvp_predict(l, cum_pr, cost,params->cost_model, params->progressive_sieve);
    double G2 = get<2>(dsvp_t);
    double G = log2(pow(2,G1cum)+pow(2,G2));
    // printf("Verified cum_pr = %e \n ", cum_pr);
    // printf("Verified G1 = %3.2f, G2 = %3.2f, dsvp = %3.2f\n", G1cum,G2,get<0>(dsvp_t));
    // printf("G = %3.2f\n", G );

    return make_pair(G1cum, cum_pr);

}
