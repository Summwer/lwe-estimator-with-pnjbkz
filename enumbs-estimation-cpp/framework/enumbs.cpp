#include "enumbs.h"


void EnumBS::set_threads(int nr)
{
    assert(nr >= 1);
    threadpool.resize(nr);
}

void EnumBS::print_strategy(vector<EnumBS::strategy> S){
    cout<<"S(beta,jump,tours):{";
    for(int i =0; i < int(S.size()); i ++){
        printf("(%3d,%3d,%3d)",S[i].beta,S[i].jump,S[i].tours);
        if(i!=int(S.size()) - 1)
            printf(",");
    }
    cout<<"}"<<endl;
}

void EnumBS::print_BS(vector<blocksize_strategy> BS){
    cout<<"BS = ";
    double dsvp;
    int dsvp_r;
    cout<<"{";
    for(int i =0; i< int(BS.size()); i++){
        //print_strategy(BS[i].S);
        dsvp = get<0>(BS[i].dsvp_t), dsvp_r = get<1>(BS[i].dsvp_t);
        printf("(dsvp = %3.2f, dsvp_r = %3d, GB_BKZ = %3.2f gate, B_BKZ = %3.2f bit cum-pr = %3.2f)", dsvp, dsvp_r, BS[i].GB_BKZ.first, BS[i].GB_BKZ.second, BS[i].cum_pr);   
        if(i!=int(BS.size()) - 1)
            printf(",\n");
    }
    cout<<"}"<<endl;
}

void EnumBS::print_bs(blocksize_strategy bs){
    cout<<"bs = ";

    // double G1 = strategy_verification(l,BS[i].S).first;
    
    printf("(G_BKZ = %3.2f gate, B_BKZ = %3.2f bit cum-pr = %3.2f, dsvp = %3.2f, dsvp_r = %3d, G_dsvp = %3.2f gate, B_dsvp = %3.2f bit, G = %3.2f gate, B = %3.2f bit)\n",  bs.GB_BKZ.first, bs.GB_BKZ.second, bs.cum_pr, get<0>(bs.dsvp_t), get<1>(bs.dsvp_t),get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), log2(pow(2,get<2>(bs.dsvp_t)) + pow(2,bs.GB_BKZ.first)), max(bs.GB_BKZ.second,get<3>(bs.dsvp_t)) );   

    print_strategy(bs.S);
}

int EnumBS::find_pos_for_dsvp(int cdsvp){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */
    int cdsvp_tmp, pos;
    int len_BS = int(BS.size());

    cdsvp_tmp = ceil(get<0>(BS[len_BS-1].dsvp_t));
    if(cdsvp_tmp >= cdsvp)
        return len_BS;
    
    for(pos = 0; pos < len_BS; pos++){
        cdsvp_tmp = ceil(get<0>(BS[pos].dsvp_t));
        if(cdsvp_tmp < cdsvp)
            return pos;
    }

    //All elements >= dsvp_
    return pos; // pos = int(BS.size())
    
}


int EnumBS::find_pos_for_dsvp(double dsvp){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */
    
    int len_BS = int(BS.size());
    int pos; 
    double dsvp_tmp = get<0>(BS[len_BS-1].dsvp_t);
    if(dsvp_tmp >= dsvp)
        return len_BS;
    
    for(pos = 0; pos < len_BS; pos++){
        dsvp_tmp = get<0>(BS[pos].dsvp_t);
        if(dsvp_tmp < dsvp)
            return pos;
    }

    //All elements >= dsvp_
    return pos; // pos = int(BS.size())
    
}

int EnumBS::binary_search_for_cdsvp(double cdsvp){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */
    
    
    int len = int(BS.size());
    double cdsvp_tmp = round(get<0>(BS[len-1].dsvp_t)*params->enumbs_prec)/params->enumbs_prec ;
    if(cdsvp_tmp >= cdsvp)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        if(round(get<0>(BS[mid].dsvp_t)*params->enumbs_prec)/params->enumbs_prec >= cdsvp) left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}

int EnumBS::binary_search_for_G2(double G2){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */
    
    
    int len = int(BS.size());
    // double G2_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_prec)/params->enumbs_prec;
    double G2_tmp = round(get<2>(BS[len-1].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
    if(G2_tmp >= G2)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_prec)/params->enumbs_prec >= G2) left = mid + 1;
        if(round(get<2>(BS[mid].dsvp_t)*params->enumbs_prec)/params->enumbs_prec >= G2) left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}




int EnumBS::binary_search_for_dsvp(double dsvp){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */

    int len = int(BS.size());
    double dsvp_tmp = get<0>(BS[len-1].dsvp_t);
    if(dsvp_tmp >= dsvp)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        if(get<0>(BS[mid].dsvp_t)  >= dsvp) left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}




// pair<int,int> EnumBS::max_beta_jump_in_S(vector<EnumBS::strategy> S){
//     int beta_max = -1, beta_tmp;

//     for(int i =0; i< int(S.size()); i++){
//         beta_tmp = S[i].beta;
//         if(beta_max < beta_tmp)
//             beta_max = beta_tmp;
//     }
//     return beta_max;
// }

vector<double> EnumBS::extract_cdsvp(){
    vector<double> cdsvps;
    cdsvps.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        cdsvps.insert(cdsvps.end(), round(get<0>(BS[i].dsvp_t)*params->enumbs_prec)/params->enumbs_prec);
    }
    return cdsvps;
}

vector<double> EnumBS::extract_dsvp(){
    vector<double> dsvps;
    dsvps.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        dsvps.insert(dsvps.end(), get<0>(BS[i].dsvp_t));
    }
    return dsvps;
}

vector<double> EnumBS::extract_G2(){
    vector<double> G2s;
    G2s.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        G2s.insert(G2s.end(), round(get<2>(BS[i].dsvp_t)*params->enumbs_prec)/params->enumbs_prec);
    }
    return G2s;
}

bool EnumBS::no_repeated_value_verification(vector<int> nums){
    return set<int>(nums.begin(),nums.end()).size()==nums.size();
}

bool EnumBS::no_repeated_value_verification(vector<double> nums){
    return set<double>(nums.begin(),nums.end()).size()==nums.size();
}


pair<int,int> EnumBS::get_max_strategy(vector<EnumBS::strategy> S){
    if(S.size()!=0){
        strategy S_pos = S[S.size()-1];
        return make_pair(S_pos.beta, S_pos.jump);
    }
    else{
        return  make_pair(0,0);
    }
}


//true: S0's insertion range is smaller than or equal to S1.
//false: S0's insertion range is larger than S1.
bool EnumBS::compare_max_strategy(vector<EnumBS::strategy> S0, vector<EnumBS::strategy> S1){
    pair<int,int> max_strategy_bs_pos0 = get_max_strategy(S0);
    pair<int,int> max_strategy_bs_pos1 = get_max_strategy(S1);
    
    if(max_strategy_bs_pos0.first > max_strategy_bs_pos1.first || (max_strategy_bs_pos0.first == max_strategy_bs_pos1.first && max_strategy_bs_pos0.second <= max_strategy_bs_pos1.second)){
        return true;
    }
    else{
        return false;
    }
   
}

//return value: to determine whether the current bs0 is changed || not.
//False: bs0 is changed.
//True: bs0 is not change.
// pair<int,bool> EnumBS::BS_add(EnumBS::blocksize_strategy bs, int k){
//     //int cdsvp = ceil(get<0>(bs.dsvp_t));
//     double dsvp = get<0>(bs.dsvp_t);
//     double G = bs.GB_BKZ.first;
//     bool flag = true;

//     //BS.size() == 0, add bs directly
//     if(BS.size() == 0){
//         // BS.resize(1);
//         // BS[0] = bs;
//         BS.insert(BS.end(),bs);
//         return make_pair(k,flag);
//     }



//     // int pos = find_pos_for_dsvp(cdsvp);
//     // int pos = find_pos_for_dsvp(dsvp);
//     int pos = binary_search_for_cdsvp(dsvp);
    

    
//     //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
//     if(pos == 0){
//         return make_pair(k,flag);
//     }

//     //BS.size() > 0, && it exits some dsvps in EnumBS >= dsvp_, then pos > 0.
//     //cdsvp_tmp == cdvsp
//     pos--;
//     // int cdsvp_pos = ceil(get<0>(BS[pos].dsvp_t));
//     double dsvp_pos = get<0>(BS[pos].dsvp_t);
    
//     double G_pos = BS[pos].GB_BKZ.first;
    

//     // while((cdsvp_pos == cdsvp && G_pos > G) || (cdsvp_pos > cdsvp && G_pos >= G )){
//     while((dsvp_pos == dsvp && G_pos > G) || (dsvp_pos > dsvp && G_pos >= G )){
//         // if(BS[pos].S.size()==0){
//         //     flag = true;
//         //     break;
//         // }
//         BS.erase(BS.begin()+pos);
        
//         if( k > pos && k >0){
//             k -= 1;
//             flag = false;
//         }
//         pos -= 1;
//         if(pos < 0)
//             break;
//         // cdsvp_pos = ceil(get<0>(BS[pos].dsvp_t));
//         dsvp_pos = get<0>(BS[pos].dsvp_t);
//         G_pos = BS[pos].GB_BKZ.first;
//     }

//     if(pos == -1 || (dsvp_pos > dsvp && G_pos <= G)){
//         BS.insert(BS.begin()+pos+1,bs);
//         if(k==pos+1)
//             flag = false;
//     }

//     if(k==pos+1)
//         flag = false;


//     //Verification
//     // vector<int> cdsvps = extract_cdsvp();
//     // vector<int> sorted_cdsvps = cdsvps;
//     // sort(sorted_cdsvps.rbegin(),sorted_cdsvps.rend());
//     // assert(cdsvps == sorted_cdsvps);
//     // assert(no_repeated_value_verification(cdsvps));

//     if(params->debug){
//         vector<double> dsvps = extract_dsvp();
//         vector<double> sorted_dsvps = dsvps;
//         // print_vector(dsvps);
//         sort(sorted_dsvps.rbegin(),sorted_dsvps.rend());
//         assert(dsvps == sorted_dsvps);
//         assert(no_repeated_value_verification(dsvps));
//     }

//     return make_pair(k,flag);

// }


//return value: to determine whether the current bs0 is changed || not.
//False: bs0 is changed.
//True: bs0 is not change.
void EnumBS::BS_add_G2(EnumBS::blocksize_strategy bs, int k){
    
    // params->enumbs_prec = PREC;
    // while(pow(2,get<2>(bs.dsvp_t)) * params->enumbs_prec < params->enumbs_bound || pow(2,bs.GB_BKZ.first)*params->enumbs_prec < params->enumbs_bound ){
    //     // print_strategy(bs.S);
    //     // cout<<pow(2,get<2>(bs.dsvp_t))<<","<<pow(2,bs.GB_BKZ.first)<<endl;
    //     params->enumbs_prec *= 10;
    // }
    

    // double G2 = round(pow(2,get<2>(bs.dsvp_t))*params->enumbs_prec)/params->enumbs_prec;
    // double G = round(pow(2,bs.GB_BKZ.first)*params->enumbs_prec)/params->enumbs_prec;
    // double G = bs.GB_BKZ.first;


    double G2 = round(get<2>(bs.dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
    double G = round(bs.GB_BKZ.first*params->enumbs_prec)/params->enumbs_prec;


    //BS.size() == 0, add bs directly
    if(BS.size() == 0){
        // BS.resize(1);
        // BS[0] = bs;
        BS.insert(BS.end(),bs);
        return;
    }

    int pos = binary_search_for_G2(G2);

    // cout<<"========================="<<endl;
    // print_strategy(bs.S);
    // cout<<get<2>(bs.dsvp_t)<<endl;
    // cout<<get<0>(bs.dsvp_t)<<endl;
    // cout<<pow(2,get<2>(bs.dsvp_t))*params->enumbs_prec<<endl;
    // vector<double> G2s = extract_G2();
    // print_vector(G2s);
    // cout<<"pos = "<<pos<<endl;
    // printf("G2 = %3.2f ", G2);

    // cout<<"\n------------------------"<<endl;

    // if(params -> debug)
        // assert( get<2>(bs.dsvp_t) < MAX_NUM);


    
    //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
    if(pos == 0){
        return;
    }

    //BS.size() > 0, && it exits some G2s in EnumBS >= G2_, then pos > 0.
    //G2_tmp == G2
    pos--;
    // double G2_pos = round(pow(2,get<2>(BS[pos].dsvp_t))*params->enumbs_prec)/params->enumbs_prec;
    // double G_pos = round(pow(2,BS[pos].GB_BKZ.first)*params->enumbs_prec)/params->enumbs_prec;
    double G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
    double G_pos = round(BS[pos].GB_BKZ.first*params->enumbs_prec)/params->enumbs_prec;
    // double G_pos = BS[pos].GB_BKZ.first;

    // cout<<"\n========================"<<endl;
    // printf("G2_pos = %3.2f, G2 = %3.2f ", G2_pos, G2);
    // printf("G_pos = %3.2f, G = %3.2f ", G_pos, G);
    // print_strategy(BS[pos].S);
    // print_strategy(bs.S);
    // printf("\n cum_pr_pos = %e, cum_pr = %e \n", BS[pos].cum_pr, bs.cum_pr);
    // cout<<BS[pos].cum_pr<<endl;
    // cout<<bs.cum_pr<<endl;
    // cout<<"\n------------------------"<<endl;

    
    // cout<<G2_pos<<","<<G2<<endl;
    
    switch (params->enumbs_add_strategy){
        case 1: // Remain all equal strategies
            while( G2_pos > G2 && G_pos > G){
            // while((G2_pos == G2 && G_pos > G) ||  (G2_pos > G2 && G_pos >= G)){
        
                BS.erase(BS.begin()+pos);
                pos -= 1;
                if(pos == -1)
                    break;

                // G2_pos = round(pow(2,get<2>(BS[pos].dsvp_t))*params->enumbs_prec)/params->enumbs_prec;
                // G_pos = round(pow(2,BS[pos].GB_BKZ.first)*params->enumbs_prec)/params->enumbs_prec;
                // G_pos = BS[pos].GB_BKZ.first;

                G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
                G_pos = round(BS[pos].GB_BKZ.first*params->enumbs_prec)/params->enumbs_prec;

            }
            

            if( pos == -1 || (G2_pos >= G2 && G_pos <= G) ){
                BS.insert(BS.begin()+pos+1,bs);
            }

            break;

        case 2: //Delete equal strategies
            //(abs(G2_pos - MAX_NUM) <0.01)  &&
            while(  ((G2_pos == G2 && G_pos > G) ||  (G2_pos > G2 && G_pos >= G))){
        
                BS.erase(BS.begin()+pos);
                pos -= 1;
                if(pos == -1)
                    break;

                G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
                G_pos = round(BS[pos].GB_BKZ.first*params->enumbs_prec)/params->enumbs_prec;

            }
            
            //||  (abs(G2_pos - MAX_NUM) <0.01)
            if( pos == -1 || (G2_pos > G2 && G_pos <= G)  ){
                BS.insert(BS.begin()+pos+1,bs);
            }

            break;

        case 3: //Delete equal strategies, except small strategies
            //While cost of previous strategy is worse than new strategy && the range of strategy to erase is smaller than || equal to new strategy.
            //(abs(G2_pos - MAX_NUM) <0.01) && 
            while( (((G2_pos == G2 && G_pos > G) && compare_max_strategy(BS[pos].S, bs.S)) || (G2_pos > G2 && G_pos == G && compare_max_strategy(BS[pos].S, bs.S)) ||  (G2_pos > G2 && G_pos > G))){
        
                BS.erase(BS.begin()+pos);
                pos -= 1;
                if(pos == -1)
                    break;

                G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
                G_pos = round(BS[pos].GB_BKZ.first*params->enumbs_prec)/params->enumbs_prec;

            }
            
            // ||  (abs(G2_pos - MAX_NUM) <0.01) 
            if( pos == -1 || (G2_pos > G2 && G_pos <= G)){
                BS.insert(BS.begin()+pos+1,bs);
            }

            break;
        default:
            cerr << "Please choose right add strategy for EnuBS: (1/2/3)" << endl;

    }
 
    if(pos <= k)
        k = -1;

    if(params->debug){
        vector<double> G2s = extract_G2();
        vector<double> sorted_G2s = G2s;
        // print_vector(G2s);
        sort(sorted_G2s.rbegin(),sorted_G2s.rend());
        assert(G2s == sorted_G2s);
        // assert(no_repeated_value_verification(G2s));
    }
}


void EnumBS::BS_add_cdsvp(EnumBS::blocksize_strategy bs, int k){
    double cdsvp = round(get<0>(bs.dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
    double G = bs.GB_BKZ.first;

    //BS.size() == 0, add bs directly
    if(BS.size() == 0){
        // BS.resize(1);
        // BS[0] = bs;
        BS.insert(BS.end(),bs);
        return;
    }



    // int pos = find_pos_for_dsvp(cdsvp);
    int pos = binary_search_for_cdsvp(cdsvp);

    
    //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
    if(pos == 0){
        return;
    }

    //BS.size() > 0, && it exits some dsvps in EnumBS >= dsvp_, then pos > 0.
    //cdsvp_tmp == cdvsp
    pos--;
    double cdsvp_pos = round(get<0>(BS[pos].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
    double G_pos = BS[pos].GB_BKZ.first;

    // printf("\n%3.2f, %3.2f\n",cdsvp_pos, cdsvp);
    // printf("\n%3.7f, %3.7f\n",G_pos, G);
    
    // if(bs.S[bs.S.size()-1].beta<80){
    //     printf("\n\n====================================\n");
    //     printf("pos = %d, BS_size = %d", pos, int(BS.size()));
    //     printf("\n%3.7f, %3.7f\n",cdsvp_pos, cdsvp);
    //     printf("\n%3.7f, %3.7f\n",G_pos, G);
    //     print_strategy(BS[pos].S);
    //     print_strategy(bs.S);
    //     printf("...........................\n");
    // }
    //
    //
    while( ((cdsvp_pos == cdsvp && G_pos > G) || (cdsvp_pos > cdsvp && G_pos >= G)) && compare_max_strategy(BS[pos].S, bs.S)){

        if( k >= pos){
            // break;
            k = 0;
        }
        BS.erase(BS.begin()+pos);
    
        pos -= 1;
        
        if(pos == -1)
            break;
        cdsvp_pos = get<0>(BS[pos].dsvp_t);//round(get<0>(BS[pos].dsvp_t)*params->enumbs_prec)/params->enumbs_prec;
        G_pos = BS[pos].GB_BKZ.first;
    }

    

    // print_bs(bs);
    if(pos == -1 || (cdsvp_pos > cdsvp && G_pos <= G)){
    // if(pos == -1 || (cdsvp_pos >= cdsvp && G_pos <= G)){
    // if(pos == -1 || (cdsvp_pos > cdsvp && G_pos <= G) || cdsvp_pos == cdsvp){
        BS.insert(BS.begin()+pos+1,bs);
    }
    

    //Verification
    if(params->debug){
        vector<double> cdsvps = extract_cdsvp();
        vector<double> sorted_cdsvps = cdsvps;
        sort(sorted_cdsvps.rbegin(),sorted_cdsvps.rend());
        print_vector(cdsvps);
        assert(cdsvps == sorted_cdsvps);
        assert(no_repeated_value_verification(cdsvps));


        // vector<double> G2s = extract_G2();
        // vector<double> sorted_G2s = G2s;
        // print_vector(G2s);
        // sort(sorted_G2s.rbegin(),sorted_G2s.rend());
        // assert(G2s == sorted_G2s);
        
    }


}


tuple<double,int,double,double> EnumBS::pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB, double &cum_pr, int beta, int jump){

    //simulate pnj-bkz more precisely
    int f, beta_, d = l.size();
    f = default_dim4free_fun(beta);
    if(jump <= 2)
        beta_ = beta;
    else if(jump >=3 && jump <=4)
        beta_ = get_beta_from_sieve_dim(beta-f,d,2);
    else if(jump>=5)
        beta_ = get_beta_from_sieve_dim(beta-f,d,1);



    double rem_pr = 1. - cum_pr;
    sim -> simulate(l,l,beta_,jump,1);

    boost::math::chi_squared chisquare(beta_);
    double pr = boost::math::cdf(chisquare,pow(2,2.*l[d-beta_]));
    

    pair<double,double> GB = cost->bkz_cost(d,beta,jump,params->cost_model);
    // if(beta == 178)
    //     printf("beta = %d, l[i]= %e, G = %e, rem_pr = %e, pr = %e\n", beta, pow(2,2.*l[d-beta]), GB.first, rem_pr, pr);
    
    cum_GB.first = log2(pow(2,cum_GB.first)+(pow(2,GB.first)*rem_pr*pr));
    // cum_GB.first = log2(pow(2,cum_GB.first)+(pow(2,GB.first)));
    cum_GB.second = max(cum_GB.second, GB.second);

    // printf("cum_G = %e\n", cum_GB.first );
    cum_pr += rem_pr * pr;
    rem_pr = 1. - cum_pr;


    return dsvp_predict(l, cum_pr, cost,params->cost_model, params->progressive_sieve);

}


void EnumBS::max_tour_for_pnjbkz_beta(int k, int beta,int jump){
    EnumBS::blocksize_strategy bs = BS[k];
    tuple<double,int,double,double> dsvp_t1;
    double G20 =  get<2>(bs.dsvp_t), G21;
    // double G20 = get<2>(bs.dsvp_t), G21;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;

    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB=bs.GB_BKZ;

    int  loop = 0;
    

    dsvp_t1 = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, jump);
    
    G21 = get<2>(dsvp_t1);

    assert(G21 >= 0.);

    while((abs(G20 - MAX_NUM)<0.01 or G20 > G21) && loop < params->max_loop ){

        loop +=1;
        G20 = G21;
        
        if(loop ==1){
            S.insert(S.end(),{beta, jump, loop});
        }
        else
            // S[len_S-1].tours = loop;
            S[S.size()-1].tours = loop;
        bs = {dsvp_t1, S, l, cum_GB, cum_pr};



        if(params->verification){
            pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
            //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
            //cerr<<"cum_G="<<cum_GB.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
            assert(abs(verified_cum_G_pr.first-cum_GB.first)<0.001);
            assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
        }
    
        EnumBS::BS_add_G2(bs, k);

        if(cum_pr >= 0.999)
            break;
      

        dsvp_t1 = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, jump);
        G21 = get<2>(dsvp_t1);
        assert(G21 >= 0.);
        // G21 = get<2>(dsvp_t1);
    }
}


// void EnumBS::max_tour_for_pnjbkz_beta_G2(int k, int beta,int jump){
//     EnumBS::blocksize_strategy bs = BS[k];

//     tuple<double,int,double,double> dsvp_t1;
//     // double G20 =  get<2>(bs.dsvp_t), G21;
//     double G20 = get<2>(bs.dsvp_t), G21;
//     vector<double> l = bs.l; 
//     vector<strategy> S = bs.S;

//     double cum_pr = bs.cum_pr;
//     pair<double,double> cum_GB=bs.GB_BKZ;

//     int  loop = 0;

//     dsvp_t1 = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, jump);
    
//     // G21 = get<2>(dsvp_t1);
//     G21 = get<2>(dsvp_t1);

//     // while(G20 - G21 >= 1){
//     while(G20 > G21 && loop < params->max_loop && cum_pr <= 0.999 ){
//         // if(loop > 5)
//         //     printf("\n%d, %d, %d, %3.2f, %3.2f\n", beta, jump, loop, G20,G21);
//         loop +=1;
//         // G20 = G21;
//         G20 = G21;
        
//         if(loop ==1){
//             // len_S+=1;
//             // S.resize(len_S);
//             // S[len_S-1] = {beta, jump, loop};
//             S.insert(S.end(),{beta, jump, loop});
          
//         }
//         else
//             // S[len_S-1].tours = loop;
//             S[S.size()-1].tours = loop;
//         bs = {dsvp_t1, S, l, cum_GB, cum_pr};

//         if(params->verification){
//             pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
//             //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
//             //cerr<<"cum_G="<<cum_GB.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
//             assert(abs(verified_cum_G_pr.first-cum_GB.first)<0.001);
//             assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
//         }
    
//         // k_flag = EnumBS::BS_add(bs, k);
//         EnumBS::BS_add_G2(bs, k);

//         dsvp_t1 = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, jump);
//         // G21 = get<2>(dsvp_t1);
//         G21 = get<2>(dsvp_t1);
//     }
// }


// pair<int,bool> EnumBS::max_tour_for_pnjbkz_beta_G2_backup(int k, int beta,int jump){
//     EnumBS::blocksize_strategy bs = BS[k];
//     pair<int,bool> k_flag = make_pair(k,true);
//     tuple<double,int,double,double> dsvp_t1;
//     // double G20 =  get<2>(bs.dsvp_t), G21;
//     double G20 = get<2>(bs.dsvp_t), G21;
//     vector<double> l = bs.l; 
//     vector<strategy> S = bs.S;

//     double cum_pr = bs.cum_pr;
//     pair<double,double> cum_GB=bs.GB_BKZ;

//     int  loop = 0;

//     dsvp_t1 = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, jump);
    
//     // G21 = get<2>(dsvp_t1);
//     G21 = get<2>(dsvp_t1);

//     // while(G20 - G21 >= 1){
//     while(G20 - G21 >= 1){
//         // if(loop > 5)
//         //     printf("\n%d, %d, %d, %3.2f, %3.2f\n", beta, jump, loop, G20,G21);
//         loop +=1;
//         // G20 = G21;
//         G20 = G21;
        
//         if(loop ==1){
//             // len_S+=1;
//             // S.resize(len_S);
//             // S[len_S-1] = {beta, jump, loop};
//             S.insert(S.end(),{beta, jump, loop});
          
//         }
//         else
//             // S[len_S-1].tours = loop;
//             S[S.size()-1].tours = loop;
//         bs = {dsvp_t1, S, l, cum_GB, cum_pr};

//         if(params->verification){
//             pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
//             //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
//             //cerr<<"cum_G="<<cum_GB.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
//             assert(abs(verified_cum_G_pr.first-cum_GB.first)<0.001);
//             assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
//         }
    
//         // k_flag = EnumBS::BS_add(bs, k);
//         // k_flag = EnumBS::BS_add_G2(bs, k);
//         EnumBS::BS_add_G2(bs, k);

//         dsvp_t1 = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, jump);
//         // G21 = get<2>(dsvp_t1);
//         G21 = get<2>(dsvp_t1);
//     }
//     return k_flag;
// }



// void EnumBS::max_tour_for_pnjbkz_beta_in_parallel( int beta_j_t_id_begin, vector<pair<int,int>> beta_j_tid,  int k){
//     // bool flag = true;

//     for(int i = 0; i< int(beta_j_tid.size()); i++){
//         EnumBS::blocksize_strategy bs = BS[k];
//         tuple<double,int,double,double> dsvp_t1;
//         double G20 =  get<2>(bs.dsvp_t), G21;
//         vector<double> l = bs.l; 
//         vector<strategy> S = bs.S;

//         double cum_pr = bs.cum_pr;
//         pair<double,double> cum_GB=bs.GB_BKZ;

//         int beta = beta_j_tid[i].first, jump = beta_j_tid[i].second;
    
//         int f = dims4free(beta);
//         if(f == 0 && jump > 1)
//             return;
//         if(f!=0 && jump >= f)
//             return;
    
//         // if(jump > beta)
//         //     return;


//         int index = beta_j_t_id_begin + i;
        
//         int  loop = 0;
        
//         dsvp_t1 = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, jump);
        
//         G21 = get<2>(dsvp_t1);

//         assert(G21 >= 0.);

//         tmpBS[index].clear();

//         while(G20 > G21 && loop < params->max_loop && cum_pr <= 0.999 ){
//             loop +=1;
//             G20 = G21;
//             if(loop ==1){
//                 S.insert(S.end(),{beta, jump, loop});
            
//             }
//             else
//                 S[S.size()-1].tours = loop;
                
//             bs = {dsvp_t1, S, l, cum_GB, cum_pr};
            
//             if(params->verification){
//                 pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
//                 assert(abs(verified_cum_G_pr.first-cum_GB.first)<0.001);
//                 assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
//             }
        
//             tmpBS[index].insert(tmpBS[index].end(),bs);

//             // if(!BS_add_determine(bs, k)){
//             //     flag = false;
//             // }


//             dsvp_t1 = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, jump);
//             G21 = get<2>(dsvp_t1);
//             assert(G21 >= 0.);
//         }

//         // if(!flag)
//         //     return flag;
        
//     }
//     // return flag;
// }




void EnumBS::max_tour_for_pnjbkz_beta_in_parallel( int beta_j_t_id_begin, vector<pair<int,int>> beta_j_tid,  int k){
    // bool flag = true;

    for(int i = 0; i< int(beta_j_tid.size()); i++){
        EnumBS::blocksize_strategy bs = BS[k];
        tuple<double,int,double,double> dsvp_t1;
        double G20 = get<2>(bs.dsvp_t), G21;
        vector<double> l = bs.l; 
        vector<strategy> S = bs.S;

        double cum_pr = bs.cum_pr;
        pair<double,double> cum_GB=bs.GB_BKZ;

        int beta = beta_j_tid[i].first, jump = beta_j_tid[i].second;
    
        int f = dims4free(beta);
        if(f == 0 && jump > 1)
            return;
        if(f!=0 && jump >= f)
            return;

        int index = beta_j_t_id_begin + i;
        
        int  loop = 0;
        
        dsvp_t1 = pnjbkz_beta_loop(l, cum_GB, cum_pr, beta, jump);
        
        G21 = get<2>(dsvp_t1);

        // assert(G21 >= 0.);

        tmpBS[index].clear();
        // cerr<< "beta = "<< beta << ", jump = "<< jump <<endl;
        // cerr<< "G20 = "<< G20 <<", G21 = " << G21 <<", MAX_NUM = " << MAX_NUM<<","<<<<endl;
        while((abs(G20 - MAX_NUM)<0.01 or G20 > G21) && loop < params->max_loop){
            loop +=1;
            G20 = G21;
            if(loop ==1){
                S.insert(S.end(),{beta, jump, loop});
            }
            else
                S[S.size()-1].tours = loop;
                
            bs = {dsvp_t1, S, l, cum_GB, cum_pr};
            
            if(params->verification){
                pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
                assert(abs(verified_cum_G_pr.first-cum_GB.first)<0.001);
                assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
            }
        
            tmpBS[index].insert(tmpBS[index].end(),bs);

            if(cum_pr >= 0.999)
                break;

            dsvp_t1 = pnjbkz_beta_loop( l, cum_GB, cum_pr, beta, jump);
            G21 = get<2>(dsvp_t1);
            // assert(G21 >= 0.);
        }

        // if(!flag)
        //     return flag;
        
    }
    // return flag;
}


void EnumBS::enumbs_est(vector<double> l){
    /*
    input: l -- gs-lengths; 
    Return: Optimal strategy to minimize solving cost.
    */
    l0 = l;
    int beta_start, k = 0, d = l0.size();
    blocksize_strategy bs;

    tuple<double,int,double,double>  dsvp0_t = dsvp_predict(l0, 0., cost, params->cost_model, params->progressive_sieve);

    BS.insert(BS.end(),{dsvp0_t, {},l0,make_pair(0.,0.),0.});

    int j_start,len_S;
    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;
    
    while( k < int(BS.size())){
        bs = BS[k];
        len_S = bs.S.size();

        if(bs.cum_pr >= 0.999)
            continue;
    
        if(len_S == 0){
            beta_start = 51;
            j_start = params->J;
        }
        else{
            j_start = bs.S[len_S-1].jump;
            beta_start = bs.S[len_S-1].beta;
            if(j_start==1){
                beta_start += 1;
                j_start = params->J;
            }
            else
                j_start -= params->J_gap;
        }
        
        // if(params->verbose)
        auto start = system_clock::now();
        //int(0.9*d)
        for(int beta = beta_start; beta < min(params->max_dim, d); beta +=params->gap){
            for(int j = j_start; j>0; j-= params->J_gap){
                // k_flag = EnumBS::max_tour_for_pnjbkz_beta(k,beta,j); 
                // k_flag = EnumBS::max_tour_for_pnjbkz_beta_G2(k,beta,j); 
                // EnumBS::max_tour_for_pnjbkz_beta_G2(k,beta,j); 

                int f = dims4free(beta);
                if((f == 0 && j > 1) or (f!=0 && j >= f))
                    continue;

                EnumBS::max_tour_for_pnjbkz_beta(k,beta,j); 
                

                // if(!k_flag.second){
                //     k = k_flag.first;
                //     if(params->verbose){
                //         auto finish = system_clock::now();
                //         duration<double> diff = finish - start;
                //         printf("\r index: %4d, (beta,j): (%4d,%4d) --> (%4d,%4d), goal index: %4d, cost = " ,k+1,beta_start,j_start,(min(params->max_dim,d)-1-beta_start)/params->gap*params->gap+beta_start,1,int(BS.size()));
                //         cerr<<setprecision(2)<<diff.count()<<'s';
                //     }
                //     goto WHILE_START;
                // }
            }
            j_start = params->J;
            
        }
        if(params->verbose){
            auto finish = system_clock::now();
            duration<double> diff = finish - start;
            printf("\r index: %4d, (beta,j): (%4d,%4d) --> (%4d,%4d), goal index: %4d, cost = " ,k+1,beta_start,j_start,(min(params->max_dim,d)-1-beta_start)/params->gap*params->gap+beta_start,1,int(BS.size()));
            cout<<setprecision(2)<<diff.count()<<'s';
        }

        k ++;
    }
    printf("\n");
    //print_BS(BS);

    //Find the optimal strategy
    double Gmin = MAX_NUM, Bmin  = MAX_NUM, G1, G2, G,  B;
    EnumBS::blocksize_strategy bsmin;
    for(int i = 0; i<int(BS.size()); i++){
        G1 = BS[i].GB_BKZ.first;
        G2 = get<2>(BS[i].dsvp_t);
        G = log2(pow(2,G1)+pow(2,G2));
        B = max(get<3>(BS[i].dsvp_t),BS[i].GB_BKZ.second);
        if(G<Gmin){
            bsmin = BS[i];
            Gmin = G;
            Bmin = B;
        }
    }
    printf("Find the optimal Strategy through EumBS!!\n");
    print_bs(bsmin);
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log2(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec), Memory Cost = %3.2f log2(bit)\n", Gmin, Bmin);
}

//return value: to determine whether the current bs0 will be changed || not.
//False: bs0 will be changed.
//True: bs0 will not change.
// bool EnumBS::BS_add_determine(EnumBS::blocksize_strategy bs, int k){
//     //int cdsvp = ceil( get<2>(bs.dsvp_t));
//     double dsvp =  get<2>(bs.dsvp_t);
//     double G = bs.GB_BKZ.first;

//     //BS.size() == 0, add bs directly
//     if(BS.size() == 0){
//         return true;
//     }



//     // int pos = find_pos_for_dsvp(cdsvp);
//     // int pos = find_pos_for_dsvp(dsvp);
//     int pos = binary_search_for_dsvp(dsvp);
    

    
//     //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
//     if(pos == 0){
//         return true;
//     }

//     //BS.size() > 0, && it exits some dsvps in EnumBS >= dsvp_, then pos > 0.
//     //cdsvp_tmp == cdvsp
//     pos--;
//     // int cdsvp_pos = ceil(get<0>(BS[pos].dsvp_t));
//     double dsvp_pos = get<0>(BS[pos].dsvp_t);
    
//     double G_pos = BS[pos].GB_BKZ.first;
    

//     // while((cdsvp_pos == cdsvp && G_pos > G) || (cdsvp_pos > cdsvp && G_pos >= G )){
//     while((dsvp_pos == dsvp && G_pos > G) || (dsvp_pos > dsvp && G_pos >= G )){
//         if( k > pos && k >0){
//             k -= 1;
//             return false;
//         }
//         pos -= 1;
//         if(pos < 0)
//             break;
//         // cdsvp_pos = ceil(get<0>(BS[pos].dsvp_t));
//         dsvp_pos = get<0>(BS[pos].dsvp_t);
//         G_pos = BS[pos].GB_BKZ.first;
//     }


//     if(pos == -1 || (dsvp_pos > dsvp && G_pos <= G)){
//         if(k==pos+1){
//             return false;
//         }
//     }
    


//     //Verification
//     // vector<int> cdsvps = extract_cdsvp();
//     // vector<int> sorted_cdsvps = cdsvps;
//     // sort(sorted_cdsvps.rbegin(),sorted_cdsvps.rend());
//     // assert(cdsvps == sorted_cdsvps);
//     // assert(no_repeated_value_verification(cdsvps));

//     if(params->debug){
//         vector<double> dsvps = extract_dsvp();
//         vector<double> sorted_dsvps = dsvps;
//         // print_vector(dsvps);
//         sort(sorted_dsvps.rbegin(),sorted_dsvps.rend());
//         assert(dsvps == sorted_dsvps);
//         assert(no_repeated_value_verification(dsvps));
//     }

//     return true;

// }



void EnumBS::enumbs_est_in_parallel(vector<double> l){
    /*
    input: l -- gs-lengths;
    Return: Optimal strategy to minimize solving cost.
    */
    set_threads(params->threads);
    

    l0 = l;
    int beta_start, k = 0, d = l0.size(),beta;
    blocksize_strategy bs;

    tuple<double,int,double,double>  dsvp0_t = dsvp_predict(l0, 0., cost,params->cost_model, params->progressive_sieve);

    
    // BS.resize(1);
    // BS[0] =  {dsvp0_t, {},l0,make_pair(0.,0.),0.};;
    BS.insert(BS.end(),{dsvp0_t, {},l0,make_pair(0.,0.),0.});

    
    int j,len_S;
    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

    while( k < int(BS.size())){
        bs = BS[k];
        len_S = bs.S.size();
 
        if(bs.cum_pr >= 0.999){
            k++;
            continue;
        }
      
    
        if(len_S == 0){
            beta_start = 51;
            j = params->J;
        }
        else{
            j = bs.S[len_S-1].jump;
            beta_start = bs.S[len_S-1].beta;
            if(j==1){
                beta_start += 1;
                j = params->J;
            }
            else
                j -= params->J_gap;
        }

        //int(0.9*d)
        int t_id = 0, len = (((j-1)/params->J_gap)+1)+(min(params->max_dim,d)-1-beta_start)/params->gap*(((params->J-1)/params->J_gap)+1);

        if(len > 0){

            int threads = min(len, params->threads);
            int block = len/threads;
            vector<vector<pair<int,int>>> beta_j;
            beta_j.resize(threads);
            tmpBS.resize(len);
            vector<int> beta_j_t_id_begins(threads,0), departs(threads,0);
            
        
            for(int t_id = 0; t_id < threads; t_id++){
                if(len%(threads)!=0 && ((t_id/(len%threads)) ? 0 : 1)){
                    departs[t_id] = block + 1;
                }
                else{
                    departs[t_id] = block;
                }
                if(t_id > 0)
                    beta_j_t_id_begins[t_id] = beta_j_t_id_begins[t_id-1]+departs[t_id-1];
            } 
        
        
            //int(0.9*d)
            for(beta = beta_start; beta < min(params->max_dim, d); beta +=params->gap){
                for(; j>0; j-= params->J_gap){
                    if(int(beta_j[t_id].size()) ==  departs[t_id] && t_id < threads-1){
                        t_id++;
                    }
                    beta_j[t_id].insert(beta_j[t_id].end(),make_pair(beta,j));
                }
                j = params->J;
            }
        
            if(params->debug){
                int Sum = 0;
                for(int t_id = 0; t_id  < int(beta_j.size()); t_id ++){
                    Sum += beta_j[t_id ].size();
                    // cerr<<"t_id = "<<t_id<<", beta_j[i].size() = "<<beta_j[t_id].size()<<endl;
                }
                // cerr<<len<<","<<Sum<<endl;
                assert(len == Sum);
            }
        
            
            auto start= system_clock::now();
            for (t_id = 0; t_id < threads; t_id++){
                threadpool.push([this, t_id, beta_j_t_id_begins, beta_j, k](){
                    max_tour_for_pnjbkz_beta_in_parallel(  beta_j_t_id_begins[t_id], beta_j[t_id], k); 
                    // cerr<<"---"<<t_id<<"--"<<endl;
                });
                // if(!flag[t_id]){
                //     bool flags = true;
                //     for(int i = 0; i<=t_id; i++){
                //         if(flag[i])
                //     }
                //     break;
                // }
                
            }
            threadpool.wait_work(); 
            // printf("==============");
            
            if(params->verbose){
                auto finish = system_clock::now();
                duration<double> diff = finish - start;
                printf("\r index: %4d, (beta,j): (%4d,%4d) --> (%4d,%4d), goal index: %4d, cost =" ,k+1,beta_j[0][0].first,beta_j[0][0].second,beta_j[threads-1][beta_j[threads-1].size()-1].first,beta_j[threads-1][beta_j[threads-1].size()-1].second,int(BS.size()));
                cerr<<setprecision(2)<<diff.count()<<'s';
            }
            //BS_add method.
            // if(params->verbose)
            //     start = clock();
            // start = system_clock::now();
           
            for(int i = 0; i< len; i++){
                // printf("\nIn BS_add: i = %d, tmpBS_size = %d \n ", i,int(tmpBS[i].size()));
                
                for(int ii = 0; ii < int(tmpBS[i].size()); ii++){
                    // print_strategy(tmpBS[i][ii].S);
                    // EnumBS::BS_add_cdsvp(tmpBS[i][ii], k);
                    EnumBS::BS_add_G2(tmpBS[i][ii], k);
                    
                }
            
            }
            
            tmpBS.clear();
            beta_j.clear();
        }
        k++;
    }
    printf("\n");

    //Find the optimal strategy
    double Gmin = MAX_NUM, Bmin  = MAX_NUM, G1, G2, G,  B;
    EnumBS::blocksize_strategy bsmin;
    for(int i = 0; i<int(BS.size()); i++){
        // G1 = strategy_verification(l,BS[i].S).first;
        G1 = BS[i].GB_BKZ.first;
        G2 = get<2>(BS[i].dsvp_t);
        G = log2(pow(2,G1)+pow(2,G2));
        B = max(get<3>(BS[i].dsvp_t),BS[i].GB_BKZ.second);
        if(G<Gmin){
            bsmin = BS[i];
            Gmin = G;
            Bmin = B;
        }
    }
    printf("Find the optimal Strategy through EumBS!!\n");
    print_bs(bsmin);
    
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec), Memory Cost = %3.2f log(bit)\n", Gmin, Bmin);
}



pair<double,double> EnumBS::strategy_verification(vector<double> l,vector<strategy> S){

    int d = l.size();
    double cum_pr = 0., rem_pr = 1., proba, G1cum=0., B1cum = 0.;
    BKZJSim* sim = new BKZJSim();
    COST* cost = new COST();
    for(int i = 0; i< int(S.size()); i++){
        EnumBS::strategy bs = S[i];
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
