#include "enumbs.h"
#include <unistd.h>

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
    cout<<"{";
    for(int i =0; i< int(BS.size()); i++){
        //print_strategy(BS[i].S);
        print_bs(BS[i]);
    }
    cout<<"}"<<endl;
}

void EnumBS::print_bs(blocksize_strategy bs){
    cout<<"bs = ";

    // double G1 = strategy_verification(l,BS[i].S).first;
    if(params->cost_model == 1)
        printf("(slope = %e, G_BKZ = %e gate, B_BKZ = %e bit cum-pr = %e, dsvp = %e, dsvp_r = %3d, G_dsvp = %e gate, B_dsvp = %e bit, G = %e gate, B = %e bit, min_GB.first = %e gate)\n",  bs.slope, bs.cum_avg_GB_BKZ.first, bs.cum_avg_GB_BKZ.second, bs.cum_pr, get<0>(bs.dsvp_t), get<1>(bs.dsvp_t),get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), log2(pow(2,get<2>(bs.dsvp_t)) + pow(2,bs.cum_avg_GB_BKZ.first)), max(bs.cum_avg_GB_BKZ.second,get<3>(bs.dsvp_t)), bs.min_GB.first);  
    if(params->cost_model == 2)
        printf("(slope = %e, G_BKZ = %e sec, B_BKZ = %e bit cum-pr = %e, dsvp = %e, dsvp_r = %3d, G_dsvp = %e sec, B_dsvp = %e bit, G = %e sec, B = %e bit,  min_GB.first = %e gate)\n",  bs.slope, bs.cum_avg_GB_BKZ.first, bs.cum_avg_GB_BKZ.second, bs.cum_pr, get<0>(bs.dsvp_t), get<1>(bs.dsvp_t),get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), log2(pow(2,get<2>(bs.dsvp_t)) + pow(2,bs.cum_avg_GB_BKZ.first)), max(bs.cum_avg_GB_BKZ.second,get<3>(bs.dsvp_t)), bs.min_GB.first);  
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
    double cdsvp_tmp = round(get<0>(BS[len-1].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec ;
    if(cdsvp_tmp >= cdsvp)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        if(round(get<0>(BS[mid].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec >= cdsvp) left = mid + 1;
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
    // double G2_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
    // double G2_tmp = round(get<2>(BS[len-1].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
    double G2_tmp = get<2>(BS[len-1].dsvp_t);

    if(G2_tmp >= G2)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        // if(round(get<2>(BS[mid].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        if( get<2>(BS[mid].dsvp_t) >= G2 ) left = mid + 1;
        
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}


int EnumBS::binary_search_for_G(double G){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */
    
    
    int len = int(BS.size());
    // double G_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
    double G_tmp = BS[len-1].GB.first;
    if(G_tmp >= G)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        if(BS[mid].GB.first >= G) left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}






int EnumBS::binary_search_for_slope(double slope){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */
    
    
    int len = int(BS.size());
    // double G2_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
    double slope_tmp = BS[len-1].slope;
    if(slope_tmp <= slope)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        if(BS[mid].slope <= slope) left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}



int EnumBS::binary_search_for_G2_slope_cum_pr(blocksize_strategy bs){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    # sort in (G2,slope,cum_pr). Find the position to put (G2,slope,cum_pr), G2-->smaller, slope --> larger, cum_pr --> larger.
    */
    
    
    int len = int(BS.size());
    // double G2_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
 
    
    if( get<2>(BS[len-1].dsvp_t) > get<2>(bs.dsvp_t) || (get<2>(BS[len-1].dsvp_t) == get<2>(bs.dsvp_t) && BS[len-1].slope < bs.slope) || (get<2>(BS[len-1].dsvp_t) == get<2>(bs.dsvp_t) && BS[len-1].slope == bs.slope && BS[len-1].cum_pr <= bs.cum_pr)) 
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        if( get<2>(BS[mid].dsvp_t) > get<2>(bs.dsvp_t) || (get<2>(BS[mid].dsvp_t) == get<2>(bs.dsvp_t) && BS[mid].slope < bs.slope) || (get<2>(BS[mid].dsvp_t) == get<2>(bs.dsvp_t) && BS[mid].slope == bs.slope && BS[mid].cum_pr <= bs.cum_pr)) 
            left = mid + 1;
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
        cdsvps.insert(cdsvps.end(), round(get<0>(BS[i].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec);
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
        // G2s.insert(G2s.end(), round(get<2>(BS[i].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec);
        G2s.insert(G2s.end(),get<2>(BS[i].dsvp_t));
    }
    return G2s;
}


vector<pair<double,double>> EnumBS::extract_G2G1(){
    vector<pair<double,double>> G2s;
    G2s.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        // G2s.insert(G2s.end(), round(get<2>(BS[i].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec);
        G2s.insert(G2s.end(), make_pair(get<2>(BS[i].dsvp_t), BS[i].cum_avg_GB_BKZ.first));
    }
    return G2s;
}


vector<double> EnumBS::extract_slope(){
    vector<double> slopes;
    slopes.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        slopes.insert(slopes.end(), BS[i].slope);
    }
    return slopes;
}


vector<double> EnumBS::extract_cum_pr(){
    vector<double> cum_prs;
    cum_prs.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        cum_prs.insert(cum_prs.end(), BS[i].cum_pr);
    }
    return cum_prs;
}

vector<tuple<double,double,double>> EnumBS::extract_G2_slope_cum_pr(){
    vector<tuple<double,double,double>> basis_quality_list;
    basis_quality_list.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        basis_quality_list.insert(basis_quality_list.end(), make_tuple(get<2>(BS[i].dsvp_t), BS[i].slope, BS[i].cum_pr));
    }
    return basis_quality_list;
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
//     double G = bs.cum_avg_GB_BKZ.first;
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
    
//     double G_pos = BS[pos].cum_avg_GB_BKZ.first;
    

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
//         G_pos = BS[pos].cum_avg_GB_BKZ.first;
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



void EnumBS::BS_add_G2(EnumBS::blocksize_strategy bs, int k){

    // print_BS(BS);
    if(BS.size() == 0){
        BS.insert(BS.end(),bs);
        return;
    }

    int pos = binary_search_for_G2(get<2>(bs.dsvp_t));



    if(pos == 0){
        return;
    }


    pos--;

    if(params->debug){
        if(bs.GB.second < params->max_RAM){
            cout<<"====Add===="<<endl;
            printf("pos=%d\n",pos);
            print_bs(bs);
            print_bs(BS[pos]);
            cout<<"......."<<endl;
            usleep(10000000);
        }
    }



    // while(pos > k && ( get<2>(bs.dsvp_t) < get<2>(BS[pos].dsvp_t) + params->enumbs_G_prec || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && get<2>(bs.dsvp_t) == get<2>(BS[pos].dsvp_t)  + params->enumbs_G_prec )) && bs.cum_avg_GB_BKZ.first < BS[pos].cum_avg_GB_BKZ.first + params->enumbs_G_prec && compare_max_strategy(BS[pos].S, bs.S)){
    while(pos > k && ( get<2>(bs.dsvp_t) < get<2>(BS[pos].dsvp_t) + params->enumbs_G_prec || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && get<2>(bs.dsvp_t) == get<2>(BS[pos].dsvp_t)  + params->enumbs_G_prec )) && bs.cum_avg_GB_BKZ.first < BS[pos].cum_avg_GB_BKZ.first + params->enumbs_G_prec){

        // if(params->debug){
        //     if(bs.S.size()>0){
        //         cout<<"====Delete1===="<<endl;
        //         printf("pos=%d\n",pos);
        //         print_bs(bs);
        //         print_bs(BS[pos]);
        //         cout<<"......."<<endl;
        //         usleep(10000000);
        //     }
        // }

        BS.erase(BS.begin()+pos);

        pos--;
    }

    
    BS.insert(BS.begin()+pos+1,bs);
    
            

    
    // while( pos+1<int(BS.size())-1 && pos+1 > k && ( get<2>(BS[pos+2].dsvp_t) < get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec  &&  get<2>(BS[pos+2].dsvp_t) == get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec)) && BS[pos+2].cum_avg_GB_BKZ.first < BS[pos+1].cum_avg_GB_BKZ.first + params->enumbs_G_prec && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){
    while( pos+1<int(BS.size())-1 && pos+1 > k && ( get<2>(BS[pos+2].dsvp_t) < get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec  &&  get<2>(BS[pos+2].dsvp_t) == get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec)) && BS[pos+2].cum_avg_GB_BKZ.first < BS[pos+1].cum_avg_GB_BKZ.first + params->enumbs_G_prec){
        // if(params->debug){
        //     if(BS[pos+1].S.size()>0){
        //         cout<<"====Delete2===="<<endl;
        //         printf("pos=%d\n",pos+1);
        //         print_bs(BS[pos+1]);
        //         print_bs(BS[pos+2]);
        //         cout<<"......."<<endl;
        //         // print_vector(extract_G2G1());
        //         usleep(10000000);
        //     }
        // }

        // while (pos+1<BS.size()-1 && BS[pos+2].slope >= BS[pos+1].slope - params->enumbs_slope_prec  && BS[pos+2].cum_avg_GB_BKZ.first <= BS[pos+1].cum_avg_GB_BKZ.first +  params->enumbs_G_prec && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){
        BS.erase(BS.begin()+pos+1);
        pos--;
    }
                
           
    

 
    if(pos <= k)
        k = -1;

    if(params->debug){
        
        vector<double> G2s = extract_G2();
        vector<double> sorted_G2s = G2s;
        sort(sorted_G2s.rbegin(),sorted_G2s.rend());
        assert(G2s == sorted_G2s);
        // print_vector(G2s);
        // assert(no_repeated_value_verification(G2s));
    }
}



// void EnumBS::BS_add_G2(EnumBS::blocksize_strategy bs, int k){
    
//     // params->enumbs_G_prec = PREC;
//     // while(pow(2,get<2>(bs.dsvp_t)) * params->enumbs_G_prec < params->enumbs_bound || pow(2,bs.cum_avg_GB_BKZ.first)*params->enumbs_G_prec < params->enumbs_bound ){
//     //     // print_strategy(bs.S);
//     //     // cout<<pow(2,get<2>(bs.dsvp_t))<<","<<pow(2,bs.cum_avg_GB_BKZ.first)<<endl;
//     //     params->enumbs_G_prec *= 10;
//     // }
    

//     // double G2 = round(pow(2,get<2>(bs.dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
//     // double G = round(pow(2,bs.cum_avg_GB_BKZ.first)*params->enumbs_G_prec)/params->enumbs_G_prec;
//     // double G = bs.cum_avg_GB_BKZ.first;


//     double G2 = get<2>(bs.dsvp_t); //round(get<2>(bs.dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//     double G = bs.cum_avg_GB_BKZ.first; //round(bs.cum_avg_GB_BKZ.first*params->enumbs_G_prec)/params->enumbs_G_prec;


//     //BS.size() == 0, add bs directly
//     if(BS.size() == 0){
//         // BS.resize(1);
//         // BS[0] = bs;
//         BS.insert(BS.end(),bs);
//         return;
//     }

//     int pos = binary_search_for_G2(G2);

//     // cout<<"========================="<<endl;
//     // // print_strategy(bs.S);
//     // print_bs(bs);
//     // if(pos != 0)
//     //     print_bs(BS[pos-1]);
//     // cout<<get<2>(bs.dsvp_t)<<endl;
//     // cout<<get<0>(bs.dsvp_t)<<endl;
//     // vector<double> G2s = extract_G2();
//     // print_vector(G2s);
//     // cout<<"pos = "<<pos<<endl;
//     // printf("G2 = %e ", G2);
//     // cout<<params->max_num<<endl;

//     // cout<<"\n------------------------"<<endl;

//     // if(params -> debug)
//         // assert( get<2>(bs.dsvp_t) < params->max_num);


    
//     //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
//     if(pos == 0){
//         return;
//     }

    

//     //BS.size() > 0, && it exits some G2s in EnumBS >= G2_, then pos > 0.
//     //G2_tmp == G2
//     pos--;
//     // double G2_pos = round(pow(2,get<2>(BS[pos].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
//     // double G_pos = round(pow(2,BS[pos].cum_avg_GB_BKZ.first)*params->enumbs_G_prec)/params->enumbs_G_prec;
//     double G2_pos = get<2>(BS[pos].dsvp_t); //round(get<2>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//     double G_pos = BS[pos].cum_avg_GB_BKZ.first; //round(BS[pos].cum_avg_GB_BKZ.first*params->enumbs_G_prec)/params->enumbs_G_prec;

//     // double slope_pos = BS[pos].slope;
//     // double slope = bs.slope;

//     // double G_pos = BS[pos].cum_avg_GB_BKZ.first;

//     // cout<<"\n========================"<<endl;
//     // printf("G2_pos = %e, G2 = %e ", G2_pos, G2);
//     // printf("G_pos = %e, G = %e ", G_pos, G);
//     // print_strategy(BS[pos].S);
//     // print_strategy(bs.S);
//     // printf("\n cum_pr_pos = %e, cum_pr = %e \n", BS[pos].cum_pr, bs.cum_pr);
//     // cout<<BS[pos].cum_pr<<endl;
//     // cout<<bs.cum_pr<<endl;
//     // cout<<"\n------------------------"<<endl;

    
//     // cout<<G2_pos<<","<<G2<<endl;

//     // if(G2_pos >= params->max_num and G2 >= params->max_num){
//     //     BS.insert(BS.begin()+pos+1,bs);
//     //     return;
//     // }
//     switch (params->enumbs_add_strategy){
//         // case 1: // Remain all equal strategies
//         //     while( G2_pos > G2 && G_pos > G){
//         //     // while((G2_pos == G2 && G_pos > G) ||  (G2_pos > G2 && G_pos >= G)){
        
//         //         BS.erase(BS.begin()+pos);
//         //         pos -= 1;
//         //         if(pos == -1)
//         //             break;

//         //         // G2_pos = round(pow(2,get<2>(BS[pos].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
//         //         // G_pos = round(pow(2,BS[pos].cum_avg_GB_BKZ.first)*params->enumbs_G_prec)/params->enumbs_G_prec;
//         //         // G_pos = BS[pos].cum_avg_GB_BKZ.first;

//         //         G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//         //         G_pos = round(BS[pos].cum_avg_GB_BKZ.first*params->enumbs_G_prec)/params->enumbs_G_prec;

//         //         slope_pos = BS[pos].slope;
//         //         slope = bs.slope;

//         //     }
            

//         //     if( pos == -1 || (G2_pos >= G2 && G_pos <= G) ){
//         //         BS.insert(BS.begin()+pos+1,bs);
//         //     }

//         //     break;

//         // case 2: //Delete equal strategies
//         //     //(abs(G2_pos - params->max_num) <0.01)  &&
//         //     while(  (( (G2_pos == G2 or (slope - slope_pos)<1e-6) && G_pos > G) ||  ( (G2_pos > G2 or (slope - slope_pos) > 1e-4)&& G_pos >= G))){
        
//         //         BS.erase(BS.begin()+pos);
//         //         pos -= 1;
//         //         if(pos == -1)
//         //             break;

//         //         G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//         //         G_pos = round(BS[pos].cum_avg_GB_BKZ.first*params->enumbs_G_prec)/params->enumbs_G_prec;

//         //     }
//         //     // cout<<"G_pos = "<<G_pos<<endl;
//         //     // cout<<"G = "<<G<<endl;
//         //     //||  (abs(G2_pos - params->max_num) <0.01)
//         //     if( pos == -1 || ((G2_pos > G2 or (slope - slope_pos) > 1e-6) && G_pos <= G) ){
//         //         BS.insert(BS.begin()+pos+1,bs);
//         //     }

//         //     break;

//         case 2: //Delete equal strategies, except small strategies
//             //While cost of previous strategy is worse than new strategy && the range of strategy to erase is smaller than || equal to new strategy.
//             //(abs(G2_pos - params->max_num) <0.01) && 
//             //while(G2_pos >= G2 + params->enumbs_G_prec && G <= G_pos +  params->enumbs_G_prec && compare_max_strategy(BS[pos].S, bs.S)){
            
//             while(G2_pos > G2 && G < G_pos){
            
//                 BS.erase(BS.begin()+pos);
//                 pos -= 1;
//                 if(pos == -1)
//                     break;

//                 G2_pos = get<2>(BS[pos].dsvp_t);
//                 G_pos = BS[pos].cum_avg_GB_BKZ.first;

//             }
            
//             BS.insert(BS.begin()+pos+1,bs);

//             //&& BS[pos+2].cum_pr >= BS[pos+1].cum_pr
//             while (pos+1<BS.size()-1 &&  get<2>(BS[pos+1].dsvp_t) >= get<2>(BS[pos+2].dsvp_t) + params->enumbs_G_prec  && BS[pos+2].cum_avg_GB_BKZ.first <= BS[pos+1].cum_avg_GB_BKZ.first +  params->enumbs_G_prec && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){
//                 BS.erase(BS.begin()+pos+1);
//             }

//             break;
//         default:
//             cerr << "Please choose right add strategy for EnuBS: (1/2/3)" << endl;

//     }
 
//     if(pos <= k)
//         k = -1;

//     if(params->debug){
//         vector<double> G2s = extract_G2();
//         vector<double> sorted_G2s = G2s;
//         // print_vector(G2s);
//         sort(sorted_G2s.rbegin(),sorted_G2s.rend());
//         assert(G2s == sorted_G2s);
//         // assert(no_repeated_value_verification(G2s));
//     }
// }





void EnumBS::BS_add(EnumBS::blocksize_strategy bs, int k){
    
    // params->enumbs_G_prec = PREC;
    // while(pow(2,get<2>(bs.dsvp_t)) * params->enumbs_G_prec < params->enumbs_bound || pow(2,bs.cum_avg_GB_BKZ.first)*params->enumbs_G_prec < params->enumbs_bound ){
    //     // print_strategy(bs.S);
    //     // cout<<pow(2,get<2>(bs.dsvp_t))<<","<<pow(2,bs.cum_avg_GB_BKZ.first)<<endl;
    //     params->enumbs_G_prec *= 10;
    // }
    

    // double G2 = round(pow(2,get<2>(bs.dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
    // double G = round(pow(2,bs.cum_avg_GB_BKZ.first)*params->enumbs_G_prec)/params->enumbs_G_prec;
    // double G = bs.cum_avg_GB_BKZ.first;


    // double slope = bs.slope;
    // double G = bs.cum_avg_GB_BKZ.first;


    //BS.size() == 0, add bs directly
    if(BS.size() == 0){
        // BS.resize(1);
        // BS[0] = bs;
        BS.insert(BS.end(),bs);
        return;
    }

    // int pos = binary_search_for_G2_slope_cum_pr(bs); //binary_search_for_slope(slope);   
    // int pos = binary_search_for_slope(bs.slope);
    int pos = binary_search_for_G(bs.GB.first);



    // and bs.S[0].beta<100
    // if(k > 1000 ){
    //     cout<<"========================="<<endl;
    //     // print_strategy(bs.S);
    //     print_bs(BS[k]);
    //     print_bs(bs);
    //     if(pos != 0)
    //         print_bs(BS[pos-1]);
    //     // cout<<get<2>(bs.dsvp_t)<<endl;
    //     // cout<<get<0>(bs.dsvp_t)<<endl;
    //     // vector<double> slopes = extract_slope();
    //     // print_vector(slopes);
    //     usleep(100000);;;
    //     cout<<"pos = "<<pos<<endl;
    //     // printf("G2 = %e ", G2);
    //     // cout<<params->max_num<<endl;

    //     cout<<"\n------------------------"<<endl;
    // }

    // if(params -> debug)
        // assert( get<2>(bs.dsvp_t) < params->max_num);


    
    //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
    if(pos == 0){
        return;
    }

    

    //BS.size() > 0, && it exits some G2s in EnumBS >= G2_, then pos > 0.
    //G2_tmp == G2
    pos--;
    // double G2_pos = round(pow(2,get<2>(BS[pos].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
    // double G_pos = round(pow(2,BS[pos].cum_avg_GB_BKZ.first)*params->enumbs_G_prec)/params->enumbs_G_prec;
    // double G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;



    // double G_pos = BS[pos].cum_avg_GB_BKZ.first;
    // double slope_pos = BS[pos].slope;
    // double cum_pr_pos = BS[pos].cum_pr;
    // double cum_pr = bs.cum_pr;

    // double G_pos = BS[pos].cum_avg_GB_BKZ.first;

    // cout<<"\n========================"<<endl;
    // printf("G2_pos = %e, G2 = %e ", G2_pos, G2);
    // printf("G_pos = %e, G = %e ", G_pos, G);
    // print_strategy(BS[pos].S);
    // print_strategy(bs.S);
    // printf("\n cum_pr_pos = %e, cum_pr = %e \n", BS[pos].cum_pr, bs.cum_pr);
    // cout<<BS[pos].cum_pr<<endl;
    // cout<<bs.cum_pr<<endl;
    // cout<<"\n------------------------"<<endl;

    
    // cout<<G2_pos<<","<<G2<<endl;

    // if(G2_pos >= params->max_num and G2 >= params->max_num){
    //     BS.insert(BS.begin()+pos+1,bs);
    //     return;
    // }
    
        // case 1: // Remain all equal strategies
        //     while( G2_pos > G2 && G_pos > G){
        //     // while((G2_pos == G2 && G_pos > G) ||  (G2_pos > G2 && G_pos >= G)){
        
        //         BS.erase(BS.begin()+pos);
        //         pos -= 1;
        //         if(pos == -1)
        //             break;

        //         // G2_pos = round(pow(2,get<2>(BS[pos].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
        //         // G_pos = round(pow(2,BS[pos].cum_avg_GB_BKZ.first)*params->enumbs_G_prec)/params->enumbs_G_prec;
        //         // G_pos = BS[pos].cum_avg_GB_BKZ.first;

        //         G2_pos = round(get<2>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
        //         G_pos = round(BS[pos].cum_avg_GB_BKZ.first*params->enumbs_G_prec)/params->enumbs_G_prec;

        //         slope_pos = BS[pos].slope;
        //         slope = bs.slope;

        //     }
            

        //     if( pos == -1 || (G2_pos >= G2 && G_pos <= G) ){
        //         BS.insert(BS.begin()+pos+1,bs);
        //     }

        //     break;

            //(abs(G2_pos - params->max_num) <0.01)  &&
            //&& G_pos > 1
            // if(slope >= slope_pos && cum_pr < cum_pr_pos){
            //     cout<<"================="<<endl;
            //     print_bs(bs);
            //     print_bs(BS[pos]);
            //     cout<<"----------------"<<endl;
            //     usleep(100000);;;
                
            // }
            // && cum_pr >= cum_pr_pos
            // if(BS.size()>10000){
            //     print_BS(BS);
            //     throw "";
            // }
            // if(params->debug){

            
            //     if(pos > -1 && (get<2>(bs.dsvp_t) <  get<2>(BS[pos].dsvp_t)  || (get<2>(bs.dsvp_t) == get<2>(BS[pos].dsvp_t) && bs.slope >= BS[pos].slope && bs.cum_pr >= BS[pos].cum_pr)) && bs.cum_avg_GB_BKZ.first < BS[pos].cum_avg_GB_BKZ.first && compare_max_strategy(BS[pos].S, bs.S))
            //     {
            //         cout<<"========"<<endl;
            //         print_bs(bs);
            //         print_bs(BS[pos]);
            //         cout<<"......."<<endl;
            //         usleep(100000);;;
            //     }
            // }


            

            // if(params->debug){
            //     if(bs.S.size()>10){
            //         //bs.slope < BS[pos].slope - params->enumbs_slope_prec ||
            //         if((pos > -1 && (get<2>(bs.dsvp_t) <  get<2>(BS[pos].dsvp_t) + params->enumbs_G_prec && ( bs.cum_pr < BS[pos].cum_pr))) && bs.cum_avg_GB_BKZ.first <= BS[pos].cum_avg_GB_BKZ.first + params->enumbs_G_prec  && compare_max_strategy(BS[pos].S, bs.S)){
            //             cout<<"====Add===="<<endl;
            //             printf("pos=%d\n",pos+1);
            //             print_bs(bs);
            //             print_bs(BS[pos]);
            //             cout<<"......."<<endl;
            //             usleep(100000);
            //         }
            //     }
            // }
            //Whether should we delete the strategy BS[pos]
            // while(pos > -1 && (get<2>(bs.dsvp_t) <=  get<2>(BS[pos].dsvp_t) + params->enumbs_G_prec && bs.slope >= BS[pos].slope - params->enumbs_slope_prec && bs.cum_pr >= BS[pos].cum_pr - params->enumbs_cumpr_prec) && bs.cum_avg_GB_BKZ.first <= BS[pos].cum_avg_GB_BKZ.first + params->enumbs_G_prec  && compare_max_strategy(BS[pos].S, bs.S)){
            
            


            // while(pos > -1 && bs.slope >= BS[pos].slope - params->enumbs_slope_prec && bs.cum_avg_GB_BKZ.first <= BS[pos].cum_avg_GB_BKZ.first + params->enumbs_G_prec && compare_max_strategy(BS[pos].S, bs.S)){

        while(pos > -1 && ( bs.GB.first < BS[pos].GB.first + params->enumbs_G_prec || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && bs.GB.first == BS[pos].GB.first + params->enumbs_G_prec) ) && compare_max_strategy(BS[pos].S, bs.S)){

            if(params->debug){
                if(bs.GB.second < params->max_RAM){
                    cout<<"====Delete1===="<<endl;
                    printf("pos=%d\n",pos);
                    print_bs(bs);
                    print_bs(BS[pos]);
                    cout<<"......."<<endl;
                    usleep(10000000);
                }
            }

            BS.erase(BS.begin()+pos);

            pos--;
        }

        BS.insert(BS.begin()+pos+1,bs);
            
        // if( pos == -1 || not ( cum_pr <= cum_pr_pos && G >= G_pos && compare_max_strategy(bs.S,BS[pos].S)) ){

        // if(params->debug){
        //     //
        //     if(get<2>(bs.dsvp_t) == get<2>(BS[pos].dsvp_t) && slope == slope_pos && bs.cum_pr == BS[pos].cum_pr && G >= G_pos && compare_max_strategy(bs.S, BS[pos].S)){
        //         cout<<"========"<<endl;
        //         print_bs(bs);
        //         print_bs(BS[pos]);
        //         cout<<"......."<<endl;
        //         usleep(100000);;;
        //     }
        // }
        // if(pos == -1 or not (get<2>(bs.dsvp_t) == get<2>(BS[pos].dsvp_t) && bs.slope == BS[pos].slope && bs.cum_pr == BS[pos].cum_pr && bs.cum_avg_GB_BKZ.first >= BS[pos].cum_avg_GB_BKZ.first && compare_max_strategy(bs.S, BS[pos].S))){
        //     BS.insert(BS.begin()+pos+1,bs);

            

        // while( pos+1<BS.size()-1 && (get<2>(BS[pos+2].dsvp_t) <=  get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec  && BS[pos+2].slope >= BS[pos+1].slope - params->enumbs_slope_prec && BS[pos+2].cum_pr >= BS[pos+1].cum_pr - params->enumbs_cumpr_prec) && BS[pos+2].cum_avg_GB_BKZ.first <= BS[pos+1].cum_avg_GB_BKZ.first + params->enumbs_G_prec && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){
        // while( pos+1<BS.size()-1 && BS[pos+2].slope >= BS[pos+1].slope - params->enumbs_slope_prec &&  BS[pos+2].cum_avg_GB_BKZ.first <= BS[pos+1].cum_avg_GB_BKZ.first + params->enumbs_G_prec && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){

        while( pos+1<int(BS.size())-1 && ( BS[pos+2].GB.first < BS[pos+1].GB.first + params->enumbs_G_prec || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec &&  BS[pos+2].GB.first == BS[pos+1].GB.first + params->enumbs_G_prec )) && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){

            if(params->debug){
                if(bs.GB.second < params->max_RAM){
                    cout<<"====Delete2===="<<endl;
                    printf("pos=%d\n",pos+1);
                    print_bs(BS[pos+1]);
                    print_bs(BS[pos+2]);
                    cout<<"......."<<endl;
                    usleep(10000000);
                }
            }

            // while (pos+1<BS.size()-1 && BS[pos+2].slope >= BS[pos+1].slope - params->enumbs_slope_prec  && BS[pos+2].cum_avg_GB_BKZ.first <= BS[pos+1].cum_avg_GB_BKZ.first +  params->enumbs_G_prec && compare_max_strategy(BS[pos+1].S, BS[pos+2].S)){
            BS.erase(BS.begin()+pos+1);
            pos--;
        }
                
           
    // if(params->debug){
    //     cout<<"========"<<endl;
    //     print_bs(bs);
    //     print_BS(BS);
    //     cout<<"......."<<endl;
    //     usleep(30000000);
    // }


 
    if(pos <= k)
        k = -1;

    // if(params->debug){
        // vector<tuple<double,double,double>> basis_quality_list = extract_G2_slope_cum_pr();
        // print_vector(basis_quality_list);

        // vector<double> G2s = extract_G2();
        // vector<double> sorted_G2s = G2s;
        // sort(sorted_G2s.rbegin(),sorted_G2s.rend());
        // assert(G2s == sorted_G2s);

        // vector<double> slopes = extract_slope();
        // vector<double> sorted_slopes = slopes;
        // sort(sorted_slopes.begin(),sorted_slopes.end());
        // assert(slopes == sorted_slopes);

        // vector<double> cum_prs = extract_cum_pr();
        // vector<double> sorted_cum_prs = cum_prs;
        // sort(sorted_cum_prs.begin(), sorted_cum_prs.end());
        // assert(cum_prs == sorted_cum_prs);
        // assert(no_repeated_value_verification(cum_prs));
    // }
}




void EnumBS::BS_add_op(EnumBS::blocksize_strategy bs, int k){
    
    if(BS.size() == 0 or bs.GB.first <= BS[k].GB.first){
        BS.insert(BS.end(),bs);
     
    }
    else{
        if(params->debug){
            if(bs.GB.second < params->max_RAM){
                cout<<"====Add===="<<endl;
                printf("k=%d\n",k);
                print_bs(bs);
                print_bs(BS[k]);
                cout<<"......."<<endl;
                usleep(10000000);
            }
        }
    }
}




// void EnumBS::BS_add_cdsvp(EnumBS::blocksize_strategy bs, int k){
//     double cdsvp = round(get<0>(bs.dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//     double G = bs.cum_avg_GB_BKZ.first;

//     //BS.size() == 0, add bs directly
//     if(BS.size() == 0){
//         // BS.resize(1);
//         // BS[0] = bs;
//         BS.insert(BS.end(),bs);
//         return;
//     }



//     // int pos = find_pos_for_dsvp(cdsvp);
//     int pos = binary_search_for_cdsvp(cdsvp);

    
//     //BS.size() > 0, but all dsvps in EnumBS are smaller than dsvp_, don't add dsvp_, then pos = 0.
//     if(pos == 0){
//         return;
//     }

//     //BS.size() > 0, && it exits some dsvps in EnumBS >= dsvp_, then pos > 0.
//     //cdsvp_tmp == cdvsp
//     pos--;
//     double cdsvp_pos = round(get<0>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//     double G_pos = BS[pos].cum_avg_GB_BKZ.first;

//     // printf("\n%e, %e\n",cdsvp_pos, cdsvp);
//     // printf("\n%e, %e\n",G_pos, G);
    
//     // if(bs.S[bs.S.size()-1].beta<80){
//     //     printf("\n\n====================================\n");
//     //     printf("pos = %d, BS_size = %d", pos, int(BS.size()));
//     //     printf("\n%e, %e\n",cdsvp_pos, cdsvp);
//     //     printf("\n%e, %e\n",G_pos, G);
//     //     print_strategy(BS[pos].S);
//     //     print_strategy(bs.S);
//     //     printf("...........................\n");
//     // }
//     //
//     //
//     while( ((cdsvp_pos == cdsvp && G_pos > G) || (cdsvp_pos > cdsvp && G_pos >= G)) && compare_max_strategy(BS[pos].S, bs.S)){

//         if( k >= pos){
//             // break;
//             k = 0;
//         }
//         BS.erase(BS.begin()+pos);
    
//         pos -= 1;
        
//         if(pos == -1)
//             break;
//         cdsvp_pos = get<0>(BS[pos].dsvp_t);//round(get<0>(BS[pos].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec;
//         G_pos = BS[pos].cum_avg_GB_BKZ.first;
//     }

    

//     // print_bs(bs);
//     if(pos == -1 || (cdsvp_pos > cdsvp && G_pos <= G)){
//     // if(pos == -1 || (cdsvp_pos >= cdsvp && G_pos <= G)){
//     // if(pos == -1 || (cdsvp_pos > cdsvp && G_pos <= G) || cdsvp_pos == cdsvp){
//         BS.insert(BS.begin()+pos+1,bs);
//     }
    

//     //Verification
//     if(params->debug){
//         vector<double> cdsvps = extract_cdsvp();
//         vector<double> sorted_cdsvps = cdsvps;
//         sort(sorted_cdsvps.rbegin(),sorted_cdsvps.rend());
//         print_vector(cdsvps);
//         assert(cdsvps == sorted_cdsvps);
//         assert(no_repeated_value_verification(cdsvps));


//         // vector<double> G2s = extract_G2();
//         // vector<double> sorted_G2s = G2s;
//         // print_vector(G2s);
//         // sort(sorted_G2s.rbegin(),sorted_G2s.rend());
//         // assert(G2s == sorted_G2s);
        
//     }


// }


bool EnumBS::pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB_BKZ, pair<double,double> &cum_avg_GB_BKZ, pair<double,double> &GB, double &cum_pr, int beta, int jump, tuple<double,int,double,double> &dsvp_t_, double &slope){

    
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

    vector<double> l_;

    sim -> simulate(l_,l,beta_,jump,1);

   
    if(l_[d-beta_] == l[d-beta_]){
        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model, params->progressive_sieve, params->worst_case);
        slope = get_current_slope(l, 0, d);
        return false;
    }
    else{
        l = l_;
        slope = get_current_slope(l, 0, d);
        boost::math::chi_squared chisquare(beta_);
        double pr = boost::math::cdf(chisquare,pow(2,2.*l[d-beta_]));
        

        pair<double,double> GB_BKZ = cost->bkz_cost(d,beta,jump,params->cost_model);
        // if(beta == 178)
        //     printf("beta = %d, l[i]= %e, G = %e, rem_pr = %e, pr = %e\n", beta, pow(2,2.*l[d-beta]), GB.first, rem_pr, pr);


        cum_GB_BKZ.first = log2(pow(2,cum_GB_BKZ.first)+pow(2,GB_BKZ.first));
        if(not params->worst_case){
            cum_avg_GB_BKZ.first = log2(pow(2,cum_avg_GB_BKZ.first)+(pow(2,cum_GB_BKZ.first)*rem_pr*pr));
        }
        else{
            cum_avg_GB_BKZ.first = cum_GB_BKZ.first;
        }
        cum_avg_GB_BKZ.second = max(cum_avg_GB_BKZ.second, GB_BKZ.second);

        // printf("cum_G = %e\n", cum_avg_GB_BKZ.first );
        cum_pr += rem_pr * pr;
        rem_pr = 1. - cum_pr;


        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model, params->progressive_sieve, params->worst_case);

        GB.first = log2(pow(2,cum_avg_GB_BKZ.first) + pow(2,get<2>(dsvp_t_)));
        GB.second = max(cum_avg_GB_BKZ.second, get<3>(dsvp_t_));

        return true;
    }

}


void EnumBS::max_tour_for_pnjbkz_beta(int k, int beta,int jump){
    EnumBS::blocksize_strategy bs = BS[k];
    tuple<double,int,double,double> dsvp_t1;
    double G21;
    // double G20 = get<2>(bs.dsvp_t), G21;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;

    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ=bs.cum_avg_GB_BKZ, GB = bs.GB, min_GB = bs.GB;

    int  loop = 0;
    
    double slope0 = bs.slope, slope1;
    bool sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, beta, jump, dsvp_t1,slope1);

    if(not sim_term)
        return;
    
    G21 = get<2>(dsvp_t1);

    assert(G21 >= 0.);

    bool leaf = false;
    while(((slope1 - slope0) > params->enumbs_slope_prec) && loop < params->max_loop ){

        loop +=1;
        slope0 = slope1;
        
        if(loop ==1){
            S.insert(S.end(),{beta, jump, loop});
        }
        else
            // S[len_S-1].tours = loop;
            S[S.size()-1].tours = loop;
        
        if(min_GB.first > GB.first or min_GB.second > params->max_RAM){
            min_GB = GB;
            leaf = false;
        }
        else if(GB.first >= min_GB.first + min_G_prec)
            leaf = true;

        bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr,slope1,min_GB};

        if(params->verification){
            pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
            //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
            //cerr<<"cum_G="<<cum_avg_GB_BKZ.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
            assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
            assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
        }
    
        // EnumBS::BS_add_G2(bs, k);
        // EnumBS::BS_add_slope(bs,k);
        // BS_add_op(bs,k);
        if(params->enum_add_G2)
            BS_add_G2(bs, k);
        else
            BS_add(bs,k);

        if(cum_pr >= 0.999 or not sim_term or leaf)
            break;
    
        sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, beta, jump, dsvp_t1,slope1);
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
//     pair<double,double> cum_avg_GB_BKZ=bs.cum_avg_GB_BKZ;

//     int  loop = 0;

//     dsvp_t1 = pnjbkz_beta_loop(l, cum_avg_GB_BKZ, cum_pr, beta, jump);
    
//     // G21 = get<2>(dsvp_t1);
//     G21 = get<2>(dsvp_t1);

//     // while(G20 - G21 >= 1){
//     while(G20 > G21 && loop < params->max_loop && cum_pr <= 0.999 ){
//         // if(loop > 5)
//         //     printf("\n%d, %d, %d, %e, %e\n", beta, jump, loop, G20,G21);
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
//         bs = {dsvp_t1, S, l, cum_avg_GB_BKZ, cum_pr};

//         if(params->verification){
//             pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
//             //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
//             //cerr<<"cum_G="<<cum_avg_GB_BKZ.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
//             assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
//             assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
//         }
    
//         // k_flag = EnumBS::BS_add(bs, k);
//         EnumBS::BS_add_G2(bs, k);

//         dsvp_t1 = pnjbkz_beta_loop( l, cum_avg_GB_BKZ, cum_pr, beta, jump);
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
//     pair<double,double> cum_avg_GB_BKZ=bs.cum_avg_GB_BKZ;

//     int  loop = 0;

//     dsvp_t1 = pnjbkz_beta_loop(l, cum_avg_GB_BKZ, cum_pr, beta, jump);
    
//     // G21 = get<2>(dsvp_t1);
//     G21 = get<2>(dsvp_t1);

//     // while(G20 - G21 >= 1){
//     while(G20 - G21 >= 1){
//         // if(loop > 5)
//         //     printf("\n%d, %d, %d, %e, %e\n", beta, jump, loop, G20,G21);
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
//         bs = {dsvp_t1, S, l, cum_avg_GB_BKZ, cum_pr};

//         if(params->verification){
//             pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
//             //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
//             //cerr<<"cum_G="<<cum_avg_GB_BKZ.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
//             assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
//             assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
//         }
    
//         // k_flag = EnumBS::BS_add(bs, k);
//         // k_flag = EnumBS::BS_add_G2(bs, k);
//         EnumBS::BS_add_G2(bs, k);

//         dsvp_t1 = pnjbkz_beta_loop( l, cum_avg_GB_BKZ, cum_pr, beta, jump);
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
//         pair<double,double> cum_avg_GB_BKZ=bs.cum_avg_GB_BKZ;

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
        
//         dsvp_t1 = pnjbkz_beta_loop(l, cum_avg_GB_BKZ, cum_pr, beta, jump);
        
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
                
//             bs = {dsvp_t1, S, l, cum_avg_GB_BKZ, cum_pr};
            
//             if(params->verification){
//                 pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
//                 assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
//                 assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
//             }
        
//             tmpBS[index].insert(tmpBS[index].end(),bs);

//             // if(!BS_add_determine(bs, k)){
//             //     flag = false;
//             // }


//             dsvp_t1 = pnjbkz_beta_loop( l, cum_avg_GB_BKZ, cum_pr, beta, jump);
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
        // double G20 = get<2>(bs.dsvp_t), G21;
        vector<double> l = bs.l; 
        vector<strategy> S = bs.S;

        double cum_pr = bs.cum_pr;
        pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ = bs.cum_avg_GB_BKZ, GB = bs.GB, min_GB = bs.GB;

        int beta = beta_j_tid[i].first, jump = beta_j_tid[i].second;
    
        int f = dims4free(beta);
        if( (f == 0 or beta < 79) && jump > 1)
            continue;
        if(f!=0 && jump >= f)
            continue;

        int index = beta_j_t_id_begin + i;
        
        int  loop = 0;


        double slope0 = bs.slope, slope1;


        bool sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, beta, jump, dsvp_t1, slope1);

        if(not sim_term)
            continue;
        
        // G21 = get<2>(dsvp_t1);

        // assert(G21 >= 0.);

        tmpBS[index].clear();

        // while((abs(G20 - params->max_num)<0.01 or G20 > G21 or (slope1 - slope0) > 1e-6) && loop < params->max_loop){

        bool leaf = false;
        while( (slope1 - slope0) > params->enumbs_slope_prec && loop < params->max_loop){
            loop +=1;
            // G20 = G21;
            slope0 = slope1;
            if(loop ==1){
                S.insert(S.end(),{beta, jump, loop});
            }
            else
                S[S.size()-1].tours = loop;
                
            if(min_GB.first > GB.first or min_GB.second > params->max_RAM){
                min_GB = GB;
                leaf = false;
            }
            else if(GB.first >= min_GB.first + min_G_prec)
                leaf = true;
            
            
            // if(params->debug){
            //     if(bs.GB.second < params->max_RAM){
            //         print("----");
            //         print_bs(bs);
            //         print("====");
            //     }
            // }
    
            bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, slope1, min_GB};
            
            if(params->verification){
                pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
                assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
                assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
            }


            
        
            int len_tmpBS = tmpBS[index].size();


            // if(params->debug){
            //     cout<<"====Delete loop===="<<endl;
            //     print_bs(bs);
            //     print_bs(tmpBS[index][len_tmpBS - 1]);
            //     cout<<"......."<<endl;
            //     usleep(10000000);
            // }

        
            if( loop > 1 && ((get<2>(bs.dsvp_t) < get<2>(tmpBS[index][len_tmpBS - 1].dsvp_t) + params->enumbs_G_prec  || (get<2>(bs.dsvp_t) == get<2>(tmpBS[index][len_tmpBS - 1].dsvp_t) + params->enumbs_G_prec && bs.slope > tmpBS[index][len_tmpBS - 1].slope - params->enumbs_slope_prec ))) && bs.cum_avg_GB_BKZ.first < tmpBS[index][len_tmpBS - 1].cum_avg_GB_BKZ.first + params->enumbs_G_prec){
                tmpBS[index][len_tmpBS-1] = bs;   
            }
            else
                tmpBS[index].insert(tmpBS[index].end(),bs);
                    
            // tmpBS[index].insert(tmpBS[index].end(),bs);
            


            // BS_add_op(bs, k);
            // if(bs.cum_avg_GB_BKZ.first <= BS[k].cum_avg_GB_BKZ.first)
            //     tmpBS[index].insert(tmpBS[index].end(),bs);
            
            // if(bs.GB.first <= BS[k].GB.first){
            //     if( not flag && ( bs.cum_avg_GB_BKZ.first < tmpBS[index][len_tmpBS - 1].cum_avg_GB_BKZ.first || (bs.slope > tmpBS[index][len_tmpBS - 1].slope && bs.cum_avg_GB_BKZ.first == tmpBS[index][len_tmpBS - 1].cum_avg_GB_BKZ.first))){
            //         tmpBS[index][len_tmpBS-1] = bs;
            //     }
            //     else{
            //         tmpBS[index].insert(tmpBS[index].end(),bs);
            //         flag = false;
            //     }
            // }


            if(cum_pr >= 0.999 or not sim_term or leaf)
                break;
            
            sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, beta, jump, dsvp_t1, slope1);
        
            // G21 = get<2>(dsvp_t1);
            // assert(G21 >= 0.);
        }

        // if(!flag)
        //     return flag;
        
    }
    // return flag;
}


void EnumBS::enumbs_est(vector<double> l0){
    /*
    input: l -- gs-lengths; 
    Return: Optimal strategy to minimize solving cost.
    */
    int beta_start = params->beta_start, k = 0, d = l0.size(),j_start,len_S;
    
    tuple<double,int,double,double>  dsvp0_t = dsvp_predict(l0, 0., cost, params->cost_model, params->progressive_sieve, params->worst_case);
    blocksize_strategy bs =  {dsvp0_t, {},l0,make_pair(0.,0.), make_pair(0.,0.), make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t)), 0., get_current_slope(l0,0,d), make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t))};


    BS.insert(BS.end(),bs);

    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

        //Add a normal two-step strategy from (beta_start,1,1) to (d,1,1)


    
    tuple<double,int,double,double> dsvp_t1;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ = bs.cum_avg_GB_BKZ, GB = bs.GB, min_GB = bs.min_GB;
    double slope1;
    bool leaf = false;
    
    for(int beta = beta_start;  beta < d; beta++){
        if(cum_pr >= 0.999 or leaf)
            break;
        pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, beta, 1, dsvp_t1, slope1);
        S.insert(S.end(),{beta,1,1}); 
        if(min_GB.first > GB.first or min_GB.second > params->max_RAM){
            min_GB = GB;
            leaf = false;
        }
        else if(GB.first >= min_GB.first + min_G_prec)
            leaf = true;
    
        bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, slope1, min_GB};

        if(params->enum_add_G2)
            BS_add_G2(bs, k);
        else
            BS_add(bs,k);
        
        k++;
    }

    
    while( k < int(BS.size())){
        bs = BS[k];

        len_S = bs.S.size();

        if(bs.cum_pr >= 0.999){
            k++;
            continue;
        }
    
        if(len_S == 0){
            beta_start = params->beta_start;
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
            printf("\r index: %8d, (beta,j): (%4d,%4d) --> (%4d,%4d), goal index: %8d, cost = " ,k+1,beta_start,j_start,(min(params->max_dim,d)-1-beta_start)/params->gap*params->gap+beta_start,1,int(BS.size()));
            cout<<setprecision(2)<<diff.count()<<'s';
        }

        k ++;
    }
    printf("\n");
    //print_BS(BS);

    //Find the optimal strategy
    double Gmin = params->max_num, Bmin  = params->max_num;
    EnumBS::blocksize_strategy bsmin;
    for(int i = 0; i<int(BS.size()); i++){
        // G1 = BS[i].cum_avg_GB_BKZ.first;
        // G2 = get<2>(BS[i].dsvp_t);
        // G = log2(pow(2,G1)+pow(2,G2));
        // B = max(get<3>(BS[i].dsvp_t),BS[i].cum_avg_GB_BKZ.second);
        
        if(BS[i].GB.first<Gmin and BS[i].GB.second <= params->max_RAM){
            bsmin = BS[i];
            Gmin = BS[i].GB.first;
            Bmin = BS[i].GB.second;
        }
    }
    printf("Find the optimal Strategy through EumBS!!\n");
    print_bs(bsmin);
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log2(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec) = %3.2f h, Memory Cost = %3.2f log2(bit) = %3.2f TB \n", Gmin, pow(2,Gmin)/3600, Bmin, pow(2,Bmin-43));
}

//return value: to determine whether the current bs0 will be changed || not.
//False: bs0 will be changed.
//True: bs0 will not change.
// bool EnumBS::BS_add_determine(EnumBS::blocksize_strategy bs, int k){
//     //int cdsvp = ceil( get<2>(bs.dsvp_t));
//     double dsvp =  get<2>(bs.dsvp_t);
//     double G = bs.cum_avg_GB_BKZ.first;

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
    
//     double G_pos = BS[pos].cum_avg_GB_BKZ.first;
    

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
//         G_pos = BS[pos].cum_avg_GB_BKZ.first;
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



void EnumBS::enumbs_est_in_parallel(vector<double> l0){
    /*
    input: l0 -- gs-lengths;
    Return: Optimal strategy to minimize solving cost.
    */
    set_threads(params->threads);
    
 
    int k = 0, d = l0.size(),beta, beta_start = params->beta_start;
    blocksize_strategy bs;

    tuple<double,int,double,double>  dsvp0_t = dsvp_predict(l0, 0., cost,params->cost_model, params->progressive_sieve, params->worst_case);

    bs = {dsvp0_t, {},l0,make_pair(0.,0.), make_pair(0.,0.),make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t)), 0., get_current_slope(l0,0,d), make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t))};

    BS.insert(BS.end(),bs);
    


    //Add a normal two-step strategy from (beta_start,1,1) to (d,1,1)


    
    tuple<double,int,double,double> dsvp_t1;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ = bs.cum_avg_GB_BKZ, GB = bs.GB, min_GB = bs.min_GB;
    double slope1;
    bool leaf = false;
    
    for(int beta = beta_start;  beta < d; beta++){
        if(cum_pr >= 0.999 or leaf)
            break;
        pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, beta, 1, dsvp_t1, slope1);
        S.insert(S.end(),{beta,1,1}); 
        if(min_GB.first > GB.first or min_GB.second > params->max_RAM){
            min_GB = GB;
            leaf = false;
        }
        else if(GB.first >= min_GB.first + min_G_prec )
            leaf = true;
    
        bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, GB, cum_pr, slope1, min_GB};

        if(params->enum_add_G2)
            BS_add_G2(bs, k);
        else
            BS_add(bs,k);
        
        k++;
    }

    
    int j,len_S;
    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

    k = 0;
    while( k < int(BS.size())){
        bs = BS[k];
        len_S = bs.S.size();
        if(params->debug){
            print_bs(bs);
        }
        if(bs.cum_pr >= 0.999){
            k++;
            continue;
        }
      
    
        if(len_S == 0){
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
            
                    // cout<<beta_j[t_id][0].first<<", "<<beta_j[t_id][0].second<<endl;
                    // cerr<<"t_id = "<<t_id<<", beta_j[i].size() = "<<beta_j[t_id].size()<<endl;
                }
                // cerr<<len<<","<<Sum<<endl;
                assert(len == Sum);
            }
        
            
            auto start= system_clock::now();
            for (t_id = 0; t_id < threads; t_id++){
                threadpool.push([this, t_id, beta_j_t_id_begins, beta_j, k](){
                    max_tour_for_pnjbkz_beta_in_parallel(  beta_j_t_id_begins[t_id], beta_j[t_id], k); 
                });            
            }
            threadpool.wait_work(); 
            
            if(params->verbose){
                auto finish = system_clock::now();
                duration<double> diff = finish - start;
                printf("\r index: %8d, (beta,j): (%4d,%4d) --> (%4d,%4d), goal index: %8d, cost =" ,k+1,beta_j[0][0].first,beta_j[0][0].second,beta_j[threads-1][beta_j[threads-1].size()-1].first,beta_j[threads-1][beta_j[threads-1].size()-1].second,int(BS.size()));
                cerr<<setprecision(2)<<diff.count()<<'s';
            }
           
            
            for(int i = 0; i< len; i++){
                for(int ii = 0; ii < int(tmpBS[i].size()); ii++){
                    if(params->enum_add_G2)
                        BS_add_G2(tmpBS[i][ii], k);
                    else
                        BS_add(tmpBS[i][ii],k);
                    // BS_add(tmpBS[i][ii], k);
                    // BS_add_op(tmpBS[i][ii], k);
                    // BS.insert(BS.end(), tmpBS[i][ii]);
                    
                }
            
            }
            
            tmpBS.clear();
            beta_j.clear();
            
        }
        k++;
    }
    printf("\n");

    //Find the optimal strategy
    double Gmin = params->max_num, Bmin  = params->max_num;
    EnumBS::blocksize_strategy bsmin;
    bool flag = false;
    for(int i = 0; i<int(BS.size()); i++){
        // G1 = BS[i].cum_avg_GB_BKZ.first;
        // G2 = get<2>(BS[i].dsvp_t);
        // G = log2(pow(2,G1)+pow(2,G2));
        // B = max(get<3>(BS[i].dsvp_t),BS[i].cum_avg_GB_BKZ.second);
        // print_bs(BS[i]);
        if(BS[i].GB.first<Gmin  and BS[i].GB.second < params->max_RAM){
            bsmin = BS[i];
            Gmin = BS[i].GB.first;
            Bmin = BS[i].GB.second;
            flag = true;
        }
    }
    if(flag){
        printf("Find the optimal Strategy through EumBS!!\n");
    }
    else{
        for(int i = 0; i<int(BS.size()); i++){
            if(BS[i].GB.second<Bmin){
                bsmin = BS[i];
                Gmin = BS[i].GB.first;
                Bmin = BS[i].GB.second;
            }
        }
        printf("There's no strategy whose memory cost (%lf log(bit)) is below %lf log(bit), the lowest memory cost strategy is:\n", Bmin, params->max_RAM);
    }
    print_bs(bsmin);
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec) = %3.2f h, Memory Cost = %3.2f log2(bit) = %3.2f TB \n", Gmin, pow(2,Gmin)/3600, Bmin, pow(2,Bmin-43));
}



pair<double,double> EnumBS::strategy_verification(vector<double> l,vector<strategy> S){

    int d = l.size();
    double cum_pr = 0., rem_pr = 1., proba, G1cum=0., B1cum = 0.;
    // BKZJSim* sim = new BKZJSim();
    // COST* cost = new COST();
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

    tuple<double,int,double,double> dsvp_t =  dsvp_predict(l, cum_pr, cost,params->cost_model, params->progressive_sieve, params->worst_case);
    double G2 = get<2>(dsvp_t);
    double G = log2(pow(2,G1cum)+pow(2,G2));
    printf("Verified cum_pr = %e \n ", cum_pr);
    printf("Verified G1 = %e, G2 = %e, dsvp = %e\n", G1cum,G2,get<0>(dsvp_t));
    printf("G = %e\n", G );

    return make_pair(G1cum, cum_pr);

}
