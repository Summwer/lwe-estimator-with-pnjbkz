#include "enumbs.h"
#include <unistd.h>

void EnumBS::set_threads(int nr)
{
    assert(nr >= 1);
    threadpool.resize(nr);
}

void EnumBS::print_strategy(vector<EnumBS::strategy> S){
    cout<<"S(beta,jump,tours):[";
    for(int i = 0; i < int(S.size()); i ++){
        printf("(%d,%d,%d)",S[i].beta,S[i].jump,S[i].tours);
        if(i!=int(S.size()) - 1)
            printf(",");
    }
    cout<<"]"<<endl;
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
    int d = int(bs.l.size());
    cout<<"bs = ";
    // int d = int(bs.l.size());
    // double G1 = strategy_verification(l,BS[i].S).first;
    if(params->cost_model == 1)
        printf("(slope = %e, G_BKZ = %e log2(gate), B_BKZ = %e log2(bit), cum-pr = %e, pump-{%d,%d,%d}, G_dsvp = %e log2(gate), B_dsvp = %e bit, avgG = %e log2(gate), avgB = %e log2(bit),  G = %e log2(gate), min_GB.first = %e log2(gate), leaf = %d)\n",  bs.slope, bs.cum_GB_BKZ.first, bs.cum_GB_BKZ.second, bs.cum_pr,  d- get<1>(bs.dsvp_t), get<1>(bs.dsvp_t), get<1>(bs.dsvp_t) - get<0>(bs.dsvp_t), get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), bs.avg_GB.first, bs.avg_GB.second, bs.GB_nopr.first, bs.min_GB.first,bs.leaf);  
    if(params->cost_model == 2)
        printf("(slope = %e, G_BKZ = %e log2(sec), B_BKZ = %e log2(bit), cum-pr = %e, pump-{%d,%d,%d}, G_dsvp = %e log2(sec), B_dsvp = %e bit, avgG = %e log2(sec), avgB = %e log2(bit), G=%e log2(sec), min_GB.first = %e log2(sec), leaf = %d)\n",  bs.slope, bs.cum_GB_BKZ.first, bs.cum_GB_BKZ.second, bs.cum_pr, d- get<1>(bs.dsvp_t), get<1>(bs.dsvp_t), get<1>(bs.dsvp_t) - get<0>(bs.dsvp_t), get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), bs.avg_GB.first, bs.avg_GB.second,  bs.GB_nopr.first, bs.min_GB.first,bs.leaf);  
    print_strategy(bs.S);
}


// void EnumBS::print_bs(blocksize_strategy bs){
//     cout<<"bs = ";

//     // double G1 = strategy_verification(l,BS[i].S).first;
//     if(params->cost_model == 1)
//         printf("(slope = %e, G_BKZ = %e gate, B_BKZ = %e bit cum-pr = %e, dsvp = %e, dsvp_r = %3d, G_dsvp = %e gate, B_dsvp = %e bit, G = %e gate, B = %e bit, min_GB.first = %e gate)\n",  bs.slope, bs.cum_avg_GB_BKZ.first, bs.cum_avg_GB_BKZ.second, bs.cum_pr, get<0>(bs.dsvp_t), get<1>(bs.dsvp_t),get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), log2(pow(2,get<2>(bs.dsvp_t)) + pow(2,bs.cum_avg_GB_BKZ.first)), max(bs.cum_avg_GB_BKZ.second,get<3>(bs.dsvp_t)), bs.min_GB.first);  
//     if(params->cost_model == 2)
//         printf("(slope = %e, G_BKZ = %e sec, B_BKZ = %e bit cum-pr = %e, dsvp = %e, dsvp_r = %3d, G_dsvp = %e sec, B_dsvp = %e bit, G = %e sec, B = %e bit,  min_GB.first = %e gate)\n",  bs.slope, bs.cum_avg_GB_BKZ.first, bs.cum_avg_GB_BKZ.second, bs.cum_pr, get<0>(bs.dsvp_t), get<1>(bs.dsvp_t),get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), log2(pow(2,get<2>(bs.dsvp_t)) + pow(2,bs.cum_avg_GB_BKZ.first)), max(bs.cum_avg_GB_BKZ.second,get<3>(bs.dsvp_t)), bs.min_GB.first);  
//     print_strategy(bs.S);
// }


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
    double G_tmp = BS[len-1].avg_GB.first;
    if(G_tmp >= G)
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        if(BS[mid].avg_GB.first >= G) left = mid + 1;
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



int EnumBS::binary_search_for_G2_slope(blocksize_strategy bs){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    # sort in (G2,slope,cum_pr). Find the position to put (G2,slope,cum_pr), G2-->smaller, slope --> larger, cum_pr --> larger.
    */
    
    
    int len = int(BS.size());
    // double G2_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
 
    //|| (get<2>(BS[len-1].dsvp_t) == get<2>(bs.dsvp_t) && BS[len-1].slope == bs.slope && BS[len-1].cum_pr <= bs.cum_pr)
    if( get<4>(BS[len-1].dsvp_t) > get<4>(bs.dsvp_t) || (get<4>(BS[len-1].dsvp_t) == get<4>(bs.dsvp_t) && BS[len-1].slope < bs.slope)) 
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        //|| (get<2>(BS[mid].dsvp_t) == get<2>(bs.dsvp_t) && BS[mid].slope == bs.slope && BS[mid].cum_pr <= bs.cum_pr)
        if( get<4>(BS[mid].dsvp_t) > get<4>(bs.dsvp_t) || (get<4>(BS[mid].dsvp_t) == get<4>(bs.dsvp_t) && BS[mid].slope < bs.slope) ) 
            left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}




int EnumBS::binary_search_for_dsvp_slope(blocksize_strategy bs){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    # sort in (G2,slope,cum_pr). Find the position to put (G2,slope,cum_pr), G2-->smaller, slope --> larger, cum_pr --> larger.
    */
    
    
    int len = int(BS.size());
    // double G2_tmp = round(pow(2,get<2>(BS[len-1].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec;
 
    //|| (get<2>(BS[len-1].dsvp_t) == get<2>(bs.dsvp_t) && BS[len-1].slope == bs.slope && BS[len-1].cum_pr <= bs.cum_pr)
    if( get<4>(BS[len-1].dsvp_t) > get<4>(bs.dsvp_t) || (get<4>(BS[len-1].dsvp_t) == get<4>(bs.dsvp_t) && BS[len-1].slope < bs.slope)) 
        return len;
    
    int left = 0, right = len - 1;
    int mid = floor((left+right)/2);
    while(left<right){
        // if(round(pow(2,get<2>(BS[mid].dsvp_t))*params->enumbs_G_prec)/params->enumbs_G_prec >= G2) left = mid + 1;
        //|| (get<2>(BS[mid].dsvp_t) == get<2>(bs.dsvp_t) && BS[mid].slope == bs.slope && BS[mid].cum_pr <= bs.cum_pr)
        if( get<0>(BS[mid].dsvp_t) > get<0>(bs.dsvp_t) || (get<0>(BS[mid].dsvp_t) == get<0>(bs.dsvp_t) && BS[mid].slope < bs.slope) ) 
            left = mid + 1;
        else right = mid;
        mid = floor((left+right)/2);
    }
    return left;
}





int EnumBS::binary_search_for_dsvp(int dsvp){
    /*
    # Input: A list whose element is from large to small
    # Return the first index of the first number < dsvp in high dimensional list
    # If it is non-existent, then return len(BS)
    # If all elements >= dsvp, then return len(BS)
    */

    int len = int(BS.size());
    int dsvp_tmp = get<0>(BS[len-1].dsvp_t);
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

// vector<double> EnumBS::extract_cdsvp(){
//     vector<double> cdsvps;
//     cdsvps.resize(0);
//     for(int i =0; i < int(BS.size()); i++){
//         cdsvps.insert(cdsvps.end(), round(get<0>(BS[i].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec);
//     }
//     return cdsvps;
// }

vector<int> EnumBS::extract_dsvp(){
    vector<int> dsvps;
    dsvps.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        dsvps.insert(dsvps.end(), get<0>(BS[i].dsvp_t));
    }
    return dsvps;
}


vector<double> EnumBS::extract_PSC(){
    vector<double> PSCs;
    PSCs.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        PSCs.insert(PSCs.end(), get<4>(BS[i].dsvp_t));
    }
    return PSCs;
}


vector<double> EnumBS::extract_G2(){
    vector<double> G2s;
    G2s.resize(0);
    for(int i =0; i < int(BS.size()); i++){
        // G2s.insert(G2s.end(), round(get<2>(BS[i].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec);
        G2s.insert(G2s.end(),get<4>(BS[i].dsvp_t));
    }
    return G2s;
}


// vector<pair<double,double>> EnumBS::extract_G2G1(){
//     vector<pair<double,double>> G2s;
//     G2s.resize(0);
//     for(int i =0; i < int(BS.size()); i++){
//         // G2s.insert(G2s.end(), round(get<2>(BS[i].dsvp_t)*params->enumbs_G_prec)/params->enumbs_G_prec);
//         G2s.insert(G2s.end(), make_pair(get<2>(BS[i].dsvp_t), BS[i].cum_avg_GB_BKZ.first));
//     }
//     return G2s;
// }


// vector<double> EnumBS::extract_slope(){
//     vector<double> slopes;
//     slopes.resize(0);
//     for(int i =0; i < int(BS.size()); i++){
//         slopes.insert(slopes.end(), BS[i].slope);
//     }
//     return slopes;
// }


// vector<double> EnumBS::extract_cum_pr(){
//     vector<double> cum_prs;
//     cum_prs.resize(0);
//     for(int i =0; i < int(BS.size()); i++){
//         cum_prs.insert(cum_prs.end(), BS[i].cum_pr);
//     }
//     return cum_prs;
// }


// vector<tuple<double,double,double>> EnumBS::extract_G2_slope_cum_pr(){
//     vector<tuple<double,double,double>> basis_quality_list;
//     basis_quality_list.resize(0);
//     for(int i =0; i < int(BS.size()); i++){
//         basis_quality_list.insert(basis_quality_list.end(), make_tuple(get<2>(BS[i].dsvp_t), BS[i].slope, BS[i].cum_pr));
//     }
//     return basis_quality_list;
// }


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



void EnumBS::BS_add(EnumBS::blocksize_strategy bs, int k){

    if(BS.size() == 0){
        BS.insert(BS.end(),bs);
        return;
    }

    int pos = binary_search_for_G2_slope(bs);
    // int pos =  binary_search_for_dsvp_slope(bs);
    // int pos = binary_search_for_G2(get<2>(bs.dsvp_t));

    // cout<<"pos = "<< pos <<endl;
    // if(k == 13)
    //     throw "";
    if(pos == 0){
        return;
    }


    pos--;



    //min avgGB
    // while(pos > 0 && pos < int(BS.size()) && ( get<2>(bs.dsvp_t) < get<2>(BS[pos].dsvp_t) + params->enumbs_G_prec || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && get<2>(bs.dsvp_t) <= get<2>(BS[pos].dsvp_t)  + params->enumbs_G_prec && get<2>(bs.dsvp_t) >= get<2>(BS[pos].dsvp_t) )) && bs.cum_avg_GB_BKZ.first < BS[pos].cum_avg_GB_BKZ.first + params->enumbs_G_prec){
    //Consider min_GB
    // while(pos > 0 && pos < int(BS.size()) && bs.avg_GB.first < BS[pos].avg_GB.first + params->enumbs_G_prec){
    
    //min cumGB
    // while(pos > 0 && pos < int(BS.size()) && ( get<2>(bs.dsvp_t) < get<2>(BS[pos].dsvp_t) + params->enumbs_G_prec || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && get<2>(bs.dsvp_t) <= get<2>(BS[pos].dsvp_t)  + params->enumbs_G_prec && get<2>(bs.dsvp_t) >= get<2>(BS[pos].dsvp_t) )) && bs.cum_GB_BKZ.first < BS[pos].cum_GB_BKZ.first + params->enumbs_G_prec){

    while(pos > 0 && pos < int(BS.size()) && ( get<4>(bs.dsvp_t) < get<4>(BS[pos].dsvp_t) + params->enumbs_G_prec || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && get<4>(bs.dsvp_t) <= get<4>(BS[pos].dsvp_t)  + params->enumbs_G_prec && get<4>(bs.dsvp_t) >= get<4>(BS[pos].dsvp_t) )) && bs.cum_GB_BKZ.first < BS[pos].cum_GB_BKZ.first + params->enumbs_G_prec){ //&& bs.S[bs.S.size()-1].beta <= BS[pos].S[BS[pos].S.size()-1].beta){

    // while(pos > 0 && pos < int(BS.size()) &&  (get<0>(bs.dsvp_t) < get<0>(BS[pos].dsvp_t) || (bs.slope > BS[pos].slope - params->enumbs_slope_prec && get<0>(bs.dsvp_t) == get<0>(BS[pos].dsvp_t)) ) && bs.cum_GB_BKZ.first < BS[pos].cum_GB_BKZ.first){

        if(params->debug){
            // if(BS[pos].S[0].beta == 50 && BS[pos].S[0].jump == 1 && BS[pos].S[0].tours == 1 && BS[pos].S.size() ==1 ){
            {
                
                // cout<<"k = "<<k<<endl;
                if(BS[pos].S.size()>=1 && BS[pos].S[0].beta == 50 and BS[pos].S[0].jump == 1){
                    printf("=======erase1=====\n");
                    printf("erase strategy:\n");
                    print_bs(BS[pos]);
                    printf("May added strategy:\n");
                    print_bs(bs);
                    printf("==================\n");
                    // throw "";
                }
            }
        }
        BS.erase(BS.begin()+pos);
        pos--;
        // if(BS[pos].cum_avg_GB_BKZ.first > 0.1){
        //     BS.erase(BS.begin()+pos);
        //     pos--;
        // }
        // else
        //     break;
    }


    if(params->debug){
        // if(bs.S[0].beta == 50 && bs.S[0].jump == 1 && bs.S[0].tours == 1 && bs.S.size() ==1 ){
        {
            if(bs.S.size()>=1 && bs.S[0].beta == 50 && bs.S[0].jump == 1){
                printf("=======add=====\n");
                cout<<"k = "<<k<<endl;
                // if( bs.S.size() <=1){
                print_bs(bs);
                // }    
                printf("==================\n");
                // sleep(1);
            }
        }
    }




    BS.insert(BS.begin()+pos+1,bs);
    


    //Consider min_GB
    // while( pos+1<int(BS.size())-1 && BS[pos+2].avg_GB.first < BS[pos+1].avg_GB.first + params->enumbs_G_prec){
    //min avgGBKZ
    // while( pos+1<int(BS.size())-1 && pos+1 > k && ( get<2>(BS[pos+2].dsvp_t) < get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec  &&  get<2>(BS[pos+2].dsvp_t) == get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec)) && BS[pos+2].cum_avg_GB_BKZ.first < BS[pos+1].cum_avg_GB_BKZ.first + params->enumbs_G_prec){
    
    
    //min cumGBKZ
    // while( pos+1<int(BS.size())-1 && ( get<2>(BS[pos+2].dsvp_t) < get<2>(BS[pos+1].dsvp_t) + params->enumbs_G_prec || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec  &&  get<2>(BS[pos+2].dsvp_t) <= get<2>(BS[pos+1].dsvp_t)  + params->enumbs_G_prec && get<2>(BS[pos+2].dsvp_t) >= get<2>(BS[pos+1].dsvp_t) )) && BS[pos+2].cum_GB_BKZ.first < BS[pos+1].cum_GB_BKZ.first + params->enumbs_G_prec){


    while( (pos+1<int(BS.size())-1 && ( get<4>(BS[pos+2].dsvp_t) < get<4>(BS[pos+1].dsvp_t) + params->enumbs_G_prec || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec  &&  get<4>(BS[pos+2].dsvp_t) <= get<4>(BS[pos+1].dsvp_t)  + params->enumbs_G_prec && get<4>(BS[pos+2].dsvp_t) >= get<4>(BS[pos+1].dsvp_t) )) && BS[pos+2].cum_GB_BKZ.first < BS[pos+1].cum_GB_BKZ.first + params->enumbs_G_prec)){  //&& BS[pos+2].S[BS[pos+2].S.size()-1].beta <= BS[pos+1].S[BS[pos+1].S.size()-1].beta ){
    
    // while( pos+1<int(BS.size())-1 && ( get<0>(BS[pos+2].dsvp_t) < get<0>(BS[pos+1].dsvp_t) || (BS[pos+2].slope > BS[pos+1].slope - params->enumbs_slope_prec  and get<0>(BS[pos+2].dsvp_t) == get<0>(BS[pos+1].dsvp_t ) )) && BS[pos+2].cum_GB_BKZ.first < BS[pos+1].cum_GB_BKZ.first ){
    
      
        if(params->debug){
            // if(BS[pos+1].S[0].beta == 50 && BS[pos+1].S[0].jump == 1 && BS[pos+1].S[0].tours == 1 ){
            {
                
                // cout<<"k = "<<k<<endl;
                if(BS[pos+1].S.size()>=1 && BS[pos+1].S[0].beta == 50 && BS[pos+1].S[0].jump == 1){
                    printf("=======erase2=====\n");
                    cout<<"delete one: ";
                    print_bs(BS[pos+1]);
                    cout<<"The next one: ";
                    print_bs(BS[pos+2]);
                    printf("==================\n");
                    // throw "";
                }
            // sleep(1);
            }
        }

        BS.erase(BS.begin()+pos+1);
        pos--;
        // if(BS[pos+1].cum_avg_GB_BKZ.first > 0.1){
        //     BS.erase(BS.begin()+pos+1);
        //     pos--;
        // }
        // else
        //     break;
        
        
    }
    

    // if(pos < k)
    //     k = pos - 1;
    if(pos <= k)
        k = -1;

    if(params->debug){
        // vector<int> dsvps = extract_dsvp();
        // vector<int> sorted_dsvps = dsvps;
        // sort(sorted_dsvps.rbegin(),sorted_dsvps.rend());
        // if(dsvps != sorted_dsvps){
        //     printf("dsvps: \n");
        //     print_vector(dsvps);
        //     printf("sorted_dsvps: \n");
        //     print_vector(sorted_dsvps);
        // }
        // assert(dsvps == sorted_dsvps);



        vector<double> PSCs = extract_PSC();
        vector<double> sorted_PSCs = PSCs;
        sort(sorted_PSCs.rbegin(),sorted_PSCs.rend());
        if(PSCs != sorted_PSCs){
            printf("PSCs: \n");
            print_vector(PSCs);
            printf("sorted_PSCs: \n");
            print_vector(sorted_PSCs);
        }
        assert(PSCs == PSCs);
        // print_vector(G2s);
        // assert(no_repeated_value_verification(G2s));
    }
    // if(params->debug){
        
    //     vector<double> G2s = extract_G2();
    //     vector<double> sorted_G2s = G2s;
    //     sort(sorted_G2s.rbegin(),sorted_G2s.rend());
    //     assert(G2s == sorted_G2s);
    //     // print_vector(G2s);
    //     // assert(no_repeated_value_verification(G2s));
    // }
}




bool EnumBS::pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB_BKZ, pair<double,double> &cum_avg_GB_BKZ, pair<double,double> &avg_GB, double &cum_pr, int beta, int jump, tuple<int,int,double,double,double> &dsvp_t_, double &slope, pair<double,double> &GB_nopr, pair<double,double> &min_GB, bool &leaf){

    double rem_pr = 1. - cum_pr;

    vector<double> l_;

    int d = l.size(), beta_ = get_beta_(params,beta,jump,d);

    

    sim -> simulate(l_,l,beta,jump,1);
   

    //Delete this condition to ensure (50,1,1) or other strategy with small blocksize can be added into BS.
    if(l_[d-beta_] == l[d-beta_]){
        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model, cum_GB_BKZ);
        slope = get_current_slope(l, 0, d);
        return false;
    }
    else{
    // if(true){
        l = l_;
        slope = get_current_slope(l, 0, d);
        boost::math::chi_squared chisquare(beta_);
        double pr = boost::math::cdf(chisquare,pow(2,2.*l[d-beta_]));
        
        
        pair<double,double> GB_BKZ = cost->bkz_cost(d,beta,jump,params->cost_model);

        cum_GB_BKZ.first = log2(pow(2,cum_GB_BKZ.first)+pow(2,GB_BKZ.first));
        cum_GB_BKZ.second = max(cum_GB_BKZ.second, GB_BKZ.second);
        if(not params->worst_case){
            cum_avg_GB_BKZ.first = log2(pow(2,cum_avg_GB_BKZ.first)+(pow(2,cum_GB_BKZ.first)*rem_pr*pr));
            cum_avg_GB_BKZ.second = log2(pow(2,cum_avg_GB_BKZ.second)+(pow(2,cum_GB_BKZ.second)*rem_pr*pr));
            // cum_avg_GB_BKZ.second = max(cum_avg_GB_BKZ.second, cum_GB_BKZ.second);
        }
        else{
            // cum_avg_GB_BKZ.first = log2(pow(2,cum_avg_GB_BKZ.first)+(pow(2,cum_GB_BKZ.first)*(pr-pre_pr)));
            // cum_avg_GB_BKZ.second = log2(pow(2,cum_avg_GB_BKZ.second)+pow(2,GB_BKZ.second) * (pr-pre_pr));
            cum_avg_GB_BKZ = cum_GB_BKZ;
        }
        

        // printf("cum_G = %e\n", cum_avg_GB_BKZ.first );
        cum_pr += rem_pr * pr;
        rem_pr = 1. - cum_pr;
        pre_pr = pr;

        dsvp_t_ = dsvp_predict(l, cum_pr, cost,params->cost_model,  cum_GB_BKZ);
        
        avg_GB.first = log2(pow(2,cum_avg_GB_BKZ.first) + pow(2,get<2>(dsvp_t_)));
        avg_GB.second = log2(pow(2,cum_avg_GB_BKZ.second)+pow(2, get<3>(dsvp_t_)));

        GB_nopr.first = log2(pow(2,cum_GB_BKZ.first) + pow(2,get<4>(dsvp_t_)));
        GB_nopr.second = max(cum_GB_BKZ.second, get<3>(dsvp_t_));



        if(params->enumbs_min_G){
            if(min_GB.first > GB_nopr.first){
                min_GB = GB_nopr;
                leaf = false;
            }
            else if(GB_nopr.first >= min_GB.first + params->min_G_prec)
                leaf = true;
        }else{
            // or min_GB.second > params->max_RAM
            if(min_GB.second > GB_nopr.second){
                min_GB = GB_nopr;
                leaf = false;
            }
            else if(GB_nopr.second >= params->max_RAM and GB_nopr.first >= min_GB.first + params->min_G_prec)
                leaf = true;
        }

       return true;
    }

}


void EnumBS::max_tour_for_pnjbkz_beta(int k, int beta,int jump){
    EnumBS::blocksize_strategy bs = BS[k];
    tuple<int,int,double,double,double> dsvp_t1;
    // double G21;
    // double G20 = get<2>(bs.dsvp_t), G21;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;

    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ=bs.cum_avg_GB_BKZ, avg_GB = bs.avg_GB, min_GB = bs.min_GB, GB_nopr = bs.GB_nopr;
    bool leaf = bs.leaf;

    int  loop = 0;
    
    double slope0 = bs.slope, slope1;
    bool sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, beta, jump, dsvp_t1,slope1, GB_nopr, min_GB, leaf);

    if(not sim_term)
        return;
    
    // G21 = get<2>(dsvp_t1);

    // assert(G21 >= 0.);

    
    while(not leaf && ((slope1 - slope0) > params->enumbs_slope_prec) && loop < params->max_loop ){

        loop +=1;
        slope0 = slope1;
        
        if(loop ==1){
            S.insert(S.end(),{beta, jump, loop});
        }
        else{
            S[S.size()-1].tours = loop;
        }
        
        // if(params->enumbs_min_G){
        //     // if(log2(pow(2,cum_GB_BKZ)+pow(2,)
        //     if(min_GB.first > GB_nopr.first){
        //         min_GB = GB_nopr;
        //         leaf = false;
        //     }
        //     else if(GB_nopr.first >= min_GB.first + params->min_G_prec)
        //         leaf = true;
        // }else{
        //     // or min_GB.second > params->max_RAM
        //     if(min_GB.second > GB_nopr.second){
        //         min_GB = GB_nopr;
        //         leaf = false;
        //     }
        //     else if(GB_nopr.second >= params->max_RAM and GB_nopr.first >= min_GB.first + params->min_G_prec)
        //         leaf = true;
        // }

        if(params->enumbs_min_G){
            if(cum_pr >= params->succ_prob or not sim_term or leaf)
                break;
        }

        bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr,slope1,min_GB,GB_nopr,leaf};

        // if(params->verification){
            // pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
            //cerr<<"cum_pr="<<cum_pr<<", verified cum_pr="<< verified_cum_G_pr.second<<endl;
            //cerr<<"cum_G="<<cum_avg_GB_BKZ.first<<", verified cum_G="<< verified_cum_G_pr.first<<endl;
        //     assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
        //     assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
        // }
    
        // EnumBS::BS_add_G2(bs, k);
        // EnumBS::BS_add_slope(bs,k);
        // BS_add_op(bs,k);
        // if(params->enum_add_G2)
        //     BS_add_G2(bs, k);
        // else
        BS_add(bs,k);

        if(not params->enumbs_min_G){
            if(cum_pr >= params->succ_prob or not sim_term or leaf)
                break;
        }
        
    
        sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, beta, jump, dsvp_t1,slope1,GB_nopr, min_GB, leaf);
        // G21 = get<2>(dsvp_t1);
        // assert(G21 >= 0.);
        // G21 = get<2>(dsvp_t1);
    }
}


void EnumBS::max_tour_for_pnjbkz_beta_in_parallel( int beta_j_t_id_begin, vector<pair<int,int>> beta_j_tid,  int k){

    // bool flag = true;

    for(int i = 0; i< int(beta_j_tid.size()); i++){
        EnumBS::blocksize_strategy bs = BS[k];
        tuple<int,int,double,double,double> dsvp_t1;
        // double G20 = get<2>(bs.dsvp_t), G21;
        vector<double> l = bs.l; 
        vector<strategy> S = bs.S;

        double cum_pr = bs.cum_pr;
        pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ = bs.cum_avg_GB_BKZ, avg_GB = bs.avg_GB, min_GB = bs.min_GB, GB_nopr = bs.GB_nopr;
        bool leaf = bs.leaf;
        int beta = beta_j_tid[i].first, jump = beta_j_tid[i].second;


        if(bs.leaf)
            continue;

        // cout<<beta<<","<<jump<<endl;
        assert(beta>0 and beta <= min(params->max_dim,int(l.size())));
        assert(jump>0 and jump<=params->J);

        
        
        int jub = jump_upper_bound(params,beta,l);
    
        //((f == 0 or beta < 79 )&& jump > 1) or
        if(jub!=0 && jump > jub)
        // if( ((f == 0 or beta < 79 )&& jump > 1) or (f!=0 && jump >= min((double) f,ceil(0.1*beta))))//ceil(0.1*beta)
            continue;

        if(beta < 55 and jump > 1)
            continue;

        // if(params->cost_model == 1 and jump >= 0.1 * beta)
        //     continue;
    

        int index = beta_j_t_id_begin + i;
        
        int  loop = 0;

        double slope0 = bs.slope, slope1;


        bool sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, beta, jump, dsvp_t1, slope1, GB_nopr, min_GB, leaf);

      


        if(not sim_term)
            continue;

        tmpBS[index].clear();

        // while((abs(G20 - params->max_num)<0.01 or G20 > G21 or (slope1 - slope0) > 1e-6) && loop < params->max_loop){

        // bool leaf = bs.leaf; //If leaf = false, we should further add the strategy; otherwise, we can skip this adding.

        
        
        
        // if(beta == 117 and jump == 8 and bs.S.size() == 0)
        //     cout<< 117<<","<<jump<<","<<pow(2,cum_avg_GB_BKZ.first)/3600 <<endl;
    

        while((slope1 - slope0) > params->enumbs_slope_prec && loop < params->max_loop){

        // while(loop < params->max_loop){
            loop +=1;
            // G20 = G21;
            slope0 = slope1;
            if(beta == 1083129856)
                throw "Error!";
            

            assert(loop>=1);

            if(loop == 1){
                S.insert(S.end(),{beta, jump, loop});
            }
            else
                S[S.size()-1].tours = loop;

    
            bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, slope1, min_GB, GB_nopr, leaf};


            // if(leaf)
            //     print_bs(bs);

        

            // if(beta  == 95 and jump == 2 and bs.S.size() == 1){
            //         cout<<"======================"<<endl;
            //         cout<<"beta = "<<beta<<", jump = "<<jump;
            //         cout<<", loop = "<<loop<<endl;
            //         // cout<<"dsvp_prime = "<< get<0>(dsvp_t1)<<", cost:"<<get<2>(dsvp_t1)<<endl;
            //         print_bs(bs);
                
            //         cout<<"params->enumbs_min_G: "<<params->enumbs_min_G <<", leaf: " <<leaf<<endl;
            //         cout<<"sim term = "<<sim_term<<endl;

            //         cout<<"cum_pr = "<<cum_pr<<endl;
                    
            //         cout<<"slope1 = "<<slope1<<", slope0 = "<<slope0<<endl;
            //         cout<<"....................."<<endl;
            //         // throw "";
            // }

            
            // if(params->verification){
            //     // pair<double,double> verified_cum_G_pr = strategy_verification(l0,S);
            //     // assert(abs(verified_cum_G_pr.first-cum_avg_GB_BKZ.first)<0.001);
            //     // assert(abs(verified_cum_G_pr.second-cum_pr)<0.001);
            // }            

            int len_tmpBS = tmpBS[index].size();

            //dsvp min 
            // if( loop > 1 && (get<0>(bs.dsvp_t) < get<0>(tmpBS[index][len_tmpBS - 1].dsvp_t)  || (get<0>(bs.dsvp_t) == get<0>(tmpBS[index][len_tmpBS - 1].dsvp_t)  &&  bs.slope > tmpBS[index][len_tmpBS - 1].slope - params->enumbs_slope_prec )) && bs.cum_GB_BKZ.first < tmpBS[index][len_tmpBS - 1].cum_GB_BKZ.first){


            if( loop > 1 && (get<4>(bs.dsvp_t) < get<4>(tmpBS[index][len_tmpBS - 1].dsvp_t) + params->enumbs_G_prec  || (get<4>(bs.dsvp_t) == get<4>(tmpBS[index][len_tmpBS - 1].dsvp_t) + params->enumbs_G_prec && bs.slope > tmpBS[index][len_tmpBS - 1].slope - params->enumbs_slope_prec )) && bs.cum_GB_BKZ.first < tmpBS[index][len_tmpBS - 1].cum_GB_BKZ.first + params->enumbs_G_prec){
            // // if( loop > 1 && ((get<0>(bs.dsvp_t) < get<0>(tmpBS[index][len_tmpBS - 1].dsvp_t) || (get<0>(bs.dsvp_t) == get<0>(tmpBS[index][len_tmpBS - 1].dsvp_t) && bs.slope > tmpBS[index][len_tmpBS - 1].slope - params->enumbs_slope_prec ))) && bs.cum_GB_BKZ.first < tmpBS[index][len_tmpBS - 1].cum_GB_BKZ.first + params->enumbs_G_prec){
                tmpBS[index][len_tmpBS-1] = bs;   
            }
            else
                tmpBS[index].insert(tmpBS[index].end(),bs);

           
            // if(params->enumbs_min_G and leaf){
            //     break;
            // }

            if(params->enumbs_min_G and leaf){
                break;
            }

            if(cum_pr >= params->succ_prob or not sim_term)
                break; //Cannot return, will omit the later blocksize strategy
            

            
            
            
            sim_term = pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, beta, jump, dsvp_t1, slope1, GB_nopr, min_GB, leaf);
        
            // G21 = get<2>(dsvp_t1);
            // assert(G21 >= 0.);
        }
    }
}


void EnumBS::enumbs_est(vector<double> l0){
    /*
    input: l -- gs-lengths; 
    Return: Optimal strategy to minimize solving cost.
    */
    int beta_start = params->beta_start, k = 0, d = l0.size(),j_start,len_S;
    
    tuple<int,int,double,double,double>  dsvp0_t = dsvp_predict(l0, -1000., cost, params->cost_model);
    blocksize_strategy bs =  {dsvp0_t, {},l0,make_pair(-1000.,-1000.), make_pair(-1000.,-1000.), make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t)), 0., get_current_slope(l0,0,d), make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t))};


    BS.insert(BS.end(),bs);

    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

        //Add a normal two-step strategy from (beta_start,1,1) to (d,1,1)


    
    tuple<int,int,double,double,double> dsvp_t1;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ = bs.cum_avg_GB_BKZ, avg_GB = bs.avg_GB, min_GB = bs.min_GB, GB_nopr = bs.GB_nopr;
    double slope1;
    bool leaf = false;
    
    for(int beta = beta_start;  beta < d; beta++){
        pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, beta, 1, dsvp_t1, slope1, GB_nopr, min_GB, leaf);
        S.insert(S.end(),{beta,1,1}); 
        // if(params->enumbs_min_G){
        //     if(min_GB.first > GB_nopr.first){
        //         min_GB = GB_nopr;
        //         leaf = false;
        //     }
        //     else if(GB_nopr.first >= min_GB.first + params->min_G_prec)
        //         leaf = true;
        // }else{
        //     // or min_GB.second > params->max_RAM
        //     if(min_GB.second > GB_nopr.second){
        //         min_GB = GB_nopr;
        //         leaf = false;
        //     }
        //     else if(GB_nopr.second >= params->max_RAM and GB_nopr.first >= min_GB.first + params->min_G_prec)
        //         leaf = true;
        // }

        // if(params->enumbs_min_G){
        //     if(cum_pr >= params->succ_prob or leaf)
        //         break;
        // }
    
        bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, slope1, min_GB, GB_nopr, leaf};
        BS_add(bs,k);

        if(not params->enumbs_min_G){
            if(cum_pr >= params->succ_prob or leaf)
                break;
        }
        k++;
    }

    k = 0;
    while( k < int(BS.size())){
        
        bs = BS[k];
        

        //skip leaf strategy.
        if(bs.leaf)
            continue;
        
        len_S = bs.S.size();

        if(bs.cum_pr >= params->succ_prob){
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
                beta_start = max( params->beta_start, beta_start - params->delta_beta);
                // beta_start += 1;
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

                int jub = jump_upper_bound(params,beta,bs.l);
                
                //((f == 0 or beta < 79 )&& j > 1) or
                if(jub!=0 && j > jub)
                // if(( (f == 0 or beta < 79 )&& j > 1) or (f!=0 && j >= min((double)f,ceil(0.1*beta))))
                    continue;

                if(beta < 55 and j > 1)
                    continue;
                // if(params->cost_model == 1 and j >= 0.1 * beta)
                //     continue;
        

                max_tour_for_pnjbkz_beta(k,beta,j); 
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
        
        if(BS[i].avg_GB.first<Gmin and BS[i].avg_GB.second <= params->max_RAM){
            bsmin = BS[i];
            Gmin = BS[i].avg_GB.first;
            Bmin = BS[i].avg_GB.second;
        }
    }
    printf("Find the optimal Strategy through EumBS!!\n");
    print_bs(bsmin);
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log2(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec) = %3.2f s, Memory Cost = %3.2f log2(bit) = %3.2f avg_GB \n", Gmin, pow(2,Gmin), Bmin, pow(2,Bmin-33));
        // printf("Min Cost = %3.2f log2(sec) = %3.2f h, Memory Cost = %3.2f log2(bit) = %3.2f TB \n", Gmin, pow(2,Gmin)/3600, Bmin, pow(2,Bmin-43));
}


void EnumBS::enumbs_est_in_parallel(vector<double> l0){
    /*
    input: l0 -- gs-lengths;
    Return: Optimal strategy to minimize solving cost.
    */
    set_threads(params->threads);

    int k = 0, d = l0.size(), beta_start = params->beta_start;
    blocksize_strategy bs;

    tuple<int,int,double,double,double>  dsvp0_t = dsvp_predict(l0, -1000., cost,params->cost_model);

    bs = {dsvp0_t, {},l0,make_pair(-1000.,-1000.), make_pair(-1000.,-1000.),make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t)), 0., get_current_slope(l0,0,d), make_pair(get<2>(dsvp0_t), get<3>(dsvp0_t))};

    BS.insert(BS.end(),bs);


    //Add a normal two-step strategy from (beta_start,1,1) to (d,1,1)


    /*
    tuple<int,int,double,double,double> dsvp_t1;
    vector<double> l = bs.l; 
    vector<strategy> S = bs.S;
    double cum_pr = bs.cum_pr;
    pair<double,double> cum_GB_BKZ = bs.cum_GB_BKZ, cum_avg_GB_BKZ = bs.cum_avg_GB_BKZ, avg_GB = bs.avg_GB, min_GB = bs.min_GB;
    double slope1;
    bool leaf = false;
    for(int beta = beta_start;  beta < d; beta++){
        pnjbkz_beta_loop(l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, beta, 1, dsvp_t1, slope1);
        S.insert(S.end(),{beta,1,1}); 
        // cout<<"cum_pr = "<<cum_pr<<endl;
        // cout<<"beta = "<<beta<<", enumbs_min_G: "<<params->enumbs_min_G<<endl;
        // cout<<get<0>(dsvp_t1)<<endl;
        // cout<<"avg_GB.first = "<<avg_GB.first<<endl;
        // cout<<"min_GB.first = "<<min_GB.first<<endl;
        if(params->enumbs_min_G){
            if(min_GB.first > avg_GB.first){
                min_GB = avg_GB;
                leaf = false;
            }
            else if(avg_GB.first >= min_GB.first + params->min_G_prec)
                leaf = true;
        }else{
            // or min_GB.second > params->max_RAM
            if(min_GB.second > avg_GB.second){
                min_GB = avg_GB;
                leaf = false;
            }
            else if(avg_GB.second >= params->max_RAM and avg_GB.first >= min_GB.first + params->min_G_prec)
                leaf = true;
        }
        if(params->enumbs_min_G and (cum_pr >= params->succ_prob or leaf))
            break;

    
        bs = {dsvp_t1, S, l, cum_GB_BKZ, cum_avg_GB_BKZ, avg_GB, cum_pr, slope1, min_GB};

        BS_add(bs,k);
        
        if(params->enumbs_min_G and (cum_pr >= params->succ_prob or leaf))
            break;
        k++;
    }
    */


    int j,len_S;
    params->J = floor((params->J-1)/params->J_gap) * params->J_gap+1;

    k = 0;
    while( k < int(BS.size())){
        // sleep(10);
        bs = BS[k];
        //skip leaf strategy.
       


        len_S = bs.S.size();




        // print_bs(bs);
        // cout<<"Current minG = "<<bs.min_GB.first<<endl;
        // if(params->debug){
        //     printf("==========BS============\n");
        //     print_BS(BS);
        //     printf("======================\n");
        //     sleep(10);
        // }

        // printf("debug = %d\n ", params->debug);
        // if(params->debug){
        //     if(bs.S.size()>400){
        //         printf("======================\n");
        //         print_bs(bs);
        //         // sleep(1);
        //         printf("======================\n");
        //     }
        // } 

        if(bs.leaf){
            printf("======================\n");
            print_bs(bs);
            printf("----------------------\n");
        }

        if(bs.cum_pr >= params->succ_prob or bs.leaf){
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
                // beta_start = max(params->beta_start, beta_start - params->delta_beta);
                j = params->J;
            }
            // else
            //     j -= params->J_gap;
            else{
                beta_start = max(params->beta_start, beta_start-jump_upper_bound(params,beta_start,bs.l));
                j = params->J;    
            } // It will be very slow if add this condition. 
        }

        //int(0.9*d)
        int t_id = 0, len = (((j-1)/params->J_gap)+1)+(min(params->max_dim,d)-1-beta_start)/params->gap*(((params->J-1)/params->J_gap)+1);

        // cout<<len<<endl;
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
            for(int beta = beta_start; beta < min(params->max_dim, d); beta +=params->gap){
                for(; j>0; j-= params->J_gap){
                    if(int(beta_j[t_id].size()) ==  departs[t_id] && t_id < threads-1){
                        t_id++;
                    }
                    beta_j[t_id].insert(beta_j[t_id].end(),make_pair(beta,j));
                }
                j = params->J;
            }
        
            // if(params->debug){
            //     int Sum = 0;
            //     for(int t_id = 0; t_id  < int(beta_j.size()); t_id ++){
            //         Sum += beta_j[t_id ].size();
            
            //         // cout<<beta_j[t_id][0].first<<", "<<beta_j[t_id][0].second<<endl;
            //         // cerr<<"t_id = "<<t_id<<", beta_j[i].size() = "<<beta_j[t_id].size()<<endl;
            //     }
            //     // cerr<<len<<","<<Sum<<endl;
            //     assert(len == Sum);
            // }

           
            auto start= system_clock::now();
            for (t_id = 0; t_id < threads; t_id++){
                threadpool.push([this, t_id, beta_j_t_id_begins, beta_j, k](){
                    max_tour_for_pnjbkz_beta_in_parallel(beta_j_t_id_begins[t_id], beta_j[t_id], k); 
                });            
            }
            threadpool.wait_work(); 
            if(params->verbose){
                auto finish = system_clock::now();
                duration<double> diff = finish - start;
                if(k > 0){
                    printf("\r index: %8d, (%4d,%4d): (%4d,%4d) --> (%4d,%4d), goal index: %8d, cost = " ,k+1,bs.S[len_S-1].beta, bs.S[len_S-1].jump, beta_j[0][0].first,beta_j[0][0].second,beta_j[threads-1][beta_j[threads-1].size()-1].first,beta_j[threads-1][beta_j[threads-1].size()-1].second,int(BS.size()));
                    cerr<<setprecision(2)<<diff.count()<<'s';
                }
            }

            
           
            for(int i = 0; i< len; i++){
                for(int ii = 0; ii < int(tmpBS[i].size()); ii++){
                    // if(params->enum_add_G2)
                    //     BS_add_G2(tmpBS[i][ii], k);
                    // else
                    // if(tmpBS[i][ii].S[0].beta==50 && tmpBS[i][ii].S[0].jump==1)
                    //     printf("beta=%d, jump=%d", tmpBS[i][ii].S[0].beta, tmpBS[i][ii].S[0].jump);
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

    //Find strategy with minimal cost
    double Gmin = params->max_num, Bmin  = params->max_num;
    EnumBS::blocksize_strategy bsmin;
    bool flag = false;
    // print_BS(BS);
    for(int i = 0; i<int(BS.size()); i++){
        if(BS[i].avg_GB.first<Gmin  and BS[i].avg_GB.second < params->max_RAM){
            bsmin = BS[i];
            Gmin = BS[i].avg_GB.first;
            Bmin = BS[i].avg_GB.second;
            flag = true;
        }
    }
    if(flag){
        printf("Find the optimal Strategy through EumBS!!\n");
    }
    else{
        for(int i = 0; i<int(BS.size()); i++){
            if(BS[i].avg_GB.second<Bmin){
                bsmin = BS[i];
                Gmin = BS[i].avg_GB.first;
                Bmin = BS[i].avg_GB.second;
            }
        }
        printf("There's no strategy whose memory cost (%lf log(bit)) is below %lf log(bit), the lowest memory cost strategy is:\n", Bmin, params->max_RAM);
    }
    print_bs(bsmin);
    if(params->cost_model == 1)
        printf("Min Cost = %3.2f log2(gate), Memory Cost = %3.2f log(bit)\n", Gmin, Bmin);
    if(params->cost_model == 2)
        printf("Min Cost = %3.2f log2(sec) = %3.2f s, Memory Cost = %3.2f log2(bit) = %3.2f avg_GB \n", Gmin, pow(2,Gmin), Bmin, pow(2,Bmin-33));
        // printf("Min Cost = %3.2f log2(sec) = %3.2f h, Memory Cost = %3.2f log2(bit) = %3.2f TB \n", Gmin, pow(2,Gmin)/3600, Bmin, pow(2,Bmin-43));
}



// pair<double,double> EnumBS::strategy_verification(vector<double> l,vector<strategy> S){

//     int d = l.size();
//     double cum_pr = 0., rem_pr = 1., proba, G1cum=0., B1cum = 0.;
//     // BKZJSim* sim = new BKZJSim(params);
//     // COST* cost = new COST();
//     for(int i = 0; i< int(S.size()); i++){
//         EnumBS::strategy bs = S[i];
//         int beta = bs.beta, jump = bs.jump, N = bs.tours;
//         for(int tour = 0; tour < N; tour++){
        
//             sim -> simulate(l,l,beta,jump,1);

//             boost::math::chi_squared chisquare(beta);
//             proba = boost::math::cdf(chisquare,pow(2,2.*l[d-beta]));
            

//             pair<double,double> avg_GB = cost -> bkz_cost(d,beta,jump,params->cost_model);
//             G1cum = log2(pow(2,G1cum) + (pow(2,avg_GB.first) * rem_pr * proba));
//             B1cum = max(B1cum,avg_GB.second);

//             cum_pr += rem_pr * proba;
//             rem_pr *= 1. - proba;
//         }
//     }

//     tuple<int,int,double,double,double> dsvp_t =  dsvp_predict(l, cum_pr, cost,params->cost_model,  params->worst_case);
//     double G2 = get<2>(dsvp_t);
//     double G = log2(pow(2,G1cum)+pow(2,G2));
//     printf("Verified cum_pr = %e \n ", cum_pr);
//     printf("Verified G1 = %e, G2 = %e, dsvp = %e\n", G1cum,G2,get<0>(dsvp_t));
//     printf("G = %e\n", G );

//     return make_pair(G1cum, cum_pr);

// }
