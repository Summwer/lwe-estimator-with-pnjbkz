#include "d_svp_prediction.h"


tuple<int,int,double,double,double> dsvp_predict(vector<double> l,  double cum_pr, COST* cost, int cost_model, pair<double,double> cum_GB_BKZ ){
    /*
    return dsvp, G, B. 
    */
    switch(cost->params->dsvp_predict_param){
        case 1:
            return progressive_dsvp_predict(l, cum_pr, cost, cost_model, cum_GB_BKZ);
        case 2: 
            return average_dsvp_predict(l, cum_pr, cost, cost_model, cum_GB_BKZ);
        case 3: 
            return average_asvp_dsvp_predict(l, cum_pr, cost, cost_model, cum_GB_BKZ);
        default:
            throw "error: dsvp_predict should be a number among 1.chi-square est,2. svp est,3. asvp est\n";
            return make_tuple(-1,-1,-1,-1,-1); 
    }

}


// tuple<int,int,double,double,double> fixed_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model, bool worst_case){
//     /*
//     return dsvp, G, B. 
//     */
//     int d = l.size();
//     double psvp, p, rp, gh, dsvp_;
//     pair<double,double> p_cost;
//     if(cum_pr >= 0.999){
//         return make_tuple(0.,0,0.,0.,);
//     }
//     for(int dsvp = 50; dsvp <= d; dsvp++ ){
//         //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
//         gh = gaussian_heuristic_log2(l,d-dsvp);
        
        
//         boost::math::chi_squared chisquare(dsvp);
//         psvp = boost::math::cdf(chisquare,gh); //Compute chi-squared value
        
//         p = cum_pr + (1. - cum_pr)* psvp;
//         rp = 1. - p;
//         if(rp < 0.001){
//             dsvp_ = dsvp + dsvp * rp;
//             p_cost = cost->pump_cost(dsvp,cost_model);
//             if(not worst_case)
//                 return  make_tuple(dsvp_, dsvp, p_cost.first * ((1-cum_pr) * psvp + rp), p_cost.second);  //Avoid too small of dsvp
//             else
//                 return  make_tuple(dsvp, dsvp, p_cost.first, p_cost.second);
//         }
//     }
//     p_cost = cost->pump_cost(d,cost_model);
//     return make_tuple(d, d, p_cost.first, p_cost.second);   
// }



//no cum G with model 1
tuple<int,int,double,double,double> progressive_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model, pair<double,double> cum_GB_BKZ ){
    /*
    return dsvp, G, B. 
    */
    int d = l.size(), f = 0 , dsvp_prime = l.size(); //, dsvp_prime_;
    double psvp, pre_psvp = 0.,  gh; //, avg_d_svp  = 0.; pre_psvp2 = 0., psvp2, 
    pair<double,double> p_cost ={-1000.,-1000.};
    double G_cum = -1000., B_cum = -1000., PSC = -1000.;
    // bool flag = false;
    

    if(cum_pr >= cost->params->succ_prob){
        return make_tuple(0.,0,-1000.,-1000.,-1000.);
    }

    bool flag = false;
    vector<double> Gcums = {} , pcums = {};
    // cout<<"cum_GB_BKZ.first = "<<cum_GB_BKZ.first <<endl;

    for(int dsvp = 30; dsvp <= d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic_log2(l,d-dsvp);
        
        boost::math::chi_squared chisquare(dsvp);
        psvp = boost::math::cdf(chisquare,gh); //Compute chi-squared value
        // if(dsvp = d){
        //     cout<<"dsvp = "<<dsvp <<", psvp = "<<psvp<<endl;
        //     cout<<"gh = "<<gh<<endl;
        //     // cout<<"sigma = "<<sigma<<endl;
        // }
        if(pre_psvp >= psvp)
            continue;
        

        
        if(not flag){
            f = get_f_for_pump(cost->params,dsvp);
            dsvp_prime = floor(dsvp - f);
            p_cost = cost->pump_cost(dsvp_prime,cost_model);
            
            G_cum = log2(pow(2,G_cum)+ (pow(2,p_cost.first) + pow(2,cum_GB_BKZ.first)) * (1- cum_pr) * (psvp-pre_psvp));
            
            B_cum = log2(pow(2,B_cum)+pow(2,p_cost.second) * (1 - cum_pr) * (psvp-pre_psvp));
            PSC = log2(pow(2,PSC)+pow(2,p_cost.first) * (psvp-pre_psvp));

        }
        if(cost->params->print_Gcums){
            if(cum_pr + (1-cum_pr) * psvp> 1e-4){
                Gcums.insert(Gcums.end(), G_cum);
                pcums.insert(pcums.end(), cum_pr + (1-cum_pr) * psvp);
            }
        }
        if(cum_pr + (1-cum_pr) * psvp >= cost->params->succ_prob)
            flag = true;
        // if(cum_pr + (1-cum_pr) * psvp >= cost->params->succ_prob)

        // printf("%d, %d, %f, %f, %f\n", dsvp_prime, dsvp, G_cum, PSC, cum_pr + (1-cum_pr) * psvp);

        if(cum_pr + (1-cum_pr) * psvp >= cost->params->succ_prob){
            if(cost->params->print_Gcums){
                printf("Pcums: ");
                print_vector(pcums);
                printf("Gcums: ");
                print_vector(Gcums);
            }
            // printf("dsvp = %d, dsvp_prime = %d, gates =  %f\n ", dsvp, dsvp_prime, p_cost.first);
            return  make_tuple(dsvp_prime, dsvp, G_cum,B_cum, PSC); 
        }
        pre_psvp = psvp;
        // pre_psvp2 = psvp2;
    }
    p_cost = cost->pump_cost(d,cost_model);
    return  make_tuple(dsvp_prime, d, G_cum,B_cum, PSC); 
}






//using the determing condition sqrt(dsvp/d)*GH(L)<= 4/3GH(L_[d-dsvp:])
tuple<int,int,double,double,double> average_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model, pair<double,double> cum_GB_BKZ ){
    /*
    return dsvp, G, B. 
    */
    int d = l.size(), f = 0 , dsvp; //, dsvp_prime_;
    double psvp, pre_psvp = 0.,  gh; //, avg_d_svp  = 0.; pre_psvp2 = 0., psvp2, 
    double G_cum = -1000., B_cum = -1000., PSC = -1000.;
    pair<double,double> p_cost ={-1000.,-1000.};
    for(dsvp = 30; dsvp <= d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
 
        if( (double) dsvp/d * gaussian_heuristic_log2(l,0) <= 4/3. * gaussian_heuristic_log2(l,d-dsvp))
            break;
    }

    p_cost = cost->pump_cost(dsvp,cost_model);
    G_cum = log2(pow(2,G_cum)+ (pow(2,p_cost.first) + pow(2,cum_GB_BKZ.first)) );
    B_cum = max(B_cum,p_cost.second) ;
    PSC = p_cost.first;
  
    return  make_tuple(dsvp, d, G_cum,B_cum, PSC); 
}






//asvp determin condition: GH(L[:dsvp])<=gamma GH and sqrt(dsvp/d)*GH(L[:svp])<= 4/3GH(L_[f:dsvp])
tuple<int,int,double,double,double> average_asvp_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model, pair<double,double> cum_GB_BKZ ){
    /*
    return dsvp, G, B. 
    */
    int d = l.size(), f, dsvp; //, dsvp_prime_;
    double psvp, pre_psvp = 0.,  gh; //, avg_d_svp  = 0.; pre_psvp2 = 0., psvp2, 
    double G_cum = -1000., B_cum = -1000., PSC = -1000.;
    pair<double,double> p_cost ={-1000.,-1000.};
    for(dsvp = 30; dsvp <= d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        for(f = dsvp-10; f>=0; --f){
            // cerr<<"f = "<<f<<endl;
            // cerr<<log2(gaussian_heuristic_log2(l,0,dsvp) ) /2. <<" "<< cost->params->log2_target_norm<<"----"<<endl;
            // throw "";
            if(cost->params->log2_target_norm == -1)
                throw "Foget to set target norm";
            else if( log2(gaussian_heuristic_log2(l,0,dsvp) ) /2. <= cost->params->log2_target_norm  -log2(1.05) &&  (double) dsvp/d * gaussian_heuristic_log2(l,0,dsvp) <= 4/3. * gaussian_heuristic_log2(l,f,dsvp))
                    break;
        }
    }

    p_cost = cost->pump_cost(dsvp,cost_model);
    G_cum = log2(pow(2,G_cum)+ (pow(2,p_cost.first) + pow(2,cum_GB_BKZ.first)) );
    B_cum = max(B_cum,p_cost.second) ;
    PSC = p_cost.first;
  
    return  make_tuple(dsvp, d, G_cum,B_cum, PSC); 
}

