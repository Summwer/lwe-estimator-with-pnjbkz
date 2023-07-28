#include "d_svp_prediction.h"


tuple<int,int,double,double> dsvp_predict(vector<double> l,  double cum_pr, COST* cost, int cost_model, bool progressive_sieve, bool worst_case){
    /*
    return dsvp, G, B. 
    */
    if(progressive_sieve){
        if(cost->params-> est_model == 1)
            return progressive_dsvp_predict1(l, cum_pr, cost, cost_model, worst_case);
        else if(cost->params-> est_model == 2)
            return progressive_dsvp_predict2(l, cum_pr, cost, cost_model);
    }
    return fixed_dsvp_predict(l, cum_pr, cost, cost_model);//, worst_case
}


tuple<int,int,double,double> fixed_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model, bool worst_case){
    /*
    return dsvp, G, B. 
    */
    int d = l.size();
    double psvp, p, rp, gh, dsvp_;
    pair<double,double> p_cost;
    if(cum_pr >= 0.999){
        return make_tuple(0.,0,0.,0.);
    }
    for(int dsvp = 50; dsvp <= d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic_log2(l,d-dsvp);
        
        
        boost::math::chi_squared chisquare(dsvp);
        psvp = boost::math::cdf(chisquare,gh); //Compute chi-squared value
        
        p = cum_pr + (1. - cum_pr)* psvp;
        rp = 1. - p;
        if(rp < 0.001){
            dsvp_ = dsvp + dsvp * rp;
            p_cost = cost->pump_cost(dsvp,cost_model);
            if(not worst_case)
                return  make_tuple(dsvp_, dsvp, p_cost.first * ((1-cum_pr) * psvp + rp), p_cost.second);  //Avoid too small of dsvp
            else
                return  make_tuple(dsvp, dsvp, p_cost.first, p_cost.second);
        }
    }
    p_cost = cost->pump_cost(d,cost_model);
    return make_tuple(d, d, p_cost.first, p_cost.second);   
}



//no cum G with model 1
tuple<int,int,double,double> progressive_dsvp_predict1(vector<double> l, double cum_pr, COST* cost, int cost_model, bool worst_case){
    /*
    return dsvp, G, B. 
    */
    int d = l.size();
    double psvp, pre_psvp = 0., p = cum_pr , rp = 1.-cum_pr, gh, avg_d_svp  = 0.;
    pair<double,double> p_cost, p_cost_cum;
    double G_cum = 0., B_cum = 0.;
    
    if(cum_pr >= 0.999){
        return make_tuple(0.,0,0.,0.);
    }
    
    for(int dsvp = 50; dsvp <= d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic_log2(l,d-dsvp);
        
        boost::math::chi_squared chisquare(dsvp);
        psvp = boost::math::cdf(chisquare,gh); //Compute chi-squared value

        int f = get_f_for_pump(cost->params,dsvp);

        int dsvp_prime = floor(dsvp - f);
        // int dsvp_prime = dsvp;
        // if(cost->params->cost_model == 2)
        //     dsvp_prime = dsvp - f;
        // if(cost->params->cost_model == 1)
        //     dsvp_prime = floor(dsvp - f);
        p_cost = cost->pump_cost(dsvp_prime,cost_model);
        if(not worst_case){
            avg_d_svp += dsvp * rp * psvp;
            G_cum = log2(pow(2,G_cum)+pow(2,p_cost.first) * rp * psvp);
            // B_cum = log2(pow(2,B_cum)+pow(2,p_cost.second) * rp * psvp);
            B_cum = p_cost.second;
        }else{
            // G_cum = p_cost.first;
            // B_cum = p_cost.second;
            G_cum = log2(pow(2,G_cum)+pow(2,p_cost.first) * (psvp-pre_psvp));
            B_cum = max(B_cum,p_cost.second);
        }
        p += rp * psvp;
        rp = 1. - p;
        // cerr<<"dsvp = "<<dsvp<<", rp = "<<rp<<endl;
        // cerr<<"G_cum = " << G_cum <<", G2 = "<<p_cost.first<<endl;
        if(not worst_case){
            if(rp < 0.001){
                avg_d_svp += dsvp * rp;
                G_cum = log2(pow(2,G_cum)+pow(2,p_cost.first) * rp);
                return  make_tuple(avg_d_svp, dsvp, G_cum,B_cum); 
                // return  make_tuple(round(avg_d_svp*PREC)/PREC, dsvp, round(G_cum*PREC)/PREC,B_cum); 
            }
        }else{
            if(1-psvp < 0.001){
                return  make_tuple(dsvp, dsvp, G_cum,B_cum); 
            }
        }
        pre_psvp = psvp;
    }
    p_cost = cost->pump_cost(d,cost_model);
    return make_tuple(d,d, p_cost.first, p_cost.second);
}



//no cum G with model 2
tuple<int,int,double,double> progressive_dsvp_predict2(vector<double> l, double cum_pr, COST* cost, int cost_model){ //, bool worst_case
    /*
    return dsvp, G, B. 
    */
    int d = l.size();
    double psvp1, psvp2, pre_psvp2 =0., gh;
    pair<double,double> p_cost;
    double G_cum = 0., B_cum = 0.;
    
    if(cum_pr >= 0.999){
        return make_tuple(0,0,0.,0.);
    }
    int dsvp = 50, dsvp1, dsvp2;
    bool flag1 = false, flag2 = false;

    while(dsvp <= d){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic_log2(l,d-dsvp);
        
        boost::math::chi_squared chisquare(dsvp);
        psvp1 = boost::math::cdf(chisquare,gh); //Compute chi-squared value
        psvp2 = boost::math::cdf(chisquare,4/3. * gh); 
        
        if(not flag2){
            p_cost = cost->pump_cost(dsvp,cost_model);
            G_cum = log2(pow(2,G_cum)+pow(2,p_cost.first) * (psvp2-pre_psvp2));
            B_cum = log2(pow(2,B_cum)+pow(2,p_cost.second) * (psvp2-pre_psvp2));
        }
       
        if(psvp1 > 0.999 and not flag1){
            dsvp1 = dsvp;
            flag1 = true;
        }
        if(psvp2 > 0.999 and not flag2){
            dsvp2 = dsvp;
            flag2 = true;
        }
        if(flag1 and flag2)
            return  make_tuple(dsvp1, dsvp2, G_cum,B_cum); 
    
        pre_psvp2 = psvp2;
        dsvp++;
    }
    p_cost = cost->pump_cost(d,cost_model);
    return make_tuple(d,d, p_cost.first, p_cost.second);
}