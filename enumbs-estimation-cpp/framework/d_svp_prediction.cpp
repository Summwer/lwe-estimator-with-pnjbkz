#include "d_svp_prediction.h"


tuple<double,int,double,double> dsvp_predict(vector<double> l,  double cum_pr, COST* cost, int cost_model, bool progressive_sieve){
    /*
    return dsvp, G, B. 
    */
    if(progressive_sieve)
        return progressive_dsvp_predict(l, cum_pr, cost, cost_model);
    else
        return fixed_dsvp_predict(l, cum_pr, cost, cost_model);
}

tuple<double,int,double,double> fixed_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model){
    /*
    return dsvp, G, B. 
    */
    int d = l.size();
    double psvp, p, rp, gh, dsvp_;
    pair<double,double> p_cost;
    if(cum_pr >= 1.){
        return make_tuple(0.,0,0.,0.);
    }
    for(int dsvp = 50; dsvp < d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic_log2(l,d-dsvp);
        
        boost::math::chi_squared chisquare(dsvp);
        psvp = boost::math::cdf(chisquare,gh); //Compute chi-squared value
        
        p = cum_pr + (1. - cum_pr)* psvp;
        // cout<<dsvp<<","<<gh<<","<<p<<endl;
        rp = 1. - p;
        if(rp < 0.001){
            dsvp_ = dsvp + dsvp * rp;
            p_cost = cost->pump_cost(dsvp_,cost_model);
            return  make_tuple(dsvp_, dsvp, p_cost.first, p_cost.second);  //Avoid too small of dsvp
            // return  make_tuple(round(dsvp_*PREC)/PREC, dsvp, round(p_cost.first*PREC)/PREC, p_cost.second);  //Avoid too small of dsvp
            // return  make_tuple(round(dsvp_*PREC)/PREC, dsvp, round(p_cost.first*PREC)/PREC, p_cost.second);  //Avoid too small of dsvp
            
        }
    }
    p_cost = cost->pump_cost(d,cost_model);
    return make_tuple(d, d, p_cost.first, p_cost.second);
}

tuple<double,int,double,double> progressive_dsvp_predict(vector<double> l, double cum_pr, COST* cost, int cost_model){
    /*
    return dsvp, G, B. 
    */
    int d = l.size();
    double psvp, p = cum_pr , rp = 1.-cum_pr, gh, avg_d_svp  = 0.;
    pair<double,double> p_cost, p_cost_cum;
    double G_cum = 0., B_cum = 0.;
    if(cum_pr >= 1.){
        return make_tuple(0.,0,0.,0.);
    }
    for(int dsvp = 50; dsvp < d; dsvp++ ){
        //2**(2 * l1[d-dsvp])==2**(2 * l1d_dsvp)==gh
        gh = gaussian_heuristic_log2(l,d-dsvp);
        
        boost::math::chi_squared chisquare(dsvp);
        psvp = boost::math::cdf(chisquare,gh); //Compute chi-squared value
        
        p += rp * psvp;
        avg_d_svp += dsvp * rp * psvp;

        p_cost = cost->pump_cost(dsvp,cost_model);
        G_cum = log2(pow(2,G_cum)+pow(2,p_cost.first) * rp * psvp);
        B_cum = max(B_cum,p_cost.second);
     
        rp = 1. - p;
        if(rp < 0.001){
            avg_d_svp += dsvp * rp;
            G_cum = log2(pow(2,G_cum)+pow(2,p_cost.first) * rp);

            return  make_tuple(avg_d_svp, dsvp, G_cum,B_cum);  //Avoid too small of dsvp
            // return  make_tuple(round(avg_d_svp*PREC)/PREC, dsvp, round(G_cum*PREC)/PREC,B_cum); 
        }
    }
    p_cost = cost->pump_cost(d,cost_model);
    return make_tuple(d,d, p_cost.first, p_cost.second);
}
