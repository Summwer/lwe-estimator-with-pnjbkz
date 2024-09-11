// #include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include "framework/est.h"

using namespace std;
// using namespace boost;



double sim_strategy(Params* params, vector<double> l, vector<tuple<int,int,int>> strategy, double sigma){
    cout<<"cost_model = "<<params->cost_model<<endl;
    int dim = int(l.size());
    BKZJSim* sim = new BKZJSim(params,dim);
    COST* cost = new COST(params);
    double Gcum = -1000., Bcum = -1000., cum_pr = 0., rem_pr = 1., GBKZ = -1000.;
    pair<double,double> G;
    for(int i = 0; i<int(strategy.size()); i++){
        int beta = get<0>(strategy[i]), jump = get<1>(strategy[i]), tours = get<2>(strategy[i]);
        for(int t = 0; t< tours; t++){
            int beta_ = get_beta_(params, beta, jump, dim);
            // cout<<"beta = "<<beta <<", beta_ = "<< beta_<<". jump = "<< jump<<endl;

            // vector<double> l1, l2;
            sim -> simulate(l,l,beta,jump,1);
            double slope = get_current_slope(l,0,dim);
            boost::math::chi_squared chisquare(beta_);
            // cout<<beta_<<","<<beta<<endl;
            double pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
            // pr = 0.;
            // FP_NR<FT> pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
            // printf("(beta,jump) = (%d,%d), pr = %e\n", beta, jump, pr);
            // vector<double> l_;
            // sim -> simulate(l2,l,beta,1,1);
            // double slope1 = get_current_slope(l,0,dim);
            // boost::math::chi_squared chisquare1(beta);
            // // cout<<beta_<<","<<beta<<endl;
            
            // double pr1 = boost::math::cdf(chisquare1,pow(2,2.*l2[dim-beta]));
            // cout<<"slope = "<<slope1<<endl;
            // printf("(beta,jump) = (%d,%d), pr = %e\n", beta, 1, pr1);

            // l = l2;

            G = cost->bkz_cost(dim,beta,jump,params->cost_model);
            // cout<<"G.first"<<endl;
            if(params->verbose)
                printf("Strategy (%d,%d,%d), slope = %lf, sim-cost = %3.7f log(sec)\n", beta,jump,t+1,slope, G.first );
            
            GBKZ = log2(pow(2,GBKZ)+pow(2,G.first));
            // cout<<"GBKZ = "<<GBKZ<<", beta = "<<beta<<endl;
            if(not params->worst_case){
                Gcum = log2(pow(2,Gcum)+pow(2,GBKZ)*rem_pr*pr);
                Bcum = log2(pow(2,Bcum)+pow(2,G.second)*rem_pr*pr);
            }
            else{
                Gcum = GBKZ;
                Bcum = max(Bcum, G.second);
            }

            cum_pr += rem_pr * pr;
            rem_pr = 1. - cum_pr;

            // cout<<"cum_pr = "<<cum_pr<<endl;
        } 
    }

    // print_vector(l,0,dim);
    for(int i = 0; i < dim; i++){
        l[i] -= log2(sigma);
    }


    tuple<int,int,double,double,double> dsvp_t_;
    if(params->worst_case)
        dsvp_t_ = dsvp_predict(l, 0., cost, params->cost_model, make_pair(Gcum, Bcum));
    else
        dsvp_t_ = dsvp_predict(l, cum_pr, cost, params->cost_model, make_pair(GBKZ, G.second));
        // cout<<"cum_pr = "<<cum_pr<<endl;
    int dsvp = get<1>(dsvp_t_);
    // int f = wrapper_default_dim4free_fun(dsvp);
    int f = get_f_for_pump(params,dsvp);
    // int f = dims4free(dsvp);
    if(params->cost_model==1)
        printf("pump-{%d,%d,%d}, sim-pump cost = %3.7f sec\n",  dim - dsvp, dsvp, f,  get<2>(dsvp_t_)); 
        
    if(params->cost_model==2)
        printf("pump-{%d,%d,%d}, sim-pump cost = %3.7f sec\n",  dim - dsvp, dsvp, f,  pow(2,get<4>(dsvp_t_))); 

    if(not params->worst_case)
        Gcum = log2(pow(2,Gcum)+pow(2,get<2>(dsvp_t_)));
    else
        Gcum = log2(pow(2,Gcum)+pow(2,get<4>(dsvp_t_)));
         
    // Gcum = log2(pow(2,Gcum)+pow(2,get<4>(dsvp_t_)));
    
    Bcum = max(Bcum, get<3>(dsvp_t_));
    if(params->cost_model==2)
        cout<<"Gcum = "<<pow(2,Gcum)<<", Bcum = "<<pow(2,Bcum)<<endl;
    if(params->cost_model==1)
        cout<<"Gcum = "<< Gcum <<", Bcum = "<< Bcum<<endl;
    cout<<"============================="<<endl;

    return Gcum;
}


//Simulate the stratey from original lwe instance 
double test_lwechal_from_original_instance(Params* params, int n, double alpha, vector<tuple<int,int,int>> strategy){
    LWEchal* lwechal = gen_lwechal_instance(n, alpha);
    int dim = lwechal->dim;
    FP_NR<FT> dvol = lwechal->dvol;
    vector<double> l = lwechal->log_rr, l_;
    double  sigma = lwechal->alpha * lwechal->q;
    // printf("No sigma normalization,");
    // sim_strategy(params, l, strategy,sigma);

    printf("After a sigma normalization,");
    for(int i = 0; i < dim; i++){
        l[i] -=  log2(sigma);
    }
    double slope = get_current_slope(l,0,dim);
    printf("slope = %f\n", slope);
    return sim_strategy(params, l, strategy,1.);
}


//Simulate the stratey from gsa-gs-lengths and original lwe instance 
double test_lwechal_from_gsa(Params* params, int dim, double dvol, vector<tuple<int,int,int>> strategy){
    printf("Generate gs-lengths by GSA assumption...\n");
    vector<double>  l = gen_simulated_gso(dim, dvol);
    double slope = get_current_slope(l,0,dim);
    cout<<"Slope of gs-lengths generated by GSA assumption: "<<slope<<endl;

    return sim_strategy(params, l, strategy, 1.);
}


double test_nist_from_gsa(Params* params,int n, int m, int q,  map<int,double> D_e, map<int,double> D_s, vector<tuple<int,int,int>> strategy){
    printf("Generate gs-lengths by GSA assumption...\n");
    LWEchal* lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params->verbose);
    vector<double>  l = gen_simulated_gso(lwechal->dim, lwechal->dvol);
    double slope = get_current_slope(l,0,lwechal->dim);
    cout<<"Slope of gs-lengths generated by GSA assumption: "<<slope<<endl;

    return sim_strategy(params, l, strategy, 1.);
}


int main(){

    vector<tuple<int,int,int>> strategy;
    map<int,double> D_e, D_s;
    int n, m , q, eta;
    double alpha;
    Params* params = new Params;

    params->worst_case = false;
    params->cost_model = 1;
    params->theo_pnjbkz_d4f = 2;
    params->theo_pump_d4f = 2;
    params->print_Gcums = true;
    params->list_decoding = "matzov22";
    

    
    n = 40, alpha = 0.025;
    strategy = {{50,1,1},{51,1,1},{52,1,1},{53,1,1},{54,1,1},{55,1,1},{56,1,1},{57,1,1},{58,1,1},{59,1,1},{60,1,1},{61,1,1},{62,1,1},{63,1,1},{64,1,1},{65,1,1},{66,1,1},{67,1,1},{68,1,1},{69,1,1},{70,1,1},{71,1,1},{72,1,1},{73,1,1},{74,1,1},{75,1,1},{76,1,1},{77,1,1},{78,1,1},{87,4,1},{92,4,1},{98,4,1}};
    printf("n = %d, alpha = %1.3f, list_decoding = matzov22\n", n, alpha);
    int dim = ;
    double dvol =  
    test_lwechal_from_gsa(params, n, alpha, strategy);

   

    n = 80, alpha = 0.005;
    strategy = {{50,1,1},{51,1,1},{52,1,1},{53,1,1},{54,1,1},{55,1,1},{56,1,1},{57,1,1},{58,1,1},{59,1,1},{60,1,1},{61,1,1},{62,1,1},{63,1,1},{64,1,1},{65,1,1},{66,1,1},{67,1,1},{78,1,5},{85,4,1},{88,4,1},{89,4,1},{97,4,2},{100,4,1},{103,4,1},{107,4,1},{111,4,1}};
    printf("n = %d, alpha = %1.3f, list_decoding = matzov22\n", n, alpha);
    test_lwechal_from_gsa(params, n, alpha, strategy);



}