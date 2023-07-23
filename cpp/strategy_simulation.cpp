// #include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include "framework/est.h"

using namespace std;
// using namespace boost;



void sim_strategy(Params* params, vector<double> l, vector<tuple<int,int,int>> strategy, double sigma){
    cout<<"cost_model = "<<params->cost_model<<endl;
    int dim = int(l.size());
    BKZJSim* sim = new BKZJSim(params,dim);
    COST* cost = new COST(params);
    double Gcum = 0., Bcum = 0., cum_pr = 0., rem_pr = 1., GBKZ = 0.;
    for(int i = 0; i<int(strategy.size()); i++){
        int beta = get<0>(strategy[i]), jump = get<1>(strategy[i]), tours = get<2>(strategy[i]);
        for(int t = 0; t< tours; t++){
            sim -> simulate(l,l,beta,jump,1);
            int beta_ = get_beta_(params, beta, jump, dim);
            double slope = get_current_slope(l,0,dim);
            boost::math::chi_squared chisquare(beta_);
            // cout<<beta_<<","<<beta<<endl;
            double pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
            // FP_NR<FT> pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
            
            pair<double,double> G = cost->bkz_cost(dim,beta,jump,params->cost_model);
            if(params->verbose)
                printf("Strategy (%d,%d,%d}, slope = %lf, sim-cost = %3.7f sec\n", beta,jump,t+1,slope, pow(2,G.first));
            
            GBKZ = log2(pow(2,GBKZ)+pow(2,G.first));
            if(not params->worst_case)
                Gcum = log2(pow(2,Gcum)+pow(2,GBKZ)*rem_pr*pr);
            else
                Gcum = log2(pow(2,Gcum)+pow(2,G.first));
            Bcum = max(Bcum, G.second);

            cum_pr += rem_pr * pr;
            rem_pr = 1. - cum_pr;

            cout<<beta_<<","<<pr<<","<<l[dim-beta_]<<","<<cum_pr<<endl;
        } 
    }

    // print_vector(l,0,dim);
    for(int i = 0; i < dim; i++){
        l[i] -= log2(sigma);
    }

    tuple<int,int,double,double> dsvp_t_ = dsvp_predict(l, 0., cost, params->cost_model, params->progressive_sieve);

    int dsvp = get<1>(dsvp_t_);
    // int f = default_dim4free_fun(dsvp);
    int f = get_f_for_pump(params,dsvp);
    // int f = dims4free(dsvp);
    printf("pump-{%d,%d,%d}, sim-pump cost = %3.7f sec\n",  dim - dsvp, dsvp, f,  pow(2,get<2>(dsvp_t_)));  
    Gcum = log2(pow(2,Gcum)+pow(2,get<2>(dsvp_t_)));
    Bcum = max(Bcum, get<3>(dsvp_t_));
    cout<<"Gcum = "<<Gcum<<", Bcum = "<<Bcum<<endl;

    cout<<"============================="<<endl;
}


//Simulate the stratey from original lwe instance 
void test_lwechal_from_original_instance(Params* params, int n, double alpha, vector<tuple<int,int,int>> strategy){

   
    LWEchal* lwechal = gen_lwechal_instance(n, alpha);
    int dim = lwechal->dim;
    FP_NR<FT> dvol = lwechal->dvol;
    vector<double> l = lwechal->log_rr, l_;
    double  sigma = lwechal->alpha * lwechal->q;
    // printf("No sigma normalization,");
    // sim_strategy(l, strategy,sigma);


    printf("After a sigma normalization,");
    for(int i = 0; i < dim; i++){
        l[i] -=  log2(sigma);
    }
    double slope = get_current_slope(l,0,dim);
    printf("slope = %f\n", slope);
    
    sim_strategy(params, l, strategy,1.);

    // throw "";
}


//Simulate the stratey from gsa-gs-lengths and original lwe instance 
void test_lwechal_from_gsa(Params* params, int dim, double dvol, vector<tuple<int,int,int>> strategy){
    printf("Generate gs-lengths by GSA assumption...\n");
    vector<double>  l = gen_simulated_gso(dim, dvol);
    double slope = get_current_slope(l,0,dim);
    cout<<"Slope of gs-lengths generated by GSA assumption: "<<slope<<endl;

    sim_strategy(params, l, strategy, 1.);
}


int main(){

    // int n = 40;
    // double alpha = 0.035;
    // vector<tuple<int,int,int>> strategy = {{83, 8, 1}, {93, 8, 1}, {108, 8, 1}, {117, 8, 1}, {119, 4, 1}, {133, 4, 1}};
    // // for(int i = 10; i < 50; i++)
    // //     strategy.insert(strategy.end(},{i,1,1});
    // test_lwechal_from_original_instance(params, n, alpha, strategy);



    // Params* params = new Params;
    // params->cost_model = 2;
    // params->practical_pnjbkz_d4f = 3;
    // params->practical_pump_d4f = 2;
    // params->worst_case = true;
    // params->verbose = true;

    // int n = 45;
    // double alpha = 0.025;
    // vector<tuple<int,int,int>> strategy = {{106, 11, 1}, {116, 12, 1}, {117, 4, 1}, {120, 4, 1}};
    // test_lwechal_from_original_instance(params, n, alpha, strategy);


    Params* params = new Params;
    params->cost_model = 1;
    // params->theo_pnjbkz_d4f = 2;
    // params->theo_pump_d4f = 2;
    params->worst_case = false;
    // int dim =  2049;
    // double dvol = 15614.219317244602;
    // vector<tuple<int,int,int>> strategy;
    // strategy.resize(0);
    // for(int i = 50; i<=427; i++)
    //     strategy.insert(strategy.end(),{i,1,1});
    // printf("============= Dilithium-I\n");
    // test_lwechal_from_gsa(params, dim, dvol, strategy);



    int dim =  3540;
    double dvol = 26623.1162463;
    vector<tuple<int,int,int>> strategy;
    strategy.resize(0);
    for(int i = 50; i<=882; i++)
        strategy.insert(strategy.end(),{i,1,1});
    printf("============= Dilithium-III\n");
    test_lwechal_from_gsa(params, dim, dvol, strategy);

}