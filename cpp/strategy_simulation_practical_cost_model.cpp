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
    double Gcum = -1000., Bcum = -1000., cum_pr = 0., rem_pr = 1., GBKZ = -1000.;
    pair<double,double> G;

    for(int i = 0; i<int(strategy.size()); i++){
        int beta = get<0>(strategy[i]), jump = get<1>(strategy[i]), tours = get<2>(strategy[i]);
        for(int t = 0; t< tours; t++){
            int beta_ = get_beta_(params, beta, jump, dim);
            sim -> simulate(l,l,beta,jump,1);

            // print_vector(l,0,l.size());

            double slope = get_current_slope(l,0,dim);
            boost::math::chi_squared chisquare(beta_);
            // cout<<beta_<<","<<beta<<endl;
            double pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
            // FP_NR<FT> pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
        
            G = cost->bkz_cost(dim,beta,jump,params->cost_model);
            // printf("beta = %d", beta);
            // cout<<"G.first: "<<G.first<<endl;
            // print_vector(l);
            
            int jub  = jump_upper_bound(params,beta,l);

            // printf("d4f upper bound for beta  = %d is %d\n", beta, jub);

            GBKZ = log2(pow(2,GBKZ)+pow(2,G.first));
            // printf("Gbkz-one: %f", );
            // cout<<"GBKZ = "<<GBKZ<<", beta = "<<beta<<endl;
            if(not params->worst_case){
                Gcum = log2(pow(2,Gcum)+pow(2,GBKZ)*rem_pr*pr);
                Bcum = log2(pow(2,Bcum)+pow(2,G.second)*rem_pr*pr);
            }
            else{
                Gcum = GBKZ;
                Bcum = max(Bcum, G.second);
            }

            if(params->verbose)
                // printf("Strategy (%d,%d,%d), slope = %lf, pr = %e, l[dim-beta_] = %f, sim-cost = %3.7f log(sec), G_cum = %3.7f log(sec), G_BKZ = %3.7f log(sec), l[0] = %f\n", beta,jump,t+1,slope, pr, l[dim-beta_],  G.first, Gcum, GBKZ , l[0] );
                printf("(%d,%d,%d) %lf %lf %d\n", beta,jump,t+1,slope,  GBKZ , beta - get_f_for_pnjbkz(params, beta));

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
        cout<<"cum_pr = "<<cum_pr<<endl;
    int dsvp = get<1>(dsvp_t_);
    // int f = wrapper_default_dim4free_fun(dsvp);
    int f = get_f_for_pump(params,dsvp);
    // int f = dims4free(dsvp);
    if(params->cost_model==1)
        printf("pump-{%d,%d,%d}, Tdsvp = %3.7f sec, PSC = %3.7f sec\n",  dim - dsvp, dsvp, f,  get<2>(dsvp_t_), get<4>(dsvp_t_)); 
        
    if(params->cost_model==2)
        printf("pump-{%d,%d,%d}, PSC = %3.7f sec\n",  dim - dsvp, dsvp, f,  get<4>(dsvp_t_)); 

    if(not params->worst_case)
        Gcum = log2(pow(2,Gcum)+pow(2,get<2>(dsvp_t_)));
    else
        Gcum = log2(pow(2,Gcum)+pow(2,get<4>(dsvp_t_)));
         
    // Gcum = log2(pow(2,Gcum)+pow(2,get<4>(dsvp_t_)));
    // printf("(slope = %e, G_BKZ = %e log2(gate), B_BKZ = %e log2(bit), cum-pr = %e, pump-{%d,%d,%d}, G_dsvp = %e log2(gate), B_dsvp = %e bit, avgG = %e log2(gate), avgB = %e log2(bit),  G = %e log2(gate), min_GB.first = %e log2(gate), leaf = %d)\n",  bs.slope, bs.cum_GB_BKZ.first, bs.cum_GB_BKZ.second, bs.cum_pr,  d- get<1>(bs.dsvp_t), get<1>(bs.dsvp_t), get<1>(bs.dsvp_t) - get<0>(bs.dsvp_t), get<2>(bs.dsvp_t),get<3>(bs.dsvp_t), bs.avg_GB.first, bs.avg_GB.second, bs.GB_nopr.first, bs.min_GB.first,bs.leaf);  
    cout<<"S(beta,jump,tours):[";
    for(int i = 0; i < int(strategy.size()); i ++){
        printf("(%d,%d,%d)",get<0>(strategy[i]),get<1>(strategy[i]),get<2>(strategy[i]));
        if(i!=int(strategy.size()) - 1)
            printf(",");
    }
    cout<<"]"<<endl;

    Bcum = max(Bcum, get<3>(dsvp_t_));
    if(params->cost_model==2)
        cout<<"Gcum = "<<pow(2,Gcum)<<", Bcum = "<<pow(2,Bcum)<<endl;
    if(params->cost_model==1)
        cout<<"Gcum = "<< Gcum <<", Bcum = "<< Bcum<<endl;
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
    // sim_strategy(params, l, strategy,sigma);

    printf("After a sigma normalization,");
    for(int i = 0; i < dim; i++){
        l[i] -=  log2(sigma);
    }
    // print_vector(l,0,dim);
    double slope = get_current_slope(l,0,dim);
    printf("slope = %f\n", slope);
    sim_strategy(params, l, strategy,1.);
}


//Simulate the stratey from gsa-gs-lengths and original lwe instance 
void test_lwechal_from_gsa(Params* params, int dim, double dvol, vector<tuple<int,int,int>> strategy){
    printf("Generate gs-lengths by GSA assumption...\n");
    vector<double>  l = gen_simulated_gso(dim, dvol);
    double slope = get_current_slope(l,0,dim);
    cout<<"Slope of gs-lengths generated by GSA assumption: "<<slope<<endl;

    sim_strategy(params, l, strategy, 1.);
}


void test_nist_from_gsa(Params* params,int n, int m, int q,  map<int,double> D_e, map<int,double> D_s, vector<tuple<int,int,int>> strategy){
    printf("Generate gs-lengths by GSA assumption...\n");
    LWEchal* lwechal = gen_LWE_instance_with_input_distribution( n, q, m, D_e, D_s, params->verbose);
    vector<double>  l = gen_simulated_gso(lwechal->dim, lwechal->dvol);
    double slope = get_current_slope(l,0,lwechal->dim);
    cout<<"Slope of gs-lengths generated by GSA assumption: "<<slope<<endl;

    sim_strategy(params, l, strategy, 1.);
}


int main(){
    map<int,double> D_e, D_s;
    int n, m , q, eta;
    double alpha;
    Params* params = new Params;
    params->cost_model = 2;
    params->practical_pnjbkz_d4f = 3;
    params->practical_pump_d4f = 3;
    params->worst_case = true;
    params->verbose = true;
    vector<tuple<int,int,int>> strategy;

    int dim=222;
    double dvol = 522.7192362;
    n = 75;
    alpha = 0.010;

    strategy =  {{76,8,1},{89,9,1},{92,9,1},{115,10,1},{117,10,1},{117,10,1},{117,4,1},{118,4,1},{128,4,1},{132,4,1},{141,4,1},{144,4,1},{150,4,1},{155,4,1}, {157,4,1}, {159,4,1}, {160,1,1},  {165,1,1}, {170,2,1},  {173,1,1}};
   
    test_lwechal_from_original_instance(params, n, alpha, strategy);


    n = 100, alpha = 0.005;
    params->worst_case = false;
    strategy =  {{89,9, 1}, {90,9, 1}, {93,9, 1}, {96,9, 1}, {117, 10, 1}, {117, 10, 1}, {117,10, 1}, {119,10, 1}, {117,4,1}, {118,4, 1}, {121,4, 1}, {128,4, 1}, {131,4, 1}, {135,4, 1}, {141,4, 1}, {141,4, 1},{145,4, 1},{148,4, 1}, {149,2, 1}, {155,4,1}, {156,2,1}, {160,2,1}, {167,2,1}, {170,2,1}, {172,2,1}, {174,2,1}}; 
    test_lwechal_from_original_instance(params, n, alpha, strategy);
    
}