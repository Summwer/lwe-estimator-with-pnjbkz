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
            double slope = get_current_slope(l,0,dim);
            boost::math::chi_squared chisquare(beta_);
            // cout<<beta_<<","<<beta<<endl;
            double pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
            // FP_NR<FT> pr = boost::math::cdf(chisquare,pow(2,2.*l[dim-beta_]));
        
            G = cost->bkz_cost(dim,beta,jump,params->cost_model);
            // cout<<"G.first"<<endl;
            // print_vector(l);
            if(params->verbose)
                printf("Strategy (%d,%d,%d), slope = %lf, pr = %e, l[dim-beta_] = %f, sim-cost = %3.7f log(sec)\n", beta,jump,t+1,slope, pr, l[dim-beta_],  G.first );
            
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
        cout<<"cum_pr = "<<cum_pr<<endl;
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

    // vector<tuple<int,int,int>> strategy;
    // int n;
    // double alpha;
    // Params* params = new Params;
    // params->cost_model = 2;
    // params->practical_pnjbkz_d4f = 3;
    // params->practical_pump_d4f = 2;
    // params->worst_case = true;
    // params->verbose = true;

    // vector<tuple<int,int,int>> strategy = {{83, 8, 1}, {93, 8, 1}, {108, 8, 1}, {117, 8, 1}, {119, 4, 1}, {133, 4, 1}};
    // test_lwechal_from_original_instance(params, n, alpha, strategy);


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


    // n = 40, alpha = 0.030;
    // strategy = {{73,8,1},{89,9,1},{117,10,1},{119,10,1}};
    // test_lwechal_from_original_instance(params, n, alpha, strategy);



    // params->worst_case = false;
    // params->cost_model = 1;
    // params->theo_pnjbkz_d4f = 2;
    // params->theo_pump_d4f = 2;
    // params->list_decoding = "agps20"; //"matzov22";
    // printf("============= Kyber-1024\n");
    // n = 1024, m = 1024, q = 3329;
    // D_s = build_centered_binomial_law(2);
    // D_e = D_s;
    // // strategy = {};
    // // for(int i = 50; i<= 903; i++){
    // //     strategy.insert(strategy.end(),{i,1,1});
    // // }
    // // test_nist_from_gsa(params, n, m, q, D_e, D_s,strategy);


    // strategy = {{50,7,1},{51,7,1},{52,7,1},{53,7,1},{78,1,1},{79,7,1},{79,4,1},{80,7,1},{80,4,1},{81,7,1},{82,7,1},{82,4,1},{83,7,1},{83,4,1},{84,7,1},{84,4,1},{85,7,1},{85,4,1},{86,7,1},{86,4,1},{87,7,1},{87,4,1},{88,7,1},{88,4,1},{89,7,1},{89,4,1},{90,7,1},{90,4,1},{91,7,1},{91,4,1},{92,7,1},{92,4,1},{93,7,1},{93,4,1},{94,7,1},{94,4,1},{95,7,1},{95,4,1},{96,7,1},{97,7,1},{97,4,1},{98,7,1},{98,4,1},{99,7,1},{99,4,1},{100,7,1},{100,4,1},{101,7,1},{101,4,1},{102,7,1},{102,4,1},{103,7,1},{103,4,1},{104,7,1},{104,4,1},{105,7,1},{105,4,1},{106,7,1},{106,4,1},{107,7,1},{107,4,1},{108,7,1},{108,4,1},{109,7,1},{109,4,1},{110,7,1},{111,7,1},{111,4,1},{112,7,1},{112,4,1},{113,7,1},{113,4,1},{114,7,1},{114,4,1},{115,7,1},{115,4,1},{116,7,1},{116,4,1},{117,7,1},{117,4,1},{118,7,1},{118,4,1},{119,7,1},{119,4,1},{120,7,1},{120,4,1},{121,7,1},{121,4,1},{122,7,1},{122,4,1},{123,7,1},{123,4,1},{124,7,1},{125,7,1},{125,4,1},{126,7,1},{126,4,1},{127,7,1},{127,4,1},{128,7,1},{128,4,1},{129,7,1},{129,4,1},{130,7,1},{130,4,1},{131,7,1},{131,4,1},{132,7,1},{132,4,1},{133,7,1},{133,4,1},{134,7,1},{134,4,1},{135,7,1},{135,4,1},{136,7,1},{136,4,1},{138,10,1},{138,7,1},{138,4,1},{139,10,1},{139,7,1},{140,10,1},{140,4,1},{141,10,1},{142,10,1},{142,4,1},{143,10,1},{144,10,1},{144,4,1},{145,10,1},{146,10,1},{146,4,1},{147,10,1},{148,10,1},{148,7,1},{148,4,1},{149,10,1},{149,7,1},{149,4,1},{150,10,1},{150,7,1},{150,4,1},{151,10,1},{151,7,1},{152,7,1},{152,4,1},{153,10,1},{153,7,1},{153,4,1},{154,10,1},{154,7,1},{154,4,1},{155,10,1},{155,7,1},{155,4,1},{156,10,1},{156,7,1},{156,4,1},{157,10,1},{157,7,1},{157,4,1},{158,10,1},{158,7,1},{158,4,1},{159,10,1},{159,7,1},{159,4,1},{160,10,1},{160,7,1},{160,4,1},{161,10,1},{161,7,1},{161,4,1},{162,10,1},{162,7,1},{162,4,1},{164,10,1},{164,7,1},{165,10,1},{166,10,1},{167,10,1},{168,10,1},{169,10,1},{170,10,1},{171,10,1},{172,10,1},{173,10,1},{173,7,1},{173,4,1},{174,10,1},{175,10,1},{175,7,1},{176,10,1},{177,10,1},{178,10,1},{179,10,1},{180,10,1},{181,10,1},{182,10,1},{182,7,1},{183,10,1},{184,10,1},{185,10,1},{186,10,1},{186,7,1},{187,10,1},{188,10,1},{189,10,1},{190,10,1},{191,10,1},{191,7,1},{193,10,1},{194,10,1},{195,10,1},{196,10,1},{197,10,1},{198,10,1},{198,7,1},{199,10,1},{199,7,1},{201,10,1},{203,10,1},{204,10,1},{205,10,1},{206,10,1},{207,10,1},{208,10,1},{209,10,1},{210,10,1},{210,7,1},{210,4,1},{211,4,1},{212,4,1},{214,4,1},{215,4,1},{216,4,1},{217,4,1},{218,4,1},{219,4,1},{220,4,1},{221,4,1},{222,4,1},{224,4,1},{225,4,1},{226,4,1},{227,4,1},{228,4,1},{229,4,1},{230,4,1},{231,4,1},{232,4,1},{233,4,1},{235,4,1},{236,4,1},{239,4,1},{240,4,1},{241,4,1},{242,4,1},{243,4,1},{245,4,1},{246,4,1},{248,4,1},{250,4,1},{251,4,1},{253,4,1},{254,4,1},{256,4,1},{257,4,1},{259,4,1},{261,4,1},{262,4,1},{264,4,1},{266,4,1},{268,4,1},{269,4,1},{271,4,1},{273,4,1},{275,4,1},{277,4,1},{279,4,1},{280,4,1},{284,4,1},{285,4,1},{286,4,1},{288,4,1},{291,4,1},{292,4,1},{294,4,1},{296,4,1},{299,4,1},{301,4,1},{303,4,1},{305,4,1},{307,4,1},{309,4,1},{312,4,1},{314,4,1},{316,4,1},{319,4,1},{321,4,1},{324,4,1},{326,4,1},{328,4,1},{331,4,1},{334,4,1},{336,4,1},{339,4,1},{342,4,1},{345,4,1},{347,4,1},{350,4,1},{353,4,1},{355,4,1},{359,4,1},{362,4,1},{365,4,1},{368,4,1},{371,4,1},{374,4,1},{377,4,1},{380,4,1},{383,4,1},{387,4,1},{390,4,1},{394,4,1},{397,4,1},{400,4,1},{404,4,1},{407,4,1},{411,4,1},{415,4,1},{419,4,1},{422,4,1},{427,4,1},{430,4,1},{434,4,1},{439,4,1},{443,4,1},{447,4,1},{451,4,1},{455,4,1},{460,4,1},{465,4,1},{469,4,1},{474,4,1},{479,4,1},{483,4,1},{489,4,1},{493,4,1},{499,4,1},{504,4,1},{509,4,1},{514,4,1},{520,4,1},{526,4,1},{532,4,1},{538,4,1},{543,4,1},{549,4,1},{555,4,1},{562,4,1},{568,4,1},{575,4,1},{581,4,1},{589,4,1},{595,4,1},{603,4,1},{610,4,1},{617,4,1},{626,4,1},{634,4,1},{641,4,1},{649,4,1},{658,4,1},{666,4,1},{676,4,1},{685,4,1},{694,4,1},{703,4,1},{714,4,1},{723,4,1},{734,4,1},{745,4,1},{755,4,1},{766,4,1},{778,4,1},{789,4,1},{801,4,1},{814,4,1},{826,4,1}};
    // test_nist_from_gsa(params, n, m, q, D_e, D_s,strategy);




    params->worst_case = false;
    params->cost_model = 1;
    params->theo_pnjbkz_d4f = 2;
    params->theo_pump_d4f = 2;
    params->list_decoding = "matzov22";
    printf("============= Dilithium-I\n");
    n = 4*256, m = 4*256, q = 8380417, eta = 2;
     D_s={},D_e={};
    for(int x=-eta; x<=eta; x++){
        D_s[x] = 1./(2*eta+1);
        D_e[x] = 1./(2*eta+1);
    }

    // strategy = {{50,1,1},{51,1,1},{52,1,1},{53,1,1},{54,1,1},{55,1,1},{56,1,1},{57,1,1},{58,1,1},{59,1,1},{60,1,1},{61,1,1},{62,1,1},{63,1,1},{64,1,1},{65,1,1},{66,1,1},{67,1,1},{68,1,1},{69,1,1},{70,1,1},{71,1,1},{72,1,1},{73,1,1},{74,1,1},{75,1,1},{76,1,1},{77,1,1},{78,1,1},{79,1,1},{80,1,1},{81,1,1},{82,1,1},{83,1,1},{84,1,1},{85,1,1},{86,1,1},{87,1,1},{88,1,1},{89,1,1},{90,1,1},{91,1,1},{92,1,1},{93,1,1},{94,1,1},{95,1,1},{96,1,1},{97,1,1},{98,1,1},{99,1,1},{100,1,1},{101,1,1},{102,1,1},{103,1,1},{104,1,1},{105,1,1},{106,1,1},{107,1,1},{108,1,1},{109,1,1},{110,1,1},{111,1,1},{112,1,1},{113,1,1},{114,1,1},{115,1,1},{116,1,1},{117,1,1},{118,1,1},{119,1,1},{120,1,1},{121,1,1},{122,1,1},{123,1,1},{124,1,1},{125,1,1},{126,1,1},{127,1,1},{128,1,1},{129,1,1},{130,1,1},{131,1,1},{132,1,1},{133,1,1},{134,1,1},{135,1,1},{136,1,1},{137,1,1},{138,1,1},{139,1,1},{140,1,1},{141,1,1},{142,1,1},{143,1,1},{144,1,1},{145,1,1},{146,1,1},{147,1,1},{148,1,1},{149,1,1},{150,1,1},{151,1,1},{152,1,1},{153,1,1},{154,1,1},{155,1,1},{156,1,1},{157,1,1},{158,1,1},{159,1,1},{160,1,1},{161,1,1},{162,1,1},{163,1,1},{164,1,1},{165,1,1},{166,1,1},{167,1,1},{168,1,1},{169,1,1},{170,1,1},{171,1,1},{172,1,1},{173,1,1},{174,1,1},{175,1,1},{176,1,1},{177,1,1},{178,1,1},{179,1,1},{180,1,1},{181,1,1},{182,1,1},{183,1,1},{184,1,1},{185,1,1},{186,1,1},{187,1,1},{188,1,1},{189,1,1},{190,1,1},{191,1,1},{192,1,1},{193,1,1},{194,1,1},{195,1,1},{196,1,1},{197,1,1},{198,1,1},{199,1,1},{200,1,1},{201,1,1},{202,1,1},{203,1,1},{204,1,1},{205,1,1},{206,1,1},{207,1,1},{208,1,1},{209,1,1},{210,1,1},{211,1,1},{212,1,1},{213,1,1},{214,1,1},{215,1,1},{216,1,1},{217,1,1},{218,1,1},{219,1,1},{220,1,1},{221,1,1},{222,1,1},{223,1,1},{224,1,1},{225,1,1},{226,1,1},{227,1,1},{228,1,1},{229,1,1},{230,1,1},{231,1,1},{232,1,1},{233,1,1},{234,1,1},{235,1,1},{236,1,1},{237,1,1},{238,1,1},{239,1,1},{240,1,1},{241,1,1},{242,1,1},{243,1,1},{244,1,1},{245,1,1},{246,1,1},{247,1,1},{248,1,1},{249,1,1},{250,1,1},{251,1,1},{252,1,1},{253,1,1},{254,1,1},{255,1,1},{256,1,1},{257,1,1},{258,1,1},{259,1,1},{260,1,1},{261,1,1},{262,1,1},{263,1,1},{264,1,1},{265,1,1},{266,1,1},{267,1,1},{268,1,1},{269,1,1},{270,1,1},{271,1,1},{272,1,1},{273,1,1},{274,1,1},{275,1,1},{276,1,1},{277,1,1},{278,1,1},{279,1,1},{280,1,1},{281,1,1},{282,1,1},{283,1,1},{284,1,1},{285,1,1},{286,1,1},{287,1,1},{288,1,1},{289,1,1},{290,1,1},{291,1,1},{292,1,1},{293,1,1},{294,1,1},{295,1,1},{296,1,1},{297,1,1},{298,1,1},{299,1,1},{300,1,1},{301,1,1},{302,1,1},{303,1,1},{304,1,1},{305,1,1},{306,1,1},{307,1,1},{308,1,1},{309,1,1},{310,1,1},{311,1,1},{312,1,1},{313,1,1},{314,1,1},{315,1,1},{316,1,1},{317,1,1},{318,1,1},{319,1,1},{320,1,1},{321,1,1},{322,1,1},{323,1,1},{324,1,1},{325,1,1},{326,1,1},{327,1,1},{328,1,1},{329,1,1},{330,1,1},{331,1,1},{332,1,1},{333,1,1},{334,1,1},{335,1,1},{336,1,1},{337,1,1},{338,1,1},{339,1,1},{340,1,1},{341,1,1},{342,1,1},{343,1,1},{344,1,1},{345,1,1},{346,1,1},{347,1,1},{348,1,1},{349,1,1},{350,1,1},{351,1,1},{352,1,1},{353,1,1},{354,1,1},{355,1,1},{356,1,1},{357,1,1},{358,1,1},{359,1,1},{360,1,1},{361,1,1},{362,1,1},{363,1,1},{364,1,1},{365,1,1},{366,1,1},{367,1,1},{368,1,1},{369,1,1},{370,1,1},{371,1,1},{372,1,1},{373,1,1},{374,1,1},{375,1,1},{376,1,1},{377,1,1},{378,1,1},{379,1,1},{380,1,1},{381,1,1},{382,1,1},{383,1,1},{384,1,1},{385,1,1},{386,1,1},{387,1,1},{388,1,1},{389,1,1},{390,1,1},{391,1,1},{392,1,1},{410,4,1},{414,4,1},{417,4,1},{421,4,1},{425,4,1},{429,4,1},{432,4,1},{436,4,1}};

    strategy = {{50,1,3},{51,1,3},{52,1,3},{53,1,3},{54,1,3},{55,1,3},{56,1,3},{57,1,3},{58,1,3},{59,1,3},{60,1,3},{61,1,3},{62,1,3},{63,1,3},{64,1,3},{65,1,3},{66,1,3},{67,1,3},{68,1,3},{69,1,3},{70,1,3},{71,1,3},{72,1,3},{73,1,3},{74,1,3},{75,1,3},{76,1,3},{77,1,3},{78,1,3},{79,1,3},{80,1,3},{81,1,3},{82,1,3},{83,1,3},{84,1,3},{85,1,3},{86,1,3},{87,1,3},{88,1,3},{89,1,3},{90,1,3},{91,1,3},{92,1,3},{93,1,3},{94,1,3},{95,1,3},{96,1,3},{97,1,3},{98,1,3},{99,1,3},{100,1,3},{101,1,3},{102,1,3},{103,1,3},{104,1,3},{105,1,3},{106,1,3},{107,1,3},{108,1,3},{109,1,3},{110,1,3},{111,1,3},{112,1,3},{113,1,3},{114,1,3},{115,1,3},{116,1,3},{117,1,3},{118,1,3},{119,1,3},{120,1,3},{121,1,3},{122,1,3},{123,1,3},{124,1,3},{125,1,3},{126,1,3},{127,1,3},{128,1,3},{129,1,3},{130,1,3},{131,1,3},{132,1,3},{133,1,3},{134,1,3},{135,1,3},{136,1,3},{137,1,3},{138,1,3},{139,1,3},{140,1,3},{141,1,3},{142,1,3},{143,1,3},{144,1,3},{145,1,3},{146,1,3},{147,1,3},{148,1,3},{149,1,3},{150,1,3},{151,1,3},{152,1,3},{153,1,3},{154,1,3},{155,1,3},{156,1,3},{157,1,3},{158,1,3},{159,1,3},{160,1,3},{161,1,3},{162,1,3},{163,1,3},{164,1,3},{165,1,3},{166,1,3},{167,1,3},{168,1,3},{169,1,3},{170,1,3},{171,1,3},{172,1,3},{173,1,3},{174,1,3},{175,1,3},{176,1,3},{177,1,3},{178,1,3},{179,1,3},{180,1,3},{181,1,3},{182,1,3},{183,1,3},{184,1,3},{185,1,3},{186,1,3},{187,1,3},{188,1,3},{189,1,3},{190,1,3},{191,1,3},{192,1,3},{193,1,3},{194,1,3},{195,1,3},{196,1,3},{197,1,3},{198,1,3},{199,1,3},{200,1,3},{201,1,3},{202,1,3},{203,1,3},{204,1,3},{205,1,3},{206,1,3},{207,1,3},{208,1,3},{209,1,3},{210,1,3},{211,1,3},{212,1,3},{213,1,3},{214,1,3},{215,1,3},{216,1,3},{217,1,3},{218,1,3},{219,1,3},{220,1,3},{221,1,3},{222,1,3},{223,1,3},{224,1,3},{225,1,3},{226,1,3},{227,1,3},{228,1,3},{229,1,3},{230,1,3},{231,1,3},{232,1,3},{233,1,3},{234,1,3},{235,1,3},{236,1,3},{237,1,3},{238,1,3},{239,1,3},{240,1,3},{241,1,3},{242,1,3},{243,1,3},{244,1,3},{245,1,3},{246,1,3},{247,1,3},{248,1,3},{249,1,3},{250,1,3},{251,1,3},{252,1,3},{253,1,3},{254,1,3},{255,1,3},{256,1,3},{257,1,3},{258,1,3},{259,1,3},{260,1,3},{261,1,3},{262,1,3},{263,1,3},{264,1,3},{265,1,3},{266,1,3},{267,1,3},{268,1,3},{269,1,3},{270,1,3},{271,1,3},{272,1,3},{273,1,3},{274,1,3},{275,1,3},{276,1,3},{277,1,3},{278,1,3},{279,1,3},{280,1,3},{281,1,3},{282,1,3},{283,1,3},{284,1,3},{285,1,3},{286,1,3},{287,1,3},{288,1,3},{289,1,3},{290,1,3},{291,1,3},{292,1,3},{293,1,3},{294,1,3},{295,1,3},{296,1,3},{297,1,3},{298,1,3},{299,1,3},{300,1,3},{301,1,3},{302,1,3},{303,1,3},{304,1,3},{305,1,3},{306,1,3},{307,1,3},{308,1,3},{309,1,3},{310,1,3},{311,1,3},{312,1,3},{313,1,3},{314,1,3},{315,1,3},{316,1,3},{317,1,3},{318,1,3},{319,1,3},{320,1,3},{321,1,3},{322,1,3},{323,1,3},{324,1,3},{325,1,3},{326,1,3},{327,1,3},{328,1,3},{329,1,3},{330,1,3},{331,1,3},{332,1,3},{333,1,3},{334,1,3},{335,1,3},{336,1,3},{337,1,3},{338,1,3},{339,1,3},{340,1,3},{341,1,3},{342,1,3},{343,1,3},{344,1,3},{345,1,3},{346,1,3},{347,1,3},{348,1,3},{349,1,3},{350,1,3},{351,1,3},{352,1,3},{353,1,3},{354,1,3},{355,1,3},{356,1,3},{357,1,3},{358,1,3},{359,1,3},{360,1,3},{361,1,3},{362,1,3},{363,1,3},{364,1,3},{365,1,3},{366,1,3},{367,1,3},{368,1,3},{369,1,3},{370,1,3},{371,1,3},{372,1,3},{373,1,3},{374,1,3},{375,1,3},{376,1,3},{377,1,3},{378,1,3},{379,1,3},{380,1,3},{381,1,3},{382,1,3},{383,1,3},{384,1,3},{385,1,3},{386,1,3},{387,1,3},{388,1,3},{389,1,3},{390,1,3},{391,1,3},{392,1,3},{410,4,3},{414,4,3},{417,4,3},{421,4,3},{425,4,3},{429,4,1},{432,4,1},{436,4,1}};

    test_nist_from_gsa(params, n, m, q, D_e, D_s,strategy);
}