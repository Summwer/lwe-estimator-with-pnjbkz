// #include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include "framework/est.h"

using namespace std;
// using namespace boost;


int main(){
    // Params* params = new Params; //J, gap, J_gap, cost_model, verbose,
    // params->max_dim = 500; // maximal selected blocksize value
    // // Kyber-I(Kyber-512) round-3 parameters
    // double cum_pr = 0., rem_pr = 1., pr;
   
    // int dim = 1004;
    // int beta = 200, jump = 1, d = dim;
    // FP_NR<FT> dvol = 3882.6780896;
    // vector<double> l = gen_simulated_gso(dim, dvol);

    // pair<double,double> GB, cum_GB;


    // BKZJSim* sim = new BKZJSim();
    // COST* cost = new COST();

    // sim -> simulate(l,l,beta,jump,1);

    // boost::math::chi_squared chisquare(beta);
    // pr = boost::math::cdf(chisquare,pow(2,2.*l[d-beta]));

    // GB = cost->bkz_cost(d,beta,jump,params->cost_model);

    // cout<<GB.first<<endl;
    // cout<<pow(2,2.*l[d-beta])<<endl;
    // cout<<pr<<endl;

    // cum_GB.first = log2(pow(2,cum_GB.first)+(pow(2,GB.first)*rem_pr*pr));

    // cout<<pr<<endl;

    // cout<<cum_GB.first<<endl;

    double t1 = 11232131.145343423;

    double t2 = round(t1*0.1)/0.1;

    cout<< t1 << endl;

    cout<< t2 <<endl;

    int d = 194;

    cout<<int(0.9*d)<<endl;

}