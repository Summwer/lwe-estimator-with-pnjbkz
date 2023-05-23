#include "../framework/bkz_with_jump_simulator.h"


void strategy_verification(vector<double> l,vector<tuple<int,int,int>> strategy, int cost_model,int progressive_sieve, bool worst_case){
    int dim =  l.size();
    double cum_pr = 0., rem_pr = 1., proba, G1cum=0., B1cum = 0.;
    BKZJSim* sim = new BKZJSim(params);
    COST* cost = new COST();

    for(int i = 0; i< int(strategy.size()); i++){
        tuple<int,int,int> bs = strategy[i];
        int beta = get<0>(bs), jump = get<1>(bs), N = get<2>(bs);
        for(int tour = 0; tour < N; tour++){
        
            sim -> simulate(l,l,beta,jump,1);

            boost::math::chi_squared chisquare(beta);
            proba = boost::math::cdf(chisquare,pow(2,2.*l[d-beta]));

            pair<double,double> GB = cost -> bkz_cost(d,beta,jump,cost_model);
            G1cum = log2(pow(2,G1cum) + (pow(2,GB.first) * rem_pr * proba));
            B1cum = max(B1cum,GB.second);

            cum_pr += rem_pr * proba;
            rem_pr *= 1. - proba;
        }
    }

    tuple<double,int,double,double> dsvp_t =  dsvp_predict(l, cum_pr, cost,cost_model, progressive_sieve, worst_case);
    double G2 = get<2>(dsvp_t);
    printf("%3.2f, %3.2f, %3.2f\n", G1cum,G2,get<0>(dsvp_t));
    printf("%3.2f\n", log2(pow(2,G1cum)+pow(2,G2)));

}

int main(){
    // int dim =  2049;
    // FP_NR<FT> dvol = 15614.219317244602;

    int dim =  287;
    FP_NR<FT> dvol = 756.43 ;

    vector<double> l = gen_simulated_gso(d, dvol);
    int cost_model = 2, progressive_sieve = true;
    bool worst_case = false;

    vector<tuple<int,int,int>> strategy = {{ 51, 15, 11},{ 51, 13, 23}};
    strategy_verification(l,strategy,cost_model, progressive_sieve,worst_case);
}