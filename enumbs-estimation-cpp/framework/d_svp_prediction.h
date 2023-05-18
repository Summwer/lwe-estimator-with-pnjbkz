#include "cost.h"


tuple<double,int,double,double> dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool progressive_sieve = true, bool worst_case = false);
tuple<double,int,double,double> fixed_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool worst_case = false);
tuple<double,int,double,double> progressive_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool wost_case = false);