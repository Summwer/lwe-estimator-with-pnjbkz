#include "cost.h"

tuple<int,int,double,double> dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool progressive_sieve = true, bool worst_case = false);
tuple<int,int,double,double> fixed_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool worst_case = false);
tuple<int,int,double,double> progressive_dsvp_predict1(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool worst_case = false);
tuple<int,int,double,double> progressive_dsvp_predict2(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1);