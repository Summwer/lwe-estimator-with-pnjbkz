#include "cost.h"





tuple<double,int,double,double> dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool progressive_sieve = true);
tuple<double,int,double,double> fixed_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1);
tuple<double,int,double,double> progressive_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1);