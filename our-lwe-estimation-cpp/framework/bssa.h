#include "bkz_with_jump_simulator.h"


class BSSA{
    public:
        //implement tours rounds of BKZ-beta-J
        struct strategy{
            int beta; 
            int jump;
            int tours;
        };
        struct blocksize_strategy{
            vector<strategy> S; //pnj-BKZ strategy: 
            vector<double> l;
            double G_BKZ;
            double cumulated_proba;          
        };

        map<int,blocksize_strategy> BS; //key: beta;  value: blocksize_strategy
        double G_min, G_tmp;
        blocksize_strategy bs0;
        int beta_tmp;
        BKZJSim* sim;

        EnumBS(){
            sim = new BKZJSim();
        }


        void bssa_generate(vector<double> l0, int beta_start, int beta_goal, int jump = 1, int gap = 1, int J_gap = 1); //Generate strategy for bssa

        pair<double,int>  max_tour_for_pnjbkz_beta(vector<double> l0,double cumulated_proba0,int beta, int jump); //maximal tours for one pnj-bkz-beta-J to reduce the basis.

        int min_tour_to_each_goal_dsvp(vector<double> l0,double cumulated_proba0,int beta, int jump, double goal_dsvp); //Find the minial tours for a pnj-BKZ-beta-J to reach a basis quality of goal_dsvp.
        
        double dsvp_predict(vector<double> l, int beta, double cumulated_proba, bool progressive_sieve = true);
        double fixed_dsvp_predict(vector<double> l, int beta, double cumulated_proba);
        double progressive_dsvp_predict(vector<double> l, int beta, double cumulated_proba);

        // void d_svp_predict_for_beta(); //Input: beta, jump and l, to determine the quality of a BKZ-beta-J reduced quality.

        void print_strategy(vector<strategy> S);

};