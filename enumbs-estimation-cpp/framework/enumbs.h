#include "bssa.h"
#include <fplll/threadpool.h>




class EnumBS{
    public:
        //implement tours rounds of BKZ-beta-J
        struct strategy{
            int beta; 
            int jump;
            int tours;
        };
        struct blocksize_strategy{
            tuple<double,int,double,double> dsvp_t;
            vector<strategy> S; //pnj-BKZ strategy: 
            vector<double> l;
            pair<double,double> cum_GB_BKZ;
            pair<double,double> cum_avg_GB_BKZ;
            pair<double,double> GB;
            double cum_pr; //cumulated probability 
            double slope;
            double min_G;
        };

        vector<blocksize_strategy> BS; //key: beta;  value: blocksize_strategy
        vector<vector<blocksize_strategy>> tmpBS; //key: beta;  value: blocksize_strategy
        
        BKZJSim* sim;
        COST* cost;
        vector<double> l0;
        bool verification = false;
        double min_G_prec = 0.1;
        thread_pool::thread_pool threadpool;
        Params* params;
        

        EnumBS(Params* params){
            sim = new BKZJSim();
            cost = new COST(params);
            this->params = params;
        }

        void set_threads(int nr);

        void print_strategy(vector<strategy> S);
        void print_BS(vector<blocksize_strategy> BS);
        void print_bs(blocksize_strategy bs);

        int find_pos_for_dsvp(int cdsvp);
        int find_pos_for_dsvp(double dsvp);
        int binary_search_for_cdsvp(double cdsvp);
        int binary_search_for_dsvp(double dsvp);
        int binary_search_for_G2(double G2);
        int binary_search_for_G(double G);
        int binary_search_for_slope(double slope);
        int binary_search_for_G2_slope_cum_pr(blocksize_strategy bs);

        pair<int,int> get_max_strategy(vector<EnumBS::strategy> S);
        bool compare_max_strategy(vector<EnumBS::strategy> S0, vector<EnumBS::strategy> S1);

        // int max_beta_in_S(vector<EnumBS::strategy> S);
        vector<double> extract_cdsvp();
        vector<double> extract_dsvp();
        vector<double> extract_G2();
        vector<pair<double,double>> extract_G2G1();
        vector<double> extract_slope();
        vector<double> extract_cum_pr();
        vector<tuple<double,double,double>> extract_G2_slope_cum_pr();
        bool no_repeated_value_verification(vector<int> nums);
        bool no_repeated_value_verification(vector<double> nums);

        void BS_add(EnumBS::blocksize_strategy bs, int k); 
        void BS_add_op(EnumBS::blocksize_strategy bs, int k);
        void BS_add_G2(EnumBS::blocksize_strategy bs, int k);
        // void BS_add_slope(EnumBS::blocksize_strategy bs, int k);
        // pair<int,bool> BS_add_G2_backup(EnumBS::blocksize_strategy bs, int k);
        
        // void BS_add_cdsvp(EnumBS::blocksize_strategy bs, int k);
        // bool BS_add_determine(EnumBS::blocksize_strategy bs, int k);


        bool pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB_BKZ,  pair<double,double> &cum_avg_GB_BKZ,  pair<double,double> &GB, double &cum_pr, int beta, int jump, tuple<double,int,double,double> &dsvp_t_, double &slope);


        void max_tour_for_pnjbkz_beta(int k, int beta,int jump);
        void max_tour_for_pnjbkz_beta_G2(int k, int beta,int jump);
        pair<int,bool> max_tour_for_pnjbkz_beta_G2_backup(int k, int beta,int jump);
        void enumbs_est(vector<double> l0); //Generate strategy for enumbs

        void max_tour_for_pnjbkz_beta_in_parallel(int  beta_j_t_id_begin, vector<pair<int,int>> beta_j_tid,int k);
        
        void enumbs_est_in_parallel(vector<double> l0); //Generate strategy for enumbs in parallel
        

        pair<double,double> strategy_verification( vector<double> l,vector<strategy> strategy);

      

};