#include <iostream>
#include <map>
#include <vector>
#include <boost/rational.hpp>
#include<bits/stdc++.h>
#include <fplll.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <chrono>



using namespace std;
using namespace boost;
using namespace fplll;
using namespace std::chrono; 

#ifndef ZT
#define ZT mpz_t
#endif

#ifndef FT
#define FT mpfr_t
#endif


#ifndef MAX_NUM
#define MAX_NUM 1000
#endif

#ifndef MAX_DIM
#define MAX_DIM 1000
#endif



#ifndef PREC
#define PREC 1e3
#endif

void print_map(map<int, rational<int>> mp);
void print_map(map<int, FP_NR<FT>>  mp);
void printf_red(const char *s);
void printf_input(int d, FP_NR<FT> dvol);
void print_vector(vector<double> v,int index_start=0, int index_end=-1);
void print_vector(vector<int> v,int index_start=0, int index_end=-1);
void print_vector(vector<Z_NR<ZT>> v,int index_start=0, int index_end=-1);
void print_vector(vector<FP_NR<FT>> v,int index_start=0, int index_end=-1);
void print_vector(vector<pair<int,int>> v,int index_start=0, int index_end=-1);

//void print_matrix(ZZ_mat<ZT> matrix);

pair<rational<int>,rational<int>> average_variance(std::map<int,rational<int>> D);
vector<int> KeySet(map<int, rational<int>> mp);
vector<double> ValueSet(map<int, rational<int>> mp);
int draw_from_distribution(std::map<int,rational<int>> D, int sample_num=1000);
vector<int> expand_samples(vector<int> samples,vector<double> prob,int sample_num);


void gen_samples(ZZ_mat<ZT> &matrix, int m, int n, int q);
void build_LWE_lattice(ZZ_mat<ZT> &matrix, ZZ_mat<ZT> A, int q);
void kannan_embedding(ZZ_mat<ZT> &matrix, vector<Z_NR<ZT>> target, int factor=1);

FP_NR<FT> compute_delta(int beta);
FP_NR<FT> bkzgsa_gso_len( FP_NR<FT> logvol, int i, int d, FP_NR<FT> delta, int beta=-1);
vector<double> gen_simulated_gso(int d, FP_NR<FT> logvol);

void simulator_test(int d, FP_NR<FT> logvol); //test simulator 

double gaussian_heuristic(vector<FP_NR<FT>> l, int index_start);
double gaussian_heuristic_log2(vector<FP_NR<FT>> l, int index_start);
double gaussian_heuristic_log2(vector<double> l, int index_start);

int dims4free(int beta); //leaky-lwe-estimator
int default_dim4free_fun(int beta);
int theo_dim4free_fun1(int beta);
int theo_dim4free_fun2(int beta);
int get_beta_from_sieve_dim(int sieve_dim, int d, int choose_dims4f_fun);

struct Params{
    int J = 15; //J -- maximal jump value;
    int gap = 1; //gap -- gap of each beta;
    int J_gap = 1; //J_gap -- gap of each jump value;
    //cost_model: 1: gate model; 2: sec model with threads=32, gpus = 2 
    int cost_model = 1; 
    bool verbose = true; //print logging or not
    //progressieve_sieve: True: progressieve sieve; False: normal sieve
    bool progressive_sieve =  true; 
    int threads = 1;
    int max_dim = MAX_DIM; //set the maximal blocksize to find the optimal strategy
    int max_loop = 10; //set the maximal loop for one blocksize to find the optimal strategy

    int method = 1; //1: enumbs estimation; 2: bssa estimation

    bool debug = false; //print debug logging or not.
    bool verification =false; //verify the correctness of strategy


    //enumbs params
    double enumbs_prec =  PREC; //1e-5; //set the precision of enumbs
    int enumbs_add_strategy = 3; //1: remain all equal strategies; 2: (Large block strategy) delete all equal strategies; 3: (Small block strategy) remain the small blocksize strategy first.


};


struct LWEchal{
    int n;
    double alpha;
    int dim;
    FP_NR<FT> dvol;
};