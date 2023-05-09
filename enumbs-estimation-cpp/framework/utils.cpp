#include "utils.h"

// return all keys in map
vector<int> KeySet(map<int, rational<int>> mp)
{
    vector<int> keys;
    for(map<int, rational<int>>::iterator it = mp.begin(); it != mp.end(); ++it){
        keys.push_back(it->first);
    }
    return keys;
}

// return all values in map
vector<double> ValueSet(map<int, rational<int>> mp)
{
    vector<double> values;

    for(map<int, rational<int>>::iterator it = mp.begin(); it != mp.end(); ++it){
        values.push_back(rational_cast<double>(it->second));
    }
    return values;
}


// Generate samples list such that number of each sample obeys a specific prob * sample_num. 
vector<int> expand_samples(vector<int> samples,vector<double> probs,int sample_num){
    vector<int> expand_samples_list;
    int index = 0;
    int upper_index = 0;

    for(int i=0; i< int(samples.size()); i++){
        if(i < int(samples.size())-1)
            upper_index += sample_num*probs[i];
        else
            upper_index = sample_num;

        expand_samples_list.insert(expand_samples_list.begin()+index, upper_index-index, samples[i]);
        index = upper_index;
    }
    return expand_samples_list;
}


void print_map(map<int, rational<int>>  mp){
    cout<<"{\t";
    for (map<int, rational<int>> ::iterator it = mp.begin();
		it != mp.end(); it++) {
            cout<<it-> first<<","<<it->second<<"\t";
	}
    cout<<"}"<<endl;
}

void print_map(map<int, FP_NR<FT>>  mp){
    cout<<"{\t";
    for (map<int, FP_NR<FT>> ::iterator it = mp.begin();
		it != mp.end(); it++) {
            cout<<it-> first<<","<<it->second<<"\t";
	}
    cout<<"}"<<endl;
}

void print_vector(vector<double> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}

void print_vector(vector<int> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}

void print_vector(vector<Z_NR<ZT>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}

void print_vector(vector<FP_NR<FT>> v, int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}

void print_vector(vector<pair<int,int>> v, int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << (v[i].first) << ", " << v[i].second << ") ";
    }
    cout<<"]"<<endl;
}


void print_vector(vector<pair<double,double>> v, int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << (v[i].first) << ", " << v[i].second << ") ";
    }
    cout<<"]"<<endl;
}


void print_vector(vector<tuple<double,double,double>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << get<0>(v[i]) << ", " << get<1>(v[i]) << ", " << get<2>(v[i])  << ") ";
    }
    cout<<"]"<<endl;
}


void print_vector(vector<tuple<int,int,int>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << "(" << get<0>(v[i]) << ", " << get<1>(v[i]) << ", " << get<2>(v[i])  << ") ";
    }
    cout<<"]"<<endl;
}

/*
void print_matrix(ZZ_mat<ZT> matrix){
    cout<< "[";
    for(auto it = matrix.begin(); it != matrix.end(); ++it) {
        cout << *it << " ";
    }
    cout<<"]"<<endl;
}
*/


void print_matrix(vector<vector<Z_NR<ZT>>> matrix){
    int row = matrix.size(), col = matrix[0].size();

    printf("[");
    for(int i = 0; i < row; i++){
        print_vector(matrix[i],0,col);
    }
    printf("]\n");
}

void printf_input(int d, FP_NR<FT> dvol)
{
    printf("\033[0m\033[1;31m dim = %3d, dvol = %3.2f \033[0m\n", d, dvol.get_d());
}


void printf_red(const char *s)
{
    printf("\033[0m\033[1;31m%s\033[0m\n", s);
}

pair<rational<int>,rational<int>> average_variance(std::map<int,rational<int>> D){
    rational<int> mu,s;

    int v;
    rational<int> p;
    
    for(auto it : D){
        v = it.first;
        p = it.second;
        mu += v * p;
        s += v * v * p;
    }

    s -= mu * mu;

    return make_pair(mu,s);
}


int draw_from_distribution(std::map<int,rational<int>> D, int sample_num){
    /*draw an element from the distribution D
    ,D, distribution in a dictionnary form
    */

    // double rational_cast<double>
    vector<int> samples = KeySet(D);
    vector<double> probs= ValueSet(D);
    vector<int> expand_samples_list = expand_samples(samples,probs,sample_num);

    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> distrib(0,sample_num-1);
    

    return expand_samples_list[distrib(gen)];

    //nc::random::choice(a, 3); //Corresponding to numpy.random.choice function in numpy.

    //X = np_random_choice([key for key in D.keys()],
    //                     1, replace=True,
    //                     p=[float(prob) for prob in D.values()])
    //return X[0]
}


void gen_samples(ZZ_mat<ZT> &matrix, int m, int n, int q){
    /*
    m, number of samples
    n, dimension for each sample
    */
    
    
    Z_NR<ZT> q2;
    q2=q;
    
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j].randm(q2);
}

void build_LWE_lattice(ZZ_mat<ZT> &matrix, ZZ_mat<ZT> A, int q){
    /*
    Builds a (m+n)*(m+n) matrix of the form
    q*I, 0
    -A^T, I
    It corresponds to the LWE matrix
    ,A, a m*n matrix
    ,q, integer
    */

    int m = A.get_rows(), n = A.get_cols(), d = m+n;

    //draw matrix B and define the lattice, B = [[qI,0],[-A^T,I]], A is lattice samples. Each col is one group of sample.
    matrix.gen_zero(d,d);
    for (int i = 0; i < m; i++)
        matrix[i][i] = q;

    for (int i = m; i < d; i++)
        matrix[i][i] = 1;

    for (int i = m; i < d; i++)
        for (int j = 0; j < m; j++)
            matrix[i][j].mul_si(A[j][i-m],-1);
   
}



void kannan_embedding(ZZ_mat<ZT> &matrix, vector<Z_NR<ZT>> target, int factor){
    /*
    Builds a (m+n)*(m+n) matrix of the form
    q*I, 0,  0
    -A^T, I, 0
     b ,  0, 1
    It corresponds to the LWE matrix
    ,A, a m*n matrix
    ,q, integer
    */

    int d = matrix.get_rows()+1;
    matrix.resize(d,d);
    
    for (int i = 0; i < d-1; i++){
        matrix[d-1][i] = target[i];
    }
    
    matrix[d-1][d-1] = factor;
    
}


FP_NR<FT> compute_delta( int beta){
    /*
    Computes delta from the block size k. Interpolation from the following
    data table,
    Source , https,//bitbucket.org/malb/lwe-estimator/
    src/9302d4204b4f4f8ceec521231c4ca62027596337/estima
    tor.py?at=master&fileviewer=file-view-default
    ,k, integer
    estimator.py table,
    */

    map<int,FP_NR<FT>> small= {{0,1e20}, {1, 1e20}, {2, 1.021900}, {3, 1.020807}, {4, 1.019713}, {5, 1.018620}, {6, 1.018128}, {7, 1.017636}, {8, 1.017144}, {9, 1.016652}, {10, 1.016160}, {11, 1.015898}, {12, 1.015636}, {13, 1.015374}, {14, 1.015112}, {15, 1.014850}, {16, 1.014720}, {17, 1.014590}, {18, 1.014460}, {19, 1.014330}, {20, 1.014200}, {21, 1.014044}, {22, 1.013888}, {23, 1.013732}, {24, 1.013576}, {25, 1.013420}, {26, 1.013383}, {27, 1.013347}, {28, 1.013310}, {29, 1.013253}, {30, 1.013197}, {31, 1.013140}, {32, 1.013084}, {33, 1.013027}, {34, 1.012970}, {35, 1.012914}, {36, 1.012857}, {37, 1.012801}, {38, 1.012744}, {39, 1.012687}, {40, 1.012631}, {41, 1.012574}, {42, 1.012518}, {43, 1.012461}, {44, 1.012404}, {45, 1.012348}, {46, 1.012291}, {47, 1.012235}, {48, 1.012178}, {49, 1.012121}, {50, 1.012065}};

    if(beta <=50)
        return small[beta];
    else{
        return pow((double)beta/(2.*M_PI*exp(1)) * pow((M_PI*(double)beta),(1./(double)beta)), (1./(2*((double)beta-1.))));
    }
}



//calculate slope in G6K
double get_current_slope(vector<double> l, int start, int end){
    /*
    Return current slope of log-gs-lengths

    :param log_rr: vector of log(squared) Gram-Schmidt norms
    :param start_row: start row
    :param stopr_row: stop row

    */
    
    int n = end - start;
    double i_mean = (n - 1) * 0.5 + start;
    double x_mean =  accumulate(l.begin(), l.end(), 0.)/n;
    double v1 = 0.0, v2 = 0.0;
    for(int i = start; i < end; i++){
        v1 += (i - i_mean) * (l[i] - x_mean);
        v2 += (i - i_mean) * (i - i_mean);
     }
    return v1 / v2 ; //round(v1 / v2 * 1e6)/1e6; 
}



FP_NR<FT> bkzgsa_gso_len( FP_NR<FT> logvol, int i, int d, FP_NR<FT> delta, int beta){
    FP_NR<FT> gso_len;
    if(delta >= 1000.)
        delta = compute_delta(beta);
    gso_len.pow_si(delta,(d-1-2*i));
    FP_NR<FT> tmp;
    tmp.exponential(logvol/(double(d)));
    gso_len.mul(gso_len,tmp);

    return gso_len;
}



// double gaussian_heuristic(vector<FP_NR<FT>> l, int index_start){
//     int n = l.size();

//     double const log_ball_square_vol = (n * log(M_PI) - 2.0 * lgamma(n / 2.0 + 1))/(log(2.));
//     FP_NR<FT> log_lattice_square_vol = 0.;
//     FP_NR<FT> gh;
//     for (int i = index_start; i < n; ++i)
//     {
//         log_lattice_square_vol += l[i];
//     }

//     //gh.exponential(1);
//     double gh = exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));
//     // cout<<log_lattice_square_vol.get_d()<<","<<log_ball_square_vol<<endl;
//     return gh.get_d();

// }

//Input: log2(gs-lengths)
double gaussian_heuristic_log2(vector<FP_NR<FT>> l, int index_start){
    int d = l.size();
    int n = d - index_start;

    double const log_ball_square_vol = (n * log(M_PI) - 2.0 * lgamma(n / 2.0 + 1));
    double log_lattice_square_vol = 0., gh;
    for (int i = index_start; i < d; ++i)
    {
        log_lattice_square_vol += (2*l[i].get_d()*log(2.));
    }
    // cout<<log_lattice_square_vol<<","<<log_ball_square_vol  <<endl;

    //FP_NR<FT> gh;
    //gh.exponential((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));
    gh = exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));

    return gh;


}

//Input: log2(gs-lengths)
double gaussian_heuristic_log2(vector<double> l, int index_start){
    int d = l.size();
    int n = d - index_start;

    double const log_ball_square_vol = (n * log(M_PI) - 2.0 * lgamma(n / 2.0 + 1));
    double log_lattice_square_vol = 0., gh;
    for (int i = index_start; i < d; ++i)
    {
        log_lattice_square_vol += (2*l[i]*log(2.));
    }
    gh = exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * n));

    return gh;


}


//Generate input simulated-gs-lengths
vector<double> gen_simulated_gso(int d, FP_NR<FT> logvol){
    FP_NR<FT> delta, gso_len;
    delta = compute_delta(2);
    
    vector<FP_NR<FT>> l;
    l.resize(d);
    for(int i = 0; i<d; i++){
        gso_len = bkzgsa_gso_len(logvol, i, d, delta);
        l[i].log(gso_len);
        l[i].div(l[i],log(2));
    }
    // print_vector(l);

    vector<double> l_;
    l_.resize(d);
    for(int i = 0; i<d; i++){
        l_[i] = l[i].get_d();
    }

    
    return l_;
    
}



//=theo_dim4free_fun2
int dims4free(int beta){
    return max(0,int(ceil((double)beta * log(4./3.) / log((double)beta/(2.*M_PI*exp(1.))))));
}

int default_dim4free_fun(int beta){
    if(beta < 40)
        return 0;
    int f = int(11.5 + 0.075*beta);
    return min(int(((double)beta - 40)/2.), f);
}

int theo_dim4free_fun1(int beta){
    if(beta < 40)
        return 0;
    int f = max(0,int(ceil((double)beta * log(4./3.) / log((double)beta/(2.*M_PI)))));
    return min(int(((double)beta - 40)/2.), f);
    
}

int theo_dim4free_fun2(int beta){
    if(beta < 40)
        return 0;
    return min(int(((double)beta - 40)/2.), dims4free(beta));
}


int get_beta_from_sieve_dim(int sieve_dim, int d, int choose_dims4f_fun){
    int f;
    for(int beta = sieve_dim; beta < d; beta++){
        if(choose_dims4f_fun == 1 )
            f = theo_dim4free_fun1(beta);
        else if(choose_dims4f_fun == 2)
            f = theo_dim4free_fun2(beta);
        if(beta - f >= sieve_dim)
            return beta;
    }
    return d;
}