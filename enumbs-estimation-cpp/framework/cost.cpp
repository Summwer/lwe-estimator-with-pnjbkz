#include "cost.h"



// Return log2 of the number of gates for FindAllPairs according to AGPS20
FP_NR<FT> COST::agps20_gates(int beta_prime){
    FP_NR<FT> k = ((double)beta_prime) / 8.;
    FP_NR<FT> d1,d2,x;
    if(k != round(k.get_d())){
        x = k - floor(k);
        d1 = agps20_gates(8*(int)floor(k).get_d());
        d2 = agps20_gates(8*((int)floor(k).get_d() + 1));
        return x * d2 + (1. - x) * d1;
    }
    return COST::agps20_gate_data[beta_prime];
}



//cost of bkz with progressive sieve
pair<double,double> COST::theo_bkz_cost(int n, int beta,int J){
    /*Return cost of bkz-beta in theoretical Gate: T/G -- time cost, B -- memory cost*/
    if(beta <=10)
        return make_pair(0.,0.);
    int beta_prime = floor(beta - dims4free(beta));
    if(beta_prime < 64 or beta < beta_prime)
        return make_pair(0.,0.);
    else if(beta_prime > 1024)
        return make_pair(MAX_NUM,MAX_NUM);
    else{
        //double gates = log2((((double)n-beta)/J)*COST::C*COST::C) + agps20_gates(beta_prime).get_d();
        double gates = log2((((double)n-beta)/J)*COST::C) + agps20_gates(beta_prime).get_d();
        double bits = log2(8*beta_prime); //+ agps20_vectors(beta_prime)
        return make_pair(gates, bits);
    }
}

pair<double,double> COST::theo_pump_cost(int beta){
    /*Return cost of pump-beta in theoretical Gate: T/G -- time cost, B -- memory cost*/
    if(beta <=10)
        return make_pair(0.,0.);
    int beta_prime = floor(beta - dims4free(beta));
    if(beta_prime < 64 or beta < beta_prime)
        return make_pair(0.,0.);
    else if(beta_prime > 1024)
        return make_pair(MAX_NUM,MAX_NUM);
    else{
        //double gates = log2(COST::C*COST::C) + agps20_gates(beta_prime).get_d();
        double gates = log2(COST::C) + agps20_gates(beta_prime).get_d();
        // double bits = log2(8.*beta_prime); // + agps20_vectors(beta_prime)
        double bits = 0.;
        return make_pair(gates, bits);
        
    }
}



//threads = 32, gpus = 2,  pnj-bkz cost
pair<double,double> COST::get_k1_k2_pnj(int beta,bool sieve){
    double k1,k2;
    if(beta >=0 and beta <10){
        k1 = 0;
        k2 = 0;
    }
    else if(beta>=10 and beta<=42 and sieve == false){
        k1 = 0.03;
        k2 = 5.188;
    }
    else if(beta<=60 and sieve == false){
        k1 = 0.19;
        k2 = -1.741;
    }
    else if(beta <= 97){
        k1 = 0.056;
        k2 = 7.85;
    }
    else if(beta <= 118){
        k1 = 0.215;
        k2 = -7.61;
    }
    else if(beta <= 128){
        k1 = 0.314;
        k2 = - 19.24;
    }
    else{
        k1 = 0.368;
        k2 = -26.15;
    }
    return make_pair(k1,k2);
}




//threads = 32, gpus = 2, pump cost paramter
pair<double,double> COST::get_k1_k2_pump(int beta){
    double k1,k2;
    if(beta >=0 and beta <10){
        k1 = 0 ;
        k2 = 0 ;
    }
    else if(beta>=10 and beta<=60){
        k1 = 0.035657;
        k2 = -2.317327;
    }
    else if(beta <= 96){
        k1 = 0.078794;
        k2 = -0.039742;
    }
    else if(beta <= 116){
        k1 = 0.231927;
        k2 = -14.713430;
    }
    else if(beta <= 128){
        k1 = 0.314;
        k2 = -24.21;
    }
    else{
        k1 = 0.368;
        k2 = -31.12;
    }
    
    return make_pair(k1,k2);
}



//get pump cost in threads = 20
double COST::practical_pump_cost(int beta){
    //make sure not use the enum cost 
    pair<double,double> k = get_k1_k2_pump(beta); // threads = 20
    double k1 = k.first, k2 = k.second;
    // k = (1/71.)*((1.33)**(beta/10.));

    return k1*((double)beta)+k2; //n_expected = beta -f , beta = d-llb
}
    

//get pnj-BKZ time test in threads = 20
double COST::practical_bkz_cost(int d,int beta,int f,int jump){
    bool sieve;
    if(beta <= 60)
        sieve = false;
    else
        sieve = true;  
    pair<double,double> k = get_k1_k2_pnj(beta,sieve); // threads = 20
    double k1 = k.first, k2 = k.second;
    double c3= 0.018, c4 = -2.24;

    return (k1*(beta-f)+k2) + log2((c3*d+c4)/jump);
}


    
pair<double,double> COST::pump_cost(int beta,int cost_model){
    if(cost_model == 1)
        return theo_pump_cost(beta);
    else if(cost_model == 2)
        return make_pair(practical_pump_cost(beta),theo_pump_cost(beta).second);
    return make_pair(MAX_NUM,MAX_NUM);
}

pair<double,double> COST::bkz_cost(int d, int beta,int J,int cost_model){
    if(cost_model == 1)
        return make_pair(theo_bkz_cost(d, beta, J).first, 0.);
    else if(cost_model == 2){
        int f = dims4free(beta);
        return make_pair(practical_bkz_cost(d,beta,f,J), 0.);//theo_bkz_cost(d, beta,J).second);
    }
    return make_pair(MAX_NUM,MAX_NUM);
}    




