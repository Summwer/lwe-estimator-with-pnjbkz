#include "bkz_with_jump_simulator.h"


void BKZJSim::init(){
    cd.resize(CD_NUM);
    FP_NR<FT> tmp_sum;
    int beta;
    for(int i = 0; i < CD_NUM; i++){
        beta = i + 1;
        if(i < 45){
            tmp_sum = 0.;
            for(int j = 45-beta; j < 45; j++){
                tmp_sum += rk[j];
            }
            cd[i] = (rk[45-beta] -  tmp_sum)/beta;

        }
        else{
            cd[i] = (lgamma( beta / 2.0 + 1)*(1./beta) - log(sqrt(M_PI)))/log(2.);
        }
    }
            
}


void BKZJSim::simulate(vector<double> &l_,vector<double> l,int beta, int jump, int N){
    if(beta >= 45)
        sim_above_45(l_,l,beta, jump, N);
    else
        sim_below_45(l_,l,beta, N);

}




void BKZJSim::sim_below_45(vector<FP_NR<FT>> &l_,vector<FP_NR<FT>> l, int beta, int N){
    /*
    Simulate gs-lengths after bkz-beta-jump with N rounds, beta < 45.
    @param beta: blocksize
    @param jump: jump value
    @param N: round number of bkz-beta-jump
    */
    
    if(N==0 or beta ==0)
        return;
    
    int d = l.size();
    bool flag;
    int beta_,f=0,ff=0;
    FP_NR<FT> sumf=0., sumk=0., logV, l_k;
    l_.resize(d);

    for(int tours = 0; tours < N; tours ++){
        flag = true;
        sumf=0.;
        sumk=0.;
        ff = 0;
        for(int k = 0; k < d; k++){
            beta_ = min(beta, d - k);
            f = min(k + beta, d);
            for(int i = ff; i<f; i++)
                sumf += l[i];
            if(k!=0)
                sumk += l_[k-1];
            logV = sumf - sumk;
            l_k = logV / beta_ + cd[beta_-1];
            ff = f;
            if(flag){
                if(l_k < l[k]){
                    l_[k] = l_k;
                    flag = false;
                }
                else
                    l_[k] = l[k];  
            }
            else{
                l_[k] = l_k;
            }           
        }
        
        l = l_;
    }            
}


void BKZJSim::sim_below_45(vector<double> &l_,vector<double> l, int beta, int N){
    /*
    Simulate gs-lengths after bkz-beta-jump with N rounds, beta < 45.
    @param beta: blocksize
    @param jump: jump value
    @param N: round number of bkz-beta-jump
    */
    
    if(N==0 or beta ==0)
        return;
    
    int d = l.size();
    bool flag;
    int beta_,f=0,ff=0;
    FP_NR<FT> sumf=0., sumk=0., logV, l_k;
    l_.resize(d);

    for(int tours = 0; tours < N; tours ++){
        flag = true;
        sumf=0.;
        sumk=0.;
        ff = 0;
        for(int k = 0; k < d; k++){
            beta_ = min(beta, d - k);
            f = min(k + beta, d);
            for(int i = ff; i<f; i++)
                sumf += l[i];
            if(k!=0)
                sumk += l_[k-1];
            logV = sumf - sumk;
            l_k = logV / beta_ + cd[beta_-1];
            ff = f;
            if(flag){
                if(l_k < l[k]){
                    l_[k] = l_k.get_d();
                    flag = false;
                }
                else
                    l_[k] = l[k];  
            }
            else{
                l_[k] = l_k.get_d();
            }           
        }
        
        l = l_;
    }            
}

void BKZJSim::sim_above_45(vector<FP_NR<FT>> &l_,vector<FP_NR<FT>> l,int beta, int jump, int N){
    /*
    Simulate gs-lengths after bkz-beta-jump with N rounds, beta >= 45.
    @param beta: blocksize
    @param jump: jump value
    @param N: round number of bkz-beta-jump
    */
    
    if(N==0 or beta ==0)
        return;
    
    int d = l.size();
    bool flag;
    int beta_,k_ = min(45,beta),f=0,ff=0;
    FP_NR<FT> sumf=0., sumk=0., logV, l_k;
    l_.resize(d);

    for(int tours = 0; tours < N; tours ++){
        flag = true;
        sumf=0.;
        sumk=0.;
        ff = 0;
        for(int k = 0; k < d - beta; k++){
            beta_ = min(beta, d - k);
            f = min(k - (k % jump) + beta, d);
            for(int i = ff; i<f; i++)
                sumf += l[i];
            if(k!=0)
                sumk += l_[k-1];
            logV = sumf - sumk;
            l_k = logV / (beta_-(k % jump)) + cd[(beta_-(k % jump))-1];
            ff = f;
            if(flag){
                if(l_k < l[k]){
                    l_[k] = l_k;
                    flag = false;
                }
                else
                    l_[k] = l[k];  
            }
            else{
                l_[k] = l_k;
            }           
        }
        
        for(int i = ff; i < d; i++)
            sumf += l[i];
        for(int k = d-beta; k < d-45; k++){
            beta_ = d - k;
            sumk += l_[k-1];

            logV = sumf - sumk;
            l_k = logV / beta_ + cd[beta_-1];
            if(flag){
                if(l_k < l[k]){
                    l_[k] = l_k;
                    flag = false;
                }
                else
                    l_[k] = l[k];
            }
            else
                l_[k] = l_k;
        }
        //last 45 norms
        logV -= l_[d-46];
        for(int k = d - k_; k < d; k++)
            l_[k] = logV/k_ + rk[k + 45 - d];
    
        l = l_;
    }
}




void BKZJSim::sim_above_45(vector<double> &l_,vector<double> l,int beta, int jump, int N){
    /*
    Simulate gs-lengths after bkz-beta-jump with N rounds, beta >=45
    @param beta: blocksize
    @param jump: jump value
    @param N: round number of bkz-beta-jump
    */
    
    if(N==0 or beta ==0)
        return;
    
    int d = l.size();
    bool flag;
    int beta_,k_ = min(45,beta),f=0,ff=0;
    FP_NR<FT> sumf=0., sumk=0., logV, l_k;
    l_.resize(d);

    for(int tours = 0; tours < N; tours ++){
        flag = true;
        sumf=0.;
        sumk=0.;
        ff = 0;
        for(int k = 0; k < d - beta; k++){
            beta_ = min(beta, d - k);
            f = min(k - (k % jump) + beta, d);
            for(int i = ff; i<f; i++)
                sumf += l[i];
            if(k!=0)
                sumk += l_[k-1];
            logV = sumf - sumk;
            l_k = logV / (beta_-(k % jump)) + cd[(beta_-(k % jump))-1];
            ff = f;
            if(flag){
                if(l_k < l[k]){
                    l_[k] = l_k.get_d();
                    flag = false;
                }
                else
                    l_[k] = l[k];  
            }
            else{
                l_[k] = l_k.get_d();
            }           
        }
        
        for(int i = ff; i < d; i++)
            sumf += l[i];
        for(int k = d-beta; k < d-45; k++){
            beta_ = d - k;
            sumk += l_[k-1];

            logV = sumf - sumk;
            l_k = logV / beta_ + cd[beta_-1];
            if(flag){
                if(l_k < l[k]){
                    l_[k] = l_k.get_d();
                    flag = false;
                }
                else
                    l_[k] = l[k];
            }
            else
                l_[k] = l_k.get_d();
        }
        //last 45 norms
        logV -= l_[d-46];
        for(int k = d - k_; k < d; k++)
            l_[k] = (logV/k_ + rk[k + 45 - d]).get_d();
    
        l = l_;
    }
}


