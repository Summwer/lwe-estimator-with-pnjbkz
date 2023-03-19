#include "../framework/cost.h"



int main(){
    COST* cost = new COST();

    int d = 323;
    for(int beta = 50; beta < 100 ; beta ++){
        for(int jump = 15; jump >0; jump -- ){
            cout<<"beta = "<<beta<<", jump = "<<jump<<": ";
            cout<<cost->theo_bkz_cost(d,beta,jump).first<<endl;
        }
    }
    

    return 1;
}