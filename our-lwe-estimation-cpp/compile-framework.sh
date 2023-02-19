cd framework
export LD_LIBRARY_PATH="enumbs-estimation-cpp"
g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -ftree-vectorize -funroll-loops -std=c++14 -o ../libest.so utils.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp enumbs.cpp cost.cpp call_enumbs.cpp est.cpp  -pthread -lgmp -lfplll -lmpfr -Wall -Wextra -fPIC -shared 
cd ..
  
