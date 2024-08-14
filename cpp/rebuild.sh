cd framework
export LD_LIBRARY_PATH="./:/usr/bin/" #"./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/local/bin/:/usr/include/:/usr/lib/"
g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -ftree-vectorize -funroll-loops -std=c++14 -o ../libest.so utils.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp enumbs.cpp  bssa.cpp cost.cpp attack.cpp est.cpp load_lwechal.cpp lwe_estimation.cpp -L. -pthread -lgmp -lmpfr  -fPIC -shared -lalglib -lcurl -lfplll  #`pkg-config --libs fplll`  #-Wall -Wextra

# g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -ftree-vectorize -funroll-loops -std=c++14 -o ../libest.so utils.cpp  lwe_estimation.cpp   -L. -pthread -lgmp -lmpfr  -fPIC -shared -lalglib -lcurl -lfplll  #`pkg-config --libs fplll`  #-Wall -Wextra


#bkz_with_jump_simulator.cpp d_svp_prediction.cpp cost.cpp load_lwechal.cpp

#load_lwechal.cpp lwe_estimation.cpp 

#lwe_estimation.cpp attack.cpp est.cpp

cd ..
  
