cd framework
export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/local/bin/:/usr/include/"
g++ -O3 -march=native -Wp,-U_FORTIFY_SOURCE -fPIC -Ofast -ftree-vectorize -funroll-loops -std=c++14 -o ../libest.so utils.cpp bkz_with_jump_simulator.cpp d_svp_prediction.cpp enumbs.cpp  bssa.cpp cost.cpp attack.cpp est.cpp  -L. -pthread -lgmp -lfplll -lmpfr -Wall -Wextra -fPIC -shared -lalglib


cd ..
  
