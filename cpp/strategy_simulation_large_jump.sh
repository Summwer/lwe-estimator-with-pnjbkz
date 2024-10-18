export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
g++ -O3 -funroll-loops -o strategy_simulation_large_jump strategy_simulation_large_jump.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest
./strategy_simulation_large_jump | tee lwechal-est-result/simulation_large_jump.log
