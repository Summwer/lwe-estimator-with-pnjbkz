export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
g++ -O3 -funroll-loops -o dilithium1-enumbs-parallel-est-gate./NIST-round3/dilithium1-enumbs-parallel-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./dilithium1-enumbs-parallel-est-gate
