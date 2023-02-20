export LD_LIBRARY_PATH="./:/usr/local/lib/"
g++ -O3 -funroll-loops -o kyber1-enumbs-parallel-est-gate ./NIST-round3/kyber1-enumbs-parallel-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./kyber1-enumbs-parallel-est-gate
