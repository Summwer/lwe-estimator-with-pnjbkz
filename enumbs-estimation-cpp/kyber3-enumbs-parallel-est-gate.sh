export LD_LIBRARY_PATH="./:/usr/local/lib/"
g++ -O3 -funroll-loops -o kyber3-enumbs-parallel-est-gate ./NIST-round3/kyber3-enumbs-parallel-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./kyber3-enumbs-parallel-est-gate
