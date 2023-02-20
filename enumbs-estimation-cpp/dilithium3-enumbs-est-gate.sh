export LD_LIBRARY_PATH="./:/usr/local/lib/"
g++ -O3 -funroll-loops -o dilithium3-enumbs-est-gate./NIST-round3/dilithium3-enumbs-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./dilithium3-enumbs-est-gate
