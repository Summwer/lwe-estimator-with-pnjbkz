export LD_LIBRARY_PATH="enumbs-estimation-cpp"
g++ -O3 -funroll-loops -o dilithium3-enumbs-est-gate./NIST-round3/dilithium3-enumbs-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./dilithium3-enumbs-est-gate
