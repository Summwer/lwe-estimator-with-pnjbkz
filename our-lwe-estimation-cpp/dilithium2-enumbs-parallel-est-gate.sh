export LD_LIBRARY_PATH="enumbs-estimation-cpp"
g++ -O3 -funroll-loops -o dilithium2-enumbs-parallel-est-gate./NIST-round3/dilithium2-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./dilithium2-enumbs-parallel-est-gate
