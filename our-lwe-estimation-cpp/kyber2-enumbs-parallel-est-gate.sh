export LD_LIBRARY_PATH="enumbs-estimation-cpp"
g++ -O3 -funroll-loops -o kyber2-enumbs-parallel-est-gate ./NIST-round3/kyber2-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./kyber2-enumbs-parallel-est-gate
