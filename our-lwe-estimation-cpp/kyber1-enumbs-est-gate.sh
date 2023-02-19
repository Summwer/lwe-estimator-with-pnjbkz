export LD_LIBRARY_PATH="enumbs-estimation-cpp"
g++ -O3 -funroll-loops -o kyber1-enumbs-est-gate ./NIST-round3/kyber1-enumbs-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./kyber1-enumbs-est-gate
