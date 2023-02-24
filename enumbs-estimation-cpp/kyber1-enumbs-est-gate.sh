export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
g++ -O3 -funroll-loops -o kyber1-enumbs-est-gate ./NIST-round3/kyber1-enumbs-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./kyber1-enumbs-est-gate
