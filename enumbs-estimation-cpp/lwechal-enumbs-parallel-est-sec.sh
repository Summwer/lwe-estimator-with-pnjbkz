export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
g++ -O3 -funroll-loops -o lwechal-enumbs-parallel-est-sec ./NIST-round3/lwechal-enumbs-parallel-est-sec.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-enumbs-parallel-est-sec
