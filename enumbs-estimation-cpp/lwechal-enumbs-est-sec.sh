export LD_LIBRARY_PATH="./:/usr/local/lib/"
g++ -O3 -funroll-loops -o lwechal-enumbs-est-sec ./NIST-round3/lwechal-enumbs-est-sec.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-enumbs-est-sec
