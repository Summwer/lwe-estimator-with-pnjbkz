export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
g++ -O3 -funroll-loops -o lwechal-bssa-est-sec ./NIST-round3/lwechal-bssa-est-sec.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-bssa-est-sec
