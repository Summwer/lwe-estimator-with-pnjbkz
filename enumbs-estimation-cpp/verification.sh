export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
g++ -O3 -funroll-loops -o verify_enumbs ./NIST-round3/verify_enumbs.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./verify_enumbs
