export LD_LIBRARY_PATH="enumbs-estimation-cpp"
g++ -O3 -funroll-loops -o verify_enumbs ./NIST-round3/verify_enumbs.cpp -pthread -lfplll -lgmp -lmpfr -lest
./verify_enumbs
