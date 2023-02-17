export LD_LIBRARY_PATH="/home/cryptothesis/summer/debug/our-lwe-estimation-cpp"
g++ -O3 -funroll-loops -o lwechal-enumbs-parallel-est-sec ./NIST-round3/lwechal-enumbs-parallel-est-sec.cpp -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-enumbs-parallel-est-sec
