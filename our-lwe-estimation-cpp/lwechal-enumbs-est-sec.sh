export LD_LIBRARY_PATH="/home/cryptothesis/summer/debug/our-lwe-estimation-cpp"
g++ -O3 -funroll-loops -o lwechal-enumbs-est-sec ./NIST-round3/lwechal-enumbs-est-sec.cpp -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-enumbs-est-sec
