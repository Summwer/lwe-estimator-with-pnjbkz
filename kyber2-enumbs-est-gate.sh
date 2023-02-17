export LD_LIBRARY_PATH="/home/cryptothesis/summer/debug/our-lwe-estimation-cpp/fplll/.lib:/home/cryptothesis/summer/debug/our-lwe-estimation-cpp"
g++ -O3 -funroll-loops -o kyber2-enumbs-est-gate ./NIST-round3/kyber2-enumbs-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./kyber2-enumbs-est-gate
