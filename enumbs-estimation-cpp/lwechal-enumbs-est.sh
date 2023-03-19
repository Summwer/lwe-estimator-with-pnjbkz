export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
g++ -O3 -funroll-loops -o lwechal-enumbs-est ./NIST-round3/lwechal-enumbs-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-enumbs-est | tee enumbs_result/lwechal-enumbs-est-sec-loop-30.log
