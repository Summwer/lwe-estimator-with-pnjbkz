export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
mkdir bssa_result
g++ -O3 -funroll-loops -o lwechal-bssa-est ./NIST-round3/lwechal-bssa-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./lwechal-bssa-est | tee bssa_result/lwechal-bssa-est-sec-loop-15-worst-case.log
