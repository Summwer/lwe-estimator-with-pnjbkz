export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
g++ -O3 -funroll-loops -o lwechal-enumbs-parallel-est ./NIST-round3/lwechal-enumbs-parallel-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest

./lwechal-enumbs-parallel-est | tee enumbs_result/beta-start-50/lwechal-enumbs-parallel-est-sec-leaf-prec3-limRAM.log
# ./lwechal-enumbs-parallel-est | tee enumbs_result/beta-start-50/lwechal-enumbs-parallel-est-sec-leaf-prec3.log

