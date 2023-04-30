export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"

mkdir bssa_result

# ./lwechal-enumbs-parallel-est-sec.sh

g++ -O3 -funroll-loops -o NIST-round3-bssa-est-gate ./NIST-round3/NIST-round3-bssa-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./NIST-round3-bssa-est-gate | tee bssa_result/NIST-round3-bssa-est.log

