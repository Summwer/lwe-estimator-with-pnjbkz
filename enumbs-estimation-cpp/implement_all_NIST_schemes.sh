export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"

mkdir enumbs_result

# ./lwechal-enumbs-parallel-est-sec.sh

g++ -O3 -funroll-loops -o NIST-round3-enumbs-parallel-est-gate ./NIST-round3/NIST-round3-enumbs-parallel-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./NIST-round3-enumbs-parallel-est-gate | tee enumbs_result/NIST-round3-enumbs-parallel-est-gate-add3.log

