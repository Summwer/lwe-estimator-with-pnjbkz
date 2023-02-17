export LD_LIBRARY_PATH="/home/cryptothesis/summer/debug/our-lwe-estimation-cpp"

g++ -O3 -funroll-loops -o dilithium1-enumbs-parallel-est-gate ./NIST-round3/dilithium1-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./dilithium1-enumbs-parallel-est-gate | tee enumbs_result/dilithium1-enumbs-parallel-est-gate.log


g++ -O3 -funroll-loops -o dilithium2-enumbs-parallel-est-gate ./NIST-round3/dilithium2-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./dilithium2-enumbs-parallel-est-gate | tee enumbs_result/dilithium2-enumbs-parallel-est-gate.log

g++ -O3 -funroll-loops -o dilithium3-enumbs-parallel-est-gate ./NIST-round3/dilithium3-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./dilithium3-enumbs-parallel-est-gate | tee enumbs_result/dilithium3-enumbs-parallel-est-gate.log


g++ -O3 -funroll-loops -o kyber1-enumbs-parallel-est-gate ./NIST-round3/kyber1-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./kyber1-enumbs-parallel-est-gate | tee enumbs_result/kyber1-enumbs-parallel-est-gate.log


g++ -O3 -funroll-loops -o kyber2-enumbs-parallel-est-gate ./NIST-round3/kyber2-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./kyber2-enumbs-parallel-est-gate | tee enumbs_result/kyber2-enumbs-parallel-est-gate.log

g++ -O3 -funroll-loops -o kyber3-enumbs-parallel-est-gate ./NIST-round3/kyber3-enumbs-parallel-est-gate.cpp -pthread -lfplll -lgmp -lmpfr -lest
./kyber3-enumbs-parallel-est-gate | tee enumbs_result/kyber3-enumbs-parallel-est-gate.log