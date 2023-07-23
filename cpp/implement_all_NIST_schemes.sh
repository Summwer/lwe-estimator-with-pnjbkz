export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"

mkdir "nist-round3-est-result"

# params input in main function in enumbs-est
# argv[0]: implemented file name
# argv[1]: method: 1 --enumbs; 2 --bssa.
# argv[2]: G2_prec =  1./argv[2]
g++ -O3 -funroll-loops -o nist-round3-est ./NIST-round3-est/NIST-round3-est-gate.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest

./nist-round3-est 1 1000 | tee nist-round3-est-result/enumbs-est2-prec3.log