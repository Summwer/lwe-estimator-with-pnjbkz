export LD_LIBRARY_PATH="./:/usr/bin/" #"./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/local/bin/:/usr/include/:/usr/lib/"

mkdir "mimabisai-est"

# params input in main function in enumbs-est
# argv[0]: implemented file name
# argv[1]: method: 1 --enumbs; 2 --bssa.
# argv[2]: list-decoding complexity 1: "agps20"; 2: "matzov22"
# argv[3]: maxloop: >=1 and is an integer
# argv[4]: maximal jump value
g++ -O3 -funroll-loops -o mimabisai-est.out ./mimabisai-est/mimabisai-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest

./mimabisai-est.out 1 2 1 100 | tee mimabisai-est/"mimabisai-est.log"