export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"

mkdir "svpchal-est-result"


# params input in main function in enumbs-est
# argv[0]: implemented file name
# argv[1]: jump value J
# argv[2]: max_loop value
# argv[3]: cost model: 1.gate model 2.practical sec model
# argv[4]: maximal dimension in enumeration
# argv[5]: threads number
# argv[6]: practical_pump_d4f
# argv[7]: start beta value.
# argv[8]: maximal RAM memory.
g++ -O3 -funroll-loops -o svpchal-enumbs-est ./svp-est/svpchal-enumbs-parallel-est.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest

# ./lwechal-enumbs-est 100 5 2 300 2 3 50 0 | tee lwechal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=100-maxloop=5.log" 


# ./svpchal-enumbs-est 100 1 2 300 2 3 50 0 | tee svpchal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=100-maxloop=1.log" 


./svpchal-enumbs-est 100 1 2 300 5 3 50 2  | tee svpchal-est-result/"enumbs(32+2gpus)-d4f-default-g6k-J=100-maxloop=1.log" 



