export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
mkdir bssa_result
g++ -O3 -funroll-loops -o lwechal-bssa-est ./NIST-round3/lwechal-bssa-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest




# params input in main function
# argv[0]: implemented file name
# argv[1]: jump value J
# argv[2]: max_loop value
# argv[3]: cost model: 1.gate model 2.practical sec model
# argv[4]: maximal dimension in enumeration
# argv[5]: enumbs_min_G: 0-- false, find minimal RAM strategy; 1--true, find minimal time cost strategy
# argv[6]: max_RAM
# argv[7]: practical_pump_d4f


# cd "bssa_result"
# mkdir "beta-start-50-last-pump-theod4f1"
# cd ..


# ./lwechal-bssa-est 8 5 2 300 1 0 1 | tee bssa_result/beta-start-50-last-pump-theod4f1/lwechal-bssa-est-sec-J-8-loop-5.log

# ./lwechal-bssa-est 8 5 2 300 0 44 1 | tee bssa_result/beta-start-50-last-pump-theod4f1/lwechal-bssa-est-sec-J-8-loop-5-limRAM.log

# ./lwechal-bssa-est 8 1 2 300 1 0 1 | tee bssa_result/beta-start-50-last-pump-theod4f1/lwechal-bssa-est-sec-J-8-loop-1.log

# ./lwechal-bssa-est 1 5 2 300 1 0 1 | tee bssa_result/beta-start-50-last-pump-theod4f1/lwechal-bssa-est-sec-J-1-loop-5.log

# ./lwechal-bssa-est 1 1 2 300 1 0 1 | tee bssa_result/beta-start-50-last-pump-theod4f1/lwechal-bssa-est-sec-J-1-loop-1.log

cd "bssa_result"
mkdir "beta-start-50-last-pump-theod4f2"
cd ..


./lwechal-bssa-est 8 5 2 300 1 0 2 | tee bssa_result/beta-start-50-last-pump-theod4f2/lwechal-bssa-est-sec-J-8-loop-5.log

./lwechal-bssa-est 8 5 2 300 0 44 2 | tee bssa_result/beta-start-50-last-pump-theod4f2/lwechal-bssa-est-sec-J-8-loop-5-limRAM.log

# ./lwechal-bssa-est 8 1 2 300 1 0 2 | tee bssa_result/beta-start-50-last-pump-theod4f2/lwechal-bssa-est-sec-J-8-loop-1.log

./lwechal-bssa-est 1 5 2 300 1 0 2 | tee bssa_result/beta-start-50-last-pump-theod4f2/lwechal-bssa-est-sec-J-1-loop-5.log

# ./lwechal-bssa-est 1 1 2 300 1 0 2 | tee bssa_result/beta-start-50-last-pump-theod4f2/lwechal-bssa-est-sec-J-1-loop-1.log



cd "bssa_result"
mkdir "beta-start-50-last-pump-default-g6k"
cd ..




# ./lwechal-bssa-est 8 5 2 300 1 0 3 | tee bssa_result/beta-start-50-last-pump-default-g6k/lwechal-bssa-est-sec-J-8-loop-5.log

# ./lwechal-bssa-est 8 5 2 300 0 44 3 | tee bssa_result/beta-start-50-last-pump-default-g6k/lwechal-bssa-est-sec-J-8-loop-5-limRAM.log

# ./lwechal-bssa-est 8 1 2 300 1 0 3 | tee bssa_result/beta-start-50-last-pump-default-g6k/lwechal-bssa-est-sec-J-8-loop-1.log

./lwechal-bssa-est 1 5 2 300 1 0 3 | tee bssa_result/beta-start-50-last-pump-default-g6k/lwechal-bssa-est-sec-J-1-loop-5.log

# ./lwechal-bssa-est 1 1 2 300 1 0 3 | tee bssa_result/beta-start-50-last-pump-default-g6k/lwechal-bssa-est-sec-J-1-loop-1.log