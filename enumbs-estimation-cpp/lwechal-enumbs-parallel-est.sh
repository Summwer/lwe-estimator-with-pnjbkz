export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
g++ -O3 -funroll-loops -o lwechal-enumbs-parallel-est ./NIST-round3/lwechal-enumbs-parallel-est.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest


# params input in main function
# argv[0]: implemented file name
# argv[1]: jump value J
# argv[2]: max_loop value
# argv[3]: cost model: 1.gate model 2.practical sec model
# argv[4]: maximal dimension in enumeration
# argv[5]: threads number
# argv[6]: practical_pump_d4f


# cd "enumbs_result"
# mkdir "beta-start-50-last-pump-theod4f1"
# cd ..


# ./lwechal-enumbs-parallel-est 8 5 2 300 10 1 #| tee enumbs_result/beta-start-50-last-pump-theod4f1/lwechal-enumbs-sec-prec3-J-8-loop-5.log

# ./lwechal-enumbs-parallel-est 8 1 2 300 10 1 | tee enumbs_result/beta-start-50-last-pump-theod4f1/lwechal-enumbs-sec-prec3-J-8-loop-1.log

# ./lwechal-enumbs-parallel-est 1 5 2 300 10 1 | tee enumbs_result/beta-start-50-last-pump-theod4f1/lwechal-enumbs-sec-prec3-J-1-loop-5.log

# ./lwechal-enumbs-parallel-est 1 1 2 300 10 1 | tee enumbs_result/beta-start-50-last-pump-theod4f1/lwechal-enumbs-sec-prec3-J-1-loop-1.log

# cd "enumbs_result"
# mkdir "beta-start-50-last-pump-theod4f2"
# cd ..

# ./lwechal-enumbs-parallel-est 8 5 2 300 10 2 #| tee enumbs_result/beta-start-50-last-pump-theod4f2/lwechal-enumbs-sec-prec3-J-8-loop-5.log

# ./lwechal-enumbs-parallel-est 8 1 2 300 10 2 | tee enumbs_result/beta-start-50-last-pump-theod4f2/lwechal-enumbs-sec-prec3-J-8-loop-1.log

# ./lwechal-enumbs-parallel-est 1 5 2 300 10 2 | tee enumbs_result/beta-start-50-last-pump-theod4f2/lwechal-enumbs-sec-prec3-J-1-loop-5.log

# ./lwechal-enumbs-parallel-est 1 1 2 300 10 2 | tee enumbs_result/beta-start-50-last-pump-theod4f2/lwechal-enumbs-sec-prec3-J-1-loop-1.log






# cd "enumbs_result"
# mkdir "beta-start-50-last-pump-default-g6k"
# cd ..

./lwechal-enumbs-parallel-est 8 5 2 300 10 3 #| tee enumbs_result/beta-start-50-last-pump-default-g6k/lwechal-enumbs-sec-prec3-J-8-loop-5.log

# ./lwechal-enumbs-parallel-est 8 1 2 300 10 3 | tee enumbs_result/beta-start-50-last-pump-default-g6k/lwechal-enumbs-sec-prec3-J-8-loop-1.log

# ./lwechal-enumbs-parallel-est 1 5 2 300 10 3 | tee enumbs_result/beta-start-50-last-pump-default-g6k/lwechal-enumbs-sec-prec3-J-1-loop-5.log

# ./lwechal-enumbs-parallel-est 1 1 2 300 10 3 #| tee enumbs_result/beta-start-50-last-pump-default-g6k/lwechal-enumbs-sec-prec3-J-1-loop-1.log
