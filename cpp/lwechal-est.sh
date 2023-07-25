export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"

mkdir "lwechal-est-result"

# params input in main function in enumbs-est
# argv[0]: implemented file name
# argv[1]: jump value J
# argv[2]: max_loop value
# argv[3]: cost model: 1.gate model 2.practical sec model
# argv[4]: maximal dimension in enumeration
# argv[5]: threads number
# argv[6]: practical_pump_d4f
g++ -O3 -funroll-loops -o lwechal-enumbs-est ./lwe-est/lwechal-enumbs-parallel-est.cpp  -L. -pthread -lfplll -lgmp -lmpfr -lest

./lwechal-enumbs-est 100 10 2 300 10 3 | tee lwechal-est-result/enumbs-est.log



# # params input in main function in bssa-est
# # argv[0]: implemented file name
# # argv[1]: jump value J
# # argv[2]: max_loop value
# # argv[3]: cost model: 1.gate model 2.practical sec model
# # argv[4]: maximal dimension in enumeration
# # argv[5]: enumbs_min_G: 0-- false, find minimal RAM strategy; 1--true, find minimal time cost strategy
# # argv[6]: max_RAM
# # argv[7]: practical_pump_d4f
g++ -O3 -funroll-loops -o lwechal-bssa-est ./lwe-est/lwechal-bssa-est.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest

./lwechal-bssa-est 100 10 2 300 1 0 3 1 | tee lwechal-est-result/bssa-est.log




