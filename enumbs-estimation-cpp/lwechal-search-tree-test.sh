export LD_LIBRARY_PATH="./:/usr/local/lib/:/usr/lib/x86_64-linux-gnu/:/usr/bin/"
g++ -O3 -funroll-loops -o test_search_tree test_search_tree.cpp -L. -pthread -lfplll -lgmp -lmpfr -lest
./test_search_tree
