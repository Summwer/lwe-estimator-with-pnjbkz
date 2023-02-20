apt-get update
apt-get install libmpfr-dev #ubuntu
apt-get install m4 #should install m4 first
apt-get install libgmp-dev
apt-get install libboost-all-dev
apt-get install autoconf

git clone https://github.com/fplll/fplll.git
cd fplll
./autogen.sh
./configure
make
make install
make check
cd ..


./compile-framework.sh
