# Guidance





### Dependences

#### Required

* GNU MP 4.2.0 or higher http://gmplib.org/ or MPIR 1.0.0 or higher [http://mpir.org](http://mpir.org/)
* MPFR 2.3.0 or higher, COMPLETE INSTALLATION http://www.mpfr.org/
* autotools 2.61 or higher
* g++ 4.9.3 or higher
* boost 7.4 or higher : https://www.boost.org/

In unbuntu, you can install the above dependences by implementing the following code directly:

```bash
./install-dependencies.sh
```



Or install the dependencies on yourserlf: 

1. C++14

2. mpfr: https://www.mpfr.org/

```bash
apt-get update
apt-get install libmpfr-dev #ubuntu
```

3. gmp: https://gmplib.org/

```bash
apt-get install m4 #should install m4 first
apt-get install libgmp-dev
```

4. boost: 

```bash
apt-get install libboost-all-dev
```

5. autoconf

```
apt-get install autoconf
```



6. download and install fplll library(https://github.com/fplll/fplll) in the folder `enumbs-estimation-cpp`: 

```bash
git clone https://github.com/fplll/fplll.git
cd fplll
./autogen.sh
./configure
make
```



### Compile codes

```
./compile-framework.sh

```





