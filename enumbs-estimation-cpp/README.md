# Guidance





### Library Dependences

To implement our EnumBS estimation code in cpp, we should first install the following dependences:

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

5. download and install fplll library(https://github.com/fplll/fplll) in the folder `enumbs-estimation-cpp`: 

```bash

git clone https://github.com/fplll/fplll.git
```

