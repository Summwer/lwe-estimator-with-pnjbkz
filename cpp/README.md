# Guidance

This project implements the EnumBS algorithm, and the security estimation of NIST-round3 schemes and lwe challenges(https://www.latticechallenge.org/lwe_challenge/challenge.php) in EnumBS and BSSA. 

All of the estimators are implemented according to the  (nist)-bkz estimator implemented in leaky-lwe-estimator regarding the norm of shortest vector (e,1) for LWE instance in primal attack as a chi-squared distribution, and computing a success probability in it. Let `worst_case = True` and `cost_model = 1`, then it is the form for computing success probability, we use this estimator to estimate the NIST schemes, corresponding to the values of Two-step column shown in Table 7 in https://eprint.iacr.org/2022/1343.pdf. Set `worst_case = False` and `cost_model = 2`, then it is an estimation used for estimating the solvability TU LWE challenge and it returns a blocksize strategy for solving TU LWE challenge by EnumBS or BSSA, corresponding to the Table 3 in https://eprint.iacr.org/2022/1343.pdf.


### Dependencies

#### Required

* GNU MP 4.2.0 or higher http://gmplib.org/ or MPIR 1.0.0 or higher [http://mpir.org](http://mpir.org/)
* MPFR 2.3.0 or higher, COMPLETE INSTALLATION http://www.mpfr.org/
* autotools 2.61 or higher
* g++ 4.9.3 or higher
* boost 7.4 or higher : https://www.boost.org/

In unbuntu, you can install the above dependences by implementing the following code directly:

```bash
./install-dependency.sh
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


6. download and install fplll library(https://github.com/fplll/fplll) in the folder `cpp`: 

```bash
git clone https://github.com/fplll/fplll.git
cd fplll
./autogen.sh
./configure
make
make install
```


### Organization of the code
There are 3 code folders in `cpp`.  `framework` contains core code for EnumBS, BSSA, pnj-BKZ simulator and cost model implementations. If developers modify the codes in it, please run `./rebuild.sh` in the main direction. 

There are two files in the folder `lwe-est`: 
- `lwechal-bssa-est.cpp` is an executable file for generating blocksize strategy by BSSA algorithm and give a cost estimation for lwe instances provided in TU LWE challenge(https://www.latticechallenge.org/lwe_challenge/challenge.php). One can run it by running code 

- `lwechal-bssa-est.cpp` is an executable file for generating blocksize strategy by EnumBS algorithm and give a cost estimation for lwe instances provided in TU LWE challenge(https://www.latticechallenge.org/lwe_challenge/challenge.php). 

One can test the blocksize strategy generation method in EnumBS and BSSA by running code 
```
./lwechal-est.sh
```
in the main directory. It will return the blocksize strategy for some of LWE challenges in the folder `lwechal-est-result`.


The code file `NIST-round3-est-gate.cpp` in `NIST-round3-est` is used for estimating all the LWE-based NIST schemes(Kyber, Dilithium and Frodo) by EnumBS estimator. One can test it by running the code 
```
./implement_all_NIST_schemes.sh
```
in the main directory. It will return the blocksize strategy for some of LWE challenges in the folder `lwechal-est-result`.

