# Guidance
Code in `sage` implements default g6k estimation, bkz-only estimation, leaky-lwe-estimator and two-step estimation for LWE challenge and NIST-round3 schemes(including Kyber, Dilithium and Frodo). It should be implemented in a sage 9.1+ environment.

All of the estimators are implemented according to the  (nist)-bkz estimator implemented in leaky-lwe-estimator regarding the norm of shortest vector (e,1) for LWE instance in primal attack as a chi-squared distribution, and computing a success probability in it. Besides, we also implement a variant based 2016 estimate for two-step mode in sage. Let `worst_case = True` and `cost_model = 1`, then it is the form for computing success probability, we use this estimator to estimate the NIST schemes, corresponding to the values of Two-step column shown in Table 7 in https://eprint.iacr.org/2022/1343.pdf. Set `worst_case = False` and `cost_model = 2`, then it is an estimation used for the solvability TU LWE challenge, corresponding to the Table 3 in https://eprint.iacr.org/2022/1343.pdf.


### Organization of the code
There are 3 code folders in `sage`.  `framework` contains core code for the (nist)-bkz estimator implemented in leaky-lwe-estimator(i.e. simple progressive BKZ estiamator, `progressive_bkz_est.sage`), default g6k estimator(`default_g6k_est.sage`), BKZ-only estimator(`bkz_only_est.sage`) and two-step estimator(`two_step_mode_est.sage`) in `est` sub-folder, pnj-BKZ simulator(`pnjbkz_simulator.sage`) and pump simulator(use for simulate default g6k process, `pump_simulator.sage`) used for BKZ and pump simulation in `simulator` sub-folder. Besides, we also the implement the lwe instance geneneration from TU LWE challenge and the cost model in `framework` folder. 

There are 4 files in the folder `lwechal-est`, they give an estimation for the lwe instance from TU LWE challenge in default g6k, BKZ-only, simple progressive BKZ and two-step mode using practical cost model(32threads + 2 nvidia 3090 gpus) declared in https://eprint.iacr.org/2022/1343.pdf.
- `lwechal-bkz-only.sage` is an executable file for BKZ-only estimation on TU LWE challenge. One can run it in sagemath:
```
sage
cd lwechal-est
load("lwechal-bkz-only.sage")
```
- `lwechal-default-g6k.sage` estimates the time and memory cost for LWE instance in default g6k mode. One can run it in sagemath:
```
sage
cd lwechal-est
load("lwechal-default-g6k.sage")
```
- `lwechal-pro-bkz.sage` estimates the time and memory cost for LWE instance in the simple progressive BKZ mode proposed in https://eprint.iacr.org/2020/292.pdf, i.e. increase blocksize one by one for BKZ. One can run it in sagemath:
```
sage
cd lwechal-est
load("lwechal-pro-bkz.sage")
```
- `lwechal-two-step.sage` estimates the time and memory cost for LWE instance in the two-step mode (simple progressive BKZ + last pump). One can run it in sagemath:
```
sage
cd lwechal-est
load("lwechal-two-step.sage")
```


The code files in `NIST-round3` is used for estimating all the LWE-based NIST schemes(Kyber, Dilithium and Frodo) by EnumBS estimator. 

- `NIST-pro-bkz.sage` calls a progressive BKZ and chi-squared distribution to estiamte NIST schemes, it is the adaptive estimation method mentioned in Kyber(round3), i.e. leaky-lwe-estimator[DDGR20] (We implement it in average cost, not average beta). One can test it by running the code 
```
sage
cd NIST-round3
load("NIST-pro-bkz.sage")
```
- `NIST-two-step.sage` calls a two-step mode(progressive BKZ + last pump) and chi-squared distribution to estiamte NIST schemes, it is the estimation method mentioned in Kyber(round3). One can test it by running the code 
```
sage
cd NIST-round3
load("NIST-two-step.sage")
```
If let `ldc_param = "MATZOV22"`, then we can change the classical list decoding cost estimation from which given in [AGPS20] to [MATZOV22].s

LWE-estimator used some conservative assumptions in two-step mode please run the code 
```
sage
cd NIST-round3
load("NIST-two-step-martin-primal-usvp.sage")
```


More detail about the above estimators, please see the article https://eprint.iacr.org/2022/1343.pdf .


Estimations for lwe challenge is in the path "sage/lwechal/"
Estimations for NIST-round3 schemes is in the path "sage/NIST-round3/", the dimension and dvol(ln(vol(lattice))) for NIST-round3 schemes are selected from the open source https://github.com/lducas/leaky-LWE-Estimator.





# The experiment in Refined-LWE-Estimator
## Figure 2
One can run the following command in the folder `sage/lwechal-est` and obatain the two-step estimation of LWE Challenge and its success probabibility during the process. It constructs the `this work(GSA for LLL)` in Fig.2. and the `Two-step + [AGPS20]`, `Two-step + [MATZOV22]` in Fig.2(a) and Fig.2(b).
```
load("lwechal-two-step.sage")    
```



## Table 2

To generate the data of Table 2, one can run the file `Lower_Bound_Estimation.py` in the folder `sage/NIST-round3/`. 


## Column "Previous" and  $S_0$ of Table 1
To obtain the result of column "Previous", one should run the file `NIST-pro-bkz.sage` in the folder `sage/NIST-round3` of [https://github.com/Summwer/lwe-estimator-with-pnjbkz/tree/refined-lwe-estimator](https://github.com/Summwer/lwe-estimator-with-pnjbkz/tree/refined-lwe-estimator). 

To obtain the result of column $S_0$, one should run the file `NIST-two-step.sage` in the folder `sage/NIST-round3` of [https://github.com/Summwer/lwe-estimator-with-pnjbkz/tree/refined-lwe-estimator](https://github.com/Summwer/lwe-estimator-with-pnjbkz/tree/refined-lwe-estimator). 

