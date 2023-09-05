# Guidance
`lwe-estimator-with-pnjbkz` provides functions for estimating the concrete security of Learning with Errors instances and simulator implementations.


Code in `cpp` implements the EnumBS algorithm, and the security estimation of NIST-round3 schemes and lwe challenges(https://www.latticechallenge.org/lwe_challenge/challenge.php) in EnumBS and BSSA. The guidance in implementing the code and the organization of the code please follows `README.md` in `cpp`.

Code in `sage` implements default g6k estimation, bkz-only estimation, leaky-lwe-estimator and two-step estimation for NIST-round3 schemes(including Kyber, Dilithium). It should be implemented in a sage 9.1+ environment. The guidance in implementing the code and the organization of the code please follows `README.md` in `sage`.


All of the estimators are implemented according to the  (nist)-bkz estimator implemented in leaky-lwe-estimator regarding the norm of shortest vector (e,1) for LWE instance in primal attack as a chi-squared distribution, and computing a success probability in it. Besides, we also implement a variant based 2016 estimate for two-step mode in sage.

Let `worst_case = True` and `cost_model = 1`, then it is the form for computing success probability, we use this estimator to estimate the NIST schemes, corresponding to the values of Two-step column shown in Table 7 in https://eprint.iacr.org/2022/1343.pdf. Set `worst_case = False` and `cost_model = 2`, then it is an estimation used for estimating the solvability TU LWE challenge, corresponding to the Table 3 in https://eprint.iacr.org/2022/1343.pdf.

More details about the above algorithims can be seen in the artical https://eprint.iacr.org/2022/1343.pdf.
