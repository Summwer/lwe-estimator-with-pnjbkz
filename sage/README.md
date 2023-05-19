# Guidance


Code in `sage` implements default g6k estimation, bkz-only estimation, pro-bkz estimation and two-step estimation. It should be implemented in a sage 9.1+ environment.

default g6k estimation: estimate the time and memory cost for approximate svp instance in default g6k mode.
bkz-only estimation: estimate the time and memory cost for approximate svp instance in bkz-only mode with a fixed bkz blocksize.
pro-bkz estimation: estimate the time and memory cost for approximate svp instance in a progressive bkz only mode.
two-step estimation: estimate the time and memory cost for approximate svp instance in two-step mode: progressive bkz + last pump.

More detail about the above estimation, please see the article https://eprint.iacr.org/2022/1343.pdf .


Estimations for lwe challenge is in the path "sage/lwechal/"
Estimations for NIST-round3 schemes is in the path "sage/NIST-round3/", the dimension and dvol(ln(vol(lattice))) for NIST-round3 schemes are selected from the open source https://github.com/lducas/leaky-LWE-Estimator.
