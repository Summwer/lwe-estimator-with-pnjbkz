import time

load("../framework/est/NIST-est.sage")


######################################
#set method parameters
'''
    
    method = 1: progressive bkz estimation.
             2: two step mode
             3: bkz-only with fixed blocksize
             4: default g6k
             "fixed-blocksize bkz-only": call fixed-blocksize bkz-only model to estimate security.
    progressive_sieve = True:  progressive sieve

    cost_model = 1: theoretical cost estimation
               = 2: experimental cost estimation
    :cumG : False: [PV21] estimate
            True: cumulated gate estimate adaptive to [PV21]
'''
#######################################
#Fixed parameters

method = 1
ldc_param = "MATZOV22" #change list decoding parameters

#------------------------------------
#[PV21] estimation
cumG = True #False #false: average beta; true: average G
kyber_est(method,cumG = cumG, ldc_param = ldc_param )
dilithium_est(method,cumG = cumG, ldc_param = ldc_param )
#hufu_est(method,cumG = cumG, ldc_param = ldc_param )
#eaglesign_est(method,cumG = cumG, ldc_param = ldc_param )
#haetae_est(method,cumG = cumG, ldc_param = ldc_param )

