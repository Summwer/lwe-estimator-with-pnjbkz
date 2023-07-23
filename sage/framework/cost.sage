# Data for the cost of sieving from:
# Estimating quantum speedups for lattice sieves
# Martin R. Albrecht and Vlad Gheorghiu and Eamonn W. Postlethwaite and John M. Schanck
# The data used is from the following article:
# https://eprint.iacr.org/2019/1161 
# This datafile can be extracted from the pdf via the linux tool `pdfdetach'.
# We use the version of Sep 2020
# https://github.com/jschanck/eprint-2019-1161/blob/main/data/cost-estimate-list_decoding-classical.csv

from mpmath import mp
from math import ceil, floor, exp, pi, pow
from math import log as ln


def log2(x):
    return ln(x)/ln(2.)

agps20_gate_data = {
          64  :42.5948446291284,   72  :44.8735917172503, 80  :47.4653141889341, 88  :50.0329479433691, 96  :52.5817667347844,  
          104 :55.1130237325179,   112 :57.6295421450947, 120 :60.133284108578,  128 :62.1470129451821, 136 :65.4744488064273, 
          144 :67.951405476229,    152 :70.0494944191399, 160 :72.50927387359,   168 :74.9619105412039, 176 :77.4100782579645, 
          184 :79.3495443657483,   192 :81.7856479853679, 200 :84.2178462414349, 208 :86.646452845262,  216 :89.0717383389617, 
          224 :91.4939375786565,   232 :93.9132560751063, 240 :96.3298751307529, 248 :98.7439563146036, 256 :101.155644837658, 
          264 :104.091650357302,   272 :106.500713866161, 280 :108.907671199501, 288 :111.312627066864, 296 :113.715679081585, 
          304 :116.11691871212,    312 :118.516432037545, 320 :120.914300351043, 328 :123.310600632063, 336 :125.705405925853, 
          344 :128.098785623819,   352 :130.490805751072, 360 :132.881529104042, 368 :135.271015458153, 376 :137.659321707881, 
          384 :140.046501985502,   392 :142.432607773479, 400 :144.817688009257, 408 :147.201789183958, 416 :149.584955436701, 
          424 :151.967228645918,   432 :154.348648518547, 440 :156.729252677678, 448 :159.109076748918, 456 :161.488154445581, 
          464 :163.866517652676,   472 :166.24419650959,  480 :168.621219491327, 488 :170.997613488119, 496 :173.373403883249, 
          504 :175.748614628914,   512 :178.123268319974, 520 :180.931640474467, 528 :183.305745118107, 536 :185.679338509895, 
          544 :188.052439374005,   552 :190.425065356218, 560 :192.797233085084, 568 :195.168958230518, 576 :197.540255559816, 
          584 :199.911138991095,   592 :202.281621644196, 600 :204.651715889082, 608 :207.02143339179,  616 :209.390785157985, 
          624 :211.759781574203,   632 :214.128432446848, 640 :216.496747039019, 648 :218.864734105257, 656 :221.232401924303, 
          664 :223.599758329925,   672 :225.96681073994,  680 :228.333566183483, 688 :230.700031326626, 696 :233.066212496418, 
          704 :235.43211570344,    712 :237.797746662944, 720 :240.163110814653, 728 :242.528213341298, 736 :244.893059185964, 
          744 :247.25765306831,    752 :249.621999499728, 760 :251.986102797502, 768 :254.349967098032, 776 :256.71359636917,  
          784 :259.076994421734,   792 :261.440164920231, 800 :263.803111392861, 808 :266.165837240825, 816 :268.528345816343, 
          824 :270.890640143248,   832 :273.252723321704, 840 :275.614598434176, 848 :277.976268306208, 856 :280.337735739304, 
          864 :282.699003457275,   872 :285.060074111424, 880 :287.420950285349, 888 :289.781634499399, 896 :292.142129214795, 
          904 :294.502436837451,   912 :296.862559721505, 920 :299.222500172584, 928 :301.582260450819, 936 :303.941842773632, 
          944 :306.301249318305,   952 :308.660482224348, 960 :311.019543595679, 968 :313.378435502636, 976 :315.737159983825, 
          984 :318.095719047813,   992 :320.454114674691, 1000:322.8123488175,   1008:325.170423403542,1016:327.52834033558, 
          1024:329.886101492934
          }
"""
agps20_gate_data = {64: 43.211883,80: 48.457505,96: 53.596866,112: 58.664746,128: 63.676389,144: 68.179263,160: 73.941003,176: 78.495338,192: 83.386586,208: 87.765161,224: 92.619433,240: 97.461642,256: 102.293241,272: 107.630129,288: 112.447623,304: 117.257149,320: 122.059455,336: 126.855209,352: 131.645008,368: 136.429391,384: 141.208847,400: 145.983818,416: 150.754701,432: 155.521855,448: 160.285602,464: 165.046233,480: 169.804005,496: 174.559148,512: 179.311870,528: 184.489459,544: 189.238986,560: 193.986524,576: 198.732208,592: 203.476157,608: 208.218480,624: 212.959270,640: 217.698612,656: 222.436582,672: 227.173249,688: 231.908671,704: 236.642905,720: 241.375999,736: 246.107999,752: 250.838945,768: 255.568875,784: 260.297823,800: 265.025820,816: 269.752895,832: 274.479077,848: 279.204391,864: 283.928860,880: 288.652508,896: 293.375355,912: 298.097424,928: 302.818732,944: 307.539300,960: 312.259145,976: 316.978285,992: 321.696736,1008: 326.414517,1024: 331.131641,1040: 336.355357,1056: 341.071849,1072: 345.787718,1088: 350.502979,1104: 355.217645,1120: 359.931730,1136: 364.645249,1152: 369.358216,1168: 374.070643,1184: 378.782544,1200: 383.493932,1216: 388.204818,1232: 392.915217,1248: 397.625139,1264: 402.334597,1280: 407.043602,1296: 411.752167,1312: 416.460301,1328: 421.168016,1344: 425.875324,1360: 430.582234,1376: 435.288756,1392: 439.994902,1408: 444.700680,1424: 449.406100,1440: 454.111173,1456: 458.815905,1472: 463.520308,1488: 468.224390,1504: 472.928159,1520: 477.631623,1536: 482.334791,1552: 487.037671,1568: 491.740270,1584: 496.442595,1600: 501.144654,1616: 505.846455,1632: 510.548003,1648: 515.249306,1664: 519.950370,1680: 524.651202,1696: 529.351807,1712: 534.052192,1728: 538.752362,1744: 543.452323,1760: 548.152080,1776: 552.851639,1792: 557.551004,1808: 562.250182,1824: 566.949176,1840: 571.647992,1856: 576.346633,1872: 581.045105,1888: 585.743410,1904: 590.441555,1920: 595.139542,1936: 599.837375,1952: 604.535058,1968: 609.232595,1984: 613.929990,2000: 618.627244,2016: 623.324363,2032: 628.021349,2048: 632.718205}

"""


# Function C from AGPS20 source code
def caps_vol(d, theta, integrate=False, prec=None):
    """
    The probability that some v from the sphere has angle at most theta with some fixed u.

    :param d: We consider spheres of dimension `d-1`
    :param theta: angle in radians
    :param: compute via explicit integration
    :param: precision to use

    EXAMPLE::

        sage: C(80, pi/3)
        mpf('1.0042233739846629e-6')

    """
    prec = prec if prec else mp.prec
    with mp.workprec(prec):
        theta = mp.mpf(theta)
        d = mp.mpf(d)
        if integrate:
            r = (
                1
                / mp.sqrt(mp.pi)
                * mp.gamma(d / 2)
                / mp.gamma((d - 1) / 2)
                * mp.quad(lambda x: mp.sin(x) ** (d - 2), (0, theta), error=True)[0]
            )
        else:
            r = mp.betainc((d - 1) / 2, 1 / 2.0, x2=mp.sin(theta) ** 2, regularized=True) / 2
        return r


# Return log2 of the number of gates for FindAllPairs according to AGPS20
def agps20_gates(beta_prime):
    k = beta_prime / 8
    if k != round(k):
        x = k - floor(k)
        d1 = agps20_gates(8*floor(k))
        d2 = agps20_gates(8*(floor(k) + 1))
        return x * d2 + (1 - x) * d1
    return agps20_gate_data[beta_prime]



# Return log2 of the number of gates for FindAllPairs according to AGPS20
def agps20_gates(beta_prime):
    k = beta_prime / 16
    if k != round(k):
        x = k - floor(k)
        d1 = agps20_gates(16*floor(k))
        d2 = agps20_gates(16*(floor(k) + 1))
        return x * d2 + (1 - x) * d1
    return agps20_gate_data[beta_prime]


# Return log2 of the number of vectors for sieving according to AGPS20
def agps20_vectors(beta_prime):
    k = round(beta_prime)
    N = 1./caps_vol(beta_prime, mp.pi/3.)
    return log2(N)


# Progressivity Overhead Factor
C = 1./(1.- 2**(-.292))

# List decoding parameter, theoretical sieve cost = 2**(a*beta_prime+b)

list_decoding_classical = {"MATZOV22": (0.29613500308205365,20.387885985467914), "AGPS20": (0.2988026130564745, 26.011121212891872)}


#theo d4f 2
def dims4free(beta):
    return ceil(beta * ln(4./3.) / ln(beta/(2*pi*exp(1))))


#cost of bkz with progressive sieve
def theo_bkz_cost(n, beta,ldc_param = "AGPS20"):
    if(beta <=10):
        return (0,0)
    f = dim4free_wrapper(theo_dim4free_fun2,beta)
    beta_prime = floor(beta - f)
    
    if(ldc_param == "AGPS20"):
        gates = log2(1.*(n-beta)*C) + (list_decoding_classical[ldc_param][0] * beta_prime + list_decoding_classical[ldc_param][1]) #agps20_gates(beta_prime) #
    else:
        gates = log2(1.*(n-beta)*C) + (list_decoding_classical[ldc_param][0] * beta_prime + list_decoding_classical[ldc_param][1])
    bits = log2(8*beta_prime) + agps20_vectors(beta_prime)
    return (gates, bits)



def theo_pump_cost(beta,ldc_param = "AGPS20"):
    if(beta <=10):
        return (0,0)
    if(ldc_param == "AGPS20"):
        gates = log2(1.*C) + (list_decoding_classical[ldc_param][0] * beta + list_decoding_classical[ldc_param][1]) #agps20_gates(beta) #
    else:
        gates = log2(1.*C) + (list_decoding_classical[ldc_param][0] * beta + list_decoding_classical[ldc_param][1])

    bits = log2(8*beta) + agps20_vectors(beta)
    return (gates, bits)

    
#Return progressive sieve cost
def pump_cost(d,beta,cost_model = 1,ldc_param = "AGPS20"):
    if(cost_model == 1):
        return theo_pump_cost(beta,ldc_param = ldc_param)
    elif(cost_model == 2):
        return practical_pump_cost(beta)



def bkz_cost(d, beta,cost_model=1,ldc_param = "AGPS20"):
    if(cost_model == 1):
        return theo_bkz_cost(d, beta,ldc_param = ldc_param)
    elif(cost_model == 2):
        beta_prime = beta - dim4free_wrapper(default_dim4free_fun,beta)
        return log2(practical_bkz_cost(d,beta,1)),practical_pump_cost(beta_prime)[1]
    


###########################
#practical cost model


# threads = 32, gpus = 2,  pnj-bkz cost
def get_k1_k2_pnj(beta, sieve):
    if beta >=0 and beta <10 and not sieve:
        k1 = 0 
        k2 = 0
    elif beta>=10 and beta<=42 and not sieve:
        k1 = 0.03
        k2 = 5.188
    elif beta < 50 and not sieve:
        k1 = 0.19
        k2 = -1.741
    elif beta <= 97 and sieve:
        k1 = 0.056
        k2 = 7.85
    elif beta <= 118 and sieve:
        k1 = 0.215
        k2 = -7.61
    elif beta <= 128 and sieve:
        k1 = 0.314
        k2 = - 19.24
    else:
        k1 = 0.368
        k2 = -26.15
    return k1,k2


# threads = 32, gpus = 2, test pump
def get_k1_k2_pump(beta):
    if beta >=0 and beta <10:
        k1 = 0 
        k2 = 0
    elif beta>=10 and beta<=60:
        k1 = 0.035657
        k2 = -2.317327
    elif beta <= 96:
        k1 = 0.078794
        k2 = -0.039742
    elif beta <= 116:
        k1 = 0.231927
        k2 = -14.713430
    elif beta <= 128:
        k1 = 0.314
        k2 = -24.21
    else:
        k1 = 0.368
        k2 = -31.12
    # else: #30 > 80
    #     k1 = 0.3642 
    #     k2 = - 24.398 
    return k1,k2



#get pump time test in threads = 20
def practical_pump_cost(beta):
    #make sure not use the enum cost 
    #f = dim4free_wrapper(default_dim4free_fun,beta)
    #f = dim4free_wrapper(theo_dim4free_fun1,beta)
    #f = dim4free_wrapper(theo_dim4free_fun2,beta)
    #beta_prime = beta - f
    beta_prime = beta
    k1, k2 = get_k1_k2_pump(beta_prime) # threads = 20
    # k = (1/71.)*((1.33)**(beta/10.))
    secs = k1*beta_prime+k2


    #unit: GB
    if( beta_prime <= 56):
    	bits = 2.0311
    elif( beta_prime >=57 and  beta_prime <=63):
    	bits = 8e-5 *  pow(beta_prime,2) - 0.0083*beta_prime + 2.2555
    elif( beta_prime >= 64 and beta_prime <= 94):
    	bits = 2.195202e-6 * pow(beta_prime,4) - 6.297613e-4 * pow(beta_prime,3) +6.803540e-2 * pow(beta_prime,2) - 3.274476 * beta_prime + 61.31963
    elif( beta_prime >= 95):
    	bits = 2.0311 + pow(2,0.1992* beta_prime  - 18.714)
    
    bits = log2(bits * pow(2,33))

    return (secs, bits)  # n_expected = beta -f , beta = d-llb

    

#get pnj-BKZ time test in threads = 20
def practical_bkz_cost(d,beta,jump):
    if beta < 50:
        sieve = False
        f = 0
    else:
        sieve = True   
        f = dim4free_wrapper(default_dim4free_fun,beta)
    k1,k2 = get_k1_k2_pnj(beta-f,sieve)
    c3, c4 = 0.018, -2.24
    T_pnj = 2**(k1*(beta-f)+k2)
    #if(beta - f <= 60):
    #    pre_pnj_time = T_pnj/jump
    #else:
    pre_pnj_time = T_pnj*(c3*d+c4)/jump
    
    return round(pre_pnj_time,4)
 
