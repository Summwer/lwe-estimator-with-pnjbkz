import os
import requests
# from fpylll import IntegerMatrix

def load_lwe_challenge(n=40, alpha=0.005):
    """
    Load LWE challenge from file or website.

    :param n: LWE dimension
    :param alpha: the *standard deviation* of the secret is alpha*q

    """
    alpha = int(round(alpha * 1000))
    start = "lwechallenge"

    if not os.path.isdir(start):
        os.mkdir(start)

    end = "{n:03d}-{alpha:03d}-challenge.txt".format(n=n, alpha=alpha)
    #end = "{n:03d}-{alpha:03d}-midmat.txt".format(n=n, alpha=alpha)
    filename = os.path.join(start, end)
    if not os.path.isfile(filename):
        url = ("https://www.latticechallenge.org/lwe_challenge/challenges/"
               "LWE_{n:d}_{alpha:03d}.txt")
        url = url.format(n=n, alpha=alpha)
        r = requests.get(url)
        m = "Cannot retrieve challenge; server response was: '%s'. \n URL was: %s" % (r.reason, url)
        if not r.status_code == 200:
            raise ValueError(m)
        fn = open(filename, "w")
        fn.write(r.text)
        fn.close()

    data = open(filename, "r").readlines()
    n, m, q = [int(x) for x in [data[0], data[1], data[2]]]

    c_index = 3 if data[3].startswith("[") else 4

    A = eval(",".join([s_.replace(" ", ", ") for s_ in data[c_index+1:]]))
    # A = IntegerMatrix.from_matrix(A)
    A = matrix(A)
    c = tuple(eval(data[c_index].replace(" ", ", ")))
    return A, c, q

