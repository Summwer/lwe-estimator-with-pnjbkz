"""
Estimator for LWR and uniform error LWE by random sampling

min_beta: beta in the first step of the two-step mode
max_beta: maximal sieving dimension the the second step of the two-step mode
sample_space, sample_vector: number of subspaces and vectors in the random sampling
type: 'uniform': both secret and error are uniform; 'gaussian': both secret and error are gaussian;
         'mixed': secret is gaussian and error is uniform

"""

load("leaky-LWE-Estimator-master/framework/utils.sage")

def estimate_LWR_instance(n,m,q,p,sigma_s,min_beta,max_beta,sample_space=10,sample_vector=1000):

    logvol=m*log(q)-n*log(sigma_s)-m*log((q*q-gcd(p,q)*gcd(p,q))/(p*p*12))
    compute_uniform_twostep(n,m,logvol,min_beta,max_beta,sample_space,sample_vector)

def estimate_uniform_LWE_instance(n,m,q,sigma_e,sigma_s,min_beta,max_beta,sample_space=10,sample_vector=1000):

    logvol=m*log(q)-n*log(sigma_s)-m*log(sigma_e)
    compute_uniform_twostep(n,m,logvol,min_beta,max_beta,sample_space,sample_vector,type='uniform')

def compute_uniform_twostep(n,m,logvol,min_beta,max_beta,sample_space,sample_vector,type='mixed'):

    d=n+m+1

    gsap=[0.0]*max_beta
    gsap1=[0.0]*max_beta
    cum=[0]*max_beta
    cum1=[0]*max_beta
    
    delta = compute_delta(2)
    l = [log(bkzgsa_gso_len(logvol, i, d, delta=delta)) / log(2) for i in range(d)]
    for beta in [x for x in range(2, max_beta)]:
        l = simBKZ(l, beta, 1)
        delta = compute_delta(beta)
        i = d - beta
        gsap[beta] = 2**(2 * l[i])
        if(beta==min_beta):
            for beta0 in range(2, max_beta):
                i0 = d - beta0
                gsap1[beta0] = 2**(2 * l[i0])

    dist = RealDistribution('gaussian',1)
    X=RealDistribution('uniform',[0.0-sqrt(3),sqrt(3)])
    aa=[0.0]*d

    #print("------------")

    for n in range(sample_space):
        #print(n)
        subs = []
        for i in range(max_beta):
            subs.append([dist.get_random_element() for _ in range(d)])
        for i in range(max_beta):
            for j in range(i):
                a=0.0
                b=0.0
                for k in range(d):
                    a=a+subs[i][k]*subs[j][k]
                    b=b+subs[j][k]*subs[j][k]
                for k in range(d):
                    subs[i][k]=subs[i][k]-subs[j][k]*a/b
            c=0.0
            for k in range(d):
                c=c+subs[i][k]*subs[i][k]
            for k in range(d):
                subs[i][k]=subs[i][k]/sqrt(c)

        for j in range(sample_vector):
            #if((j%100)==0):
            #    print(j)
            if(type=='gaussian'):
                for i in range(m):
                    aa[i]=dist.get_random_element()
            else:
                for i in range(m):
                    aa[i]=X.get_random_element()
            if(type=='uniform'):
                for i in range(m,m+n):
                    aa[i]=X.get_random_element()
                for i in range(m,m+n):
                    aa[i]=dist.get_random_element()
            aa[d-1]=1
            len=0.0
            for k in range(max_beta-1):
                proj=0.0
                for i in range(d):
                    proj = proj + subs[k][i]*aa[i]
                len=len+proj*proj
                if(len<gsap[k+1]):
                    cum[k+1]=cum[k+1]+1
                if(len<gsap1[k+1]):
                    cum1[k+1]=cum1[k+1]+1

    average_beta=0.
    cumulated_proba=0.
    remaining_proba=1.
    for beta in range(2,max_beta):
        average_beta += beta * remaining_proba * cum[beta] / sample_vector / sample_space
        cumulated_proba += remaining_proba * cum[beta] / sample_vector / sample_space
        remaining_proba *= 1. - cum[beta] / sample_vector / sample_space

    print("leaky-LWR-estimator: average-beta:")
    print(average_beta)
    print(remaining_proba)
    print("---------------")

    cumulated_time=0.
    cumulated_proba=0.
    remaining_proba=1.
    for beta in range(2,max_beta):
        cumulated_time += (1.5 ** (0.5*beta)) * (cum1[beta] / sample_vector / sample_space - cumulated_proba)
        cumulated_proba = cum1[beta] / sample_vector / sample_space
        remaining_proba = 1. - cum1[beta] / sample_vector / sample_space

    print("two-step-LWR-estimator: average-d_svp:")
    print(log(cumulated_time)/log(sqrt(1.5)))
    print(remaining_proba)