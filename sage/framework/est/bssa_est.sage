#exits error

class BSSA():
    def __init__(self,verbose=False):
        self.BS = {}
        self.chisquared_table = {i: None for i in range(MAX_DIM)}

        for i in range(5000):
            self.chisquared_table[i] = RealDistribution('chisquared', i)

    #bssa: find minimal loop for each pnj-bkz-beta-jump
    def max_tour_for_pnjbkz_beta(self,l,cumulated_proba,beta,jump=1, dsvp_star = -1, progressive_sieve = False, cost_model = 1):
        l0 = deepcopy(l)
        cumulated_proba0 = deepcopy(cumulated_proba)
        remaining_proba0 = 1. - cumulated_proba0
        d = len(l0)
        (_,_,dsvp0,dsvp0_max) = d_svp_prediction(l0,beta,d,cumulated_proba0, progressive_sieve = progressive_sieve,cost_model=cost_model)
        
        
        loop = 0
        G1 = 0.

        l1= simulate_pnjBKZ(l0, beta, jump, 1)
        
        
        proba = 1.
        proba *= chisquared_table[beta].cum_distribution_function(2**(2 * l1[d-beta]))
        cumulated_proba1 = cumulated_proba0 + remaining_proba0 * proba
        remaining_proba1 = 1. - cumulated_proba1
        (_,_,dsvp1,dsvp1_max) = d_svp_prediction(l1,beta,d,cumulated_proba1,progressive_sieve = progressive_sieve,cost_model=cost_model)

        (G,_)= bkz_cost(d, beta, jump,cost_model=cost_model)
        G1 = log(2**G1 + ((2**G) * remaining_proba0 * proba))

        #print()
        #print(dsvp0, dsvp1,  dsvp_star, beta, jump, progressive_sieve)

        if dsvp1 <= dsvp_star:
            loop += 1
            return  l1, dsvp1, loop, G1, cumulated_proba1
        
        while dsvp0 - dsvp1 >= 1:
            loop += 1
            l0 = [_ for _ in l1]
            dsvp0 = dsvp1
            cumulated_proba0 = cumulated_proba1
            remaining_proba0 = remaining_proba1
   
            l1 = simulate_pnjBKZ(l0, beta, jump, 1)

            proba = 1.
            proba *= chisquared_table[beta].cum_distribution_function(2**(2 * l1[d-beta]))

            cumulated_proba1 = cumulated_proba0 + remaining_proba0 * proba
            remaining_proba1 = 1. - cumulated_proba1
            (_,_,dsvp1,dsvp1_max) = d_svp_prediction(l1,beta,d,cumulated_proba1,progressive_sieve=progressive_sieve,cost_model=cost_model)
            
            (G,_)= bkz_cost(d, beta, jump, cost_model)
            G1 = log(2**G1 + ((2**G) * remaining_proba0 * proba))

            if dsvp1 <= dsvp_star:
                loop += 1
                return  l1, dsvp1, loop, G1, cumulated_proba1

        G = float("inf")

        #print(beta,dsvp0)

        return  l0, dsvp0, loop, G, cumulated_proba0


    def BSSAGen(self,l0, beta_start, beta_goal, J=1, gap = 1, J_gap=1, progressive_sieve = False,cost_model = 1):
        d = len(l0)
        for beta in range(beta_start, beta_goal + 1,gap):
            Gmin = float("inf")
            for beta_sstart in range(beta_start + max(0,(floor((beta-beta_start)/gap)-10)*gap), beta,gap): #speed up the blocksize generation rate
                if beta_sstart == beta_start:
                    (S,l,G,cumulated_proba) = ([],l0,0.,0.)
                elif beta_sstart > beta_start:
                    if beta_sstart not in self.BS:
                        self.BSSAGen(l0, beta_start, beta_sstart,progressive_sieve=progressive_sieve,cost_model=cost_model)
                    (S,l,G,cumulated_proba) = self.BS[beta_sstart]
                
                G_tmp,beta_tmp = float("inf"),d
                #print()
                #print("===============")
                #print(beta_sstart,l)
                _, dsvp_star, _, _, _= self.max_tour_for_pnjbkz_beta(l,cumulated_proba,beta,progressive_sieve = progressive_sieve,cost_model=cost_model)
                #print(dsvp_star,beta)
                #print("--------------")
                #raise ""
                
                
                for beta_alg in range(beta+1, min(MAX_DIM,d)):
                    print("\r Blocksize strategy selection process: %3d --> %3d --> (%3d) --> %3d --> %3d" %(beta_start,beta_sstart,beta_alg,beta,beta_goal),end='')
                    #print()
                    #print(beta_tmp,beta_alg)
                    if (beta_alg >= beta_tmp + 3):
                        break
                    
                    for j in range(J,0,-J_gap):
                        if G_tmp == 0:
                            break;
                        l_, _, loop_, G_, cumulated_proba_= self.max_tour_for_pnjbkz_beta(l,cumulated_proba,beta,jump=j,dsvp_star = dsvp_star,progressive_sieve = progressive_sieve,cost_model=cost_model)
                        
                        if G_tmp > G_:
                            S_tmp, l_tmp, G_tmp, cumulated_proba_tmp,beta_tmp =[(beta_alg,j) for _ in range(loop_)], l_,  G_,cumulated_proba_, beta_alg
                        

                        #print("cost(bkz-%d,J=%d): %3.2f gate, current cost: %3.2f gate, min cost: %3.2f gate, min strategy: " %(beta_alg,J,bkz_cost(d, beta_alg, J, cost_model)[1],G_,G_tmp),end='')
                        #print(S_tmp)
                    
                    if G_tmp == 0:
                        break;

                G_tmp = log(2**G_tmp + 2**G)

                if(Gmin > G_tmp):
                    Gmin = G_tmp
                    self.BS[beta] = (S + S_tmp, l_tmp, G_tmp, cumulated_proba_tmp)

                    #(G_sieve,_,dsvptmp) = d_svp_prediction(l_tmp,beta,d,cumulated_proba_tmp,progressive_sieve = progressive_sieve,cost_model=cost_model)
                    #print("\n beta = %d, dsvp = %d, G1 = %3.2f, G_sieve = %3.2f, G = %3.2f, cum_prob = %3.2f" %(beta,dsvptmp,Gmin,G_sieve,log(2**Gmin + 2**G_sieve),cumulated_proba_tmp))
                    #if(cumulated_proba_tmp>=1.):
                    #    raise""
        

    #LWE estimation: progressive pnj-BKZs + Pump (BSSA_est)
    def BSSA_est(self,d, logvol, l, verbose=False, progressive_sieve = False,J=1,gap=1,J_gap=1,cost_model=1):
        """
        Computes the beta value for given dimension and volumes
        It is assumed that the instance has been normalized and sphericized, 
        i.e. that the covariance matrices of the secret is the identity
        :d: integer
        :vol: float
        """
        
        # Keep increasing beta to be sure to catch the second intersection
    
        # delta = compute_delta(2)
        # l = [log(bkzgsa_gso_len(logvol, i, d, delta=delta)) / log(2) for i in range(d)]
        J = ((J-1)//J_gap)*J_gap + 1
        self.BSSAGen(l,50,min(round(0.9*d),MAX_DIM),progressive_sieve=progressive_sieve,J=J,gap=gap,J_gap = J_gap,cost_model=cost_model)
        
        print()
        Gmin,Bmin,dsvpmin,dsvpmaxmin, Smin = float("inf"), float("inf"), d,d,[]
        G1min = float("inf")
        for beta in self.BS:
            (S,l,G1,cumulated_proba) = self.BS[beta]
            
            (G_sieve,_,dsvp,dsvpmax) = d_svp_prediction(l,beta,d,cumulated_proba, progressive_sieve = progressive_sieve,cost_model=cost_model)
            
            b = max([_[0] for _ in S]+[dsvp])

            G = log(2**G1 + 2**G_sieve)
            _,B = pump_cost(d,b,cost_model=cost_model)
            if(G!= float("inf") and B!= float("inf")):
                if(G < Gmin and B < Bmin):
                    Gmin, Bmin, dsvpmin, dsvpmaxmin, Smin = G, B, dsvp,dsvpmax, S
                    G1min = G1
        print("Min Strategy generated by BSSA is", Smin)
        print("dsvp = %d" %dsvpmin)
        return Smin, dsvpmin, dsvpmaxmin, Gmin, G1min, Bmin, None
