#exits error..

class EnumBS():
    def __init__(self,verbose=False):
        self.BS = []
        self.chisquared_table = {i: None for i in range(MAX_DIM)}

        for i in range(5000):
            self.chisquared_table[i] = RealDistribution('chisquared', i)
    
    # Input: A list whose element is from large to small
    # Return the first index of the largest number <= x in high dimensional list
    # If it is non-existent, then return -1
    def binarySearch(self, arr, x): 
        # Basic Determination
        l, r = 0, len(arr)-1
        while r > l: 
            mid = (l + r) // 2
    
            # If x < mid element, then compare it to arr[l:mid]
            if arr[mid] <= x: 
                r = mid
    
            # If x >= mid element, then compare it to arr[mid+1:r+1]
            else: 
                l = mid + 1
        if(arr[l] > x):
            # Non-existent
            return -1    
        return l
    
    
    def find_pos(self,arr,x):
        if arr[0] < x:
            return  -1 #all elements < x
        for pos in range(len(arr)):
            if(arr[pos]<x):
                return pos-1
        return len(arr) #all elements in arr >=x


    def BS_add(self,bs_):
        (S_,avgG_,dsvp_,dsvpmax_,l_,cumulated_proba_) = bs_
        #print(self.BS[0][:3])
        #print([(S_,bs[2],bs[1]) for bs in self.BS])
        #print(dsvp_)
        
        #pos = self.binarySearch([bs[2] for bs in self.BS], dsvp_)
        pos = self.find_pos([bs[2] for bs in self.BS], dsvp_) #Return the minimal number >=dsvp_
        
       
        #print(pos,dsvp_,[bs[2] for bs in self.BS][pos])
        #print([bs[2] for bs in self.BS][pos],dsvp_)
        #print(len(self.BS),pos,dsvp_)
        try:
            if(pos != len(self.BS) and pos != -1):
                #print(pos,dsvp_,[bs[2] for bs in self.BS][pos])
                assert([bs[2] for bs in self.BS][pos]>=dsvp_)
        except AssertionError:
            print()
            print("Find an error position for dsvp_!")
            print("pos=",pos,[bs[2] for bs in self.BS][pos],"<",dsvp_,"!")
            raise ""
        
        #print()
        #print("index=%4d, dsvp=%4d, dsvp_=%.4e" %(pos,self.BS[pos][2],dsvp_))
        

        
        #If dsvp <= dsvp_ and avgG <= avgG_, then cannot add the target strategy
        if pos+1 < len(self.BS): #It exists a dsvp s.t. dsvp <= dsvp_
            avgG, dsvp = self.BS[pos+1][1], self.BS[pos+1][2]
            if avgG_ >= avgG:
                return

        #print()
        #print("pos = %d, BS[pos] = %d,  dsvp_ = %d, avgG_ = %3.2f" %(pos,[bs[2] for bs in self.BS][pos],dsvp_,avgG_))
        #print([(bs[2],bs[1]) for bs in self.BS])

        
        #Otherwise, we should remove all bad strategy and add the target strategy
        if(pos == len(self.BS)):
            pos -= 1
        avgG, dsvp = self.BS[pos][1], self.BS[pos][2]

        
        #print(pos,len(self.BS))
        #all elements in dsvp >= dsvp_
        if(dsvp_ == dsvp and avgG_ >= avgG):
            return
        
        while((dsvp_ <= dsvp and  avgG_ <= avgG) and pos >= 0):
            #print(pos,dsvp_,dsvp,avgG_,avgG)
            #print("=================")
            self.BS.remove(self.BS[pos])
            pos -=1
            if(pos==-1 or self.BS==[]):
                break
            avgG, dsvp = self.BS[pos][1], self.BS[pos][2]
            
        self.BS.insert(pos+1,(S_,avgG_,dsvp_,dsvpmax_,l_,cumulated_proba_))
        

        
        #print(pos,dsvp_)
        #print([_[2] for _ in self.BS] )
        try:
            assert([_[2] for _ in self.BS] ==  sorted([_[2] for _ in self.BS],reverse =True)) #Verify order of dvsp (large --> small)
            assert([_[1] for _ in self.BS] ==  sorted([_[1] for _ in self.BS])) #Verify order of Gbkz( small --> large)
        except AssertionError:
            print()
            print("Insert in error!")
            print([_[2] for _ in self.BS])
            print(pos,dsvp_)
            raise ""

        #print([bs[2] for bs in self.BS][:10])
        
        #List = [bs[2] for bs in self.BS]
        #if(List!=sorted(List,reverse=True)):
        #    raise ("error")

        
        
    def max_tour_for_pnjbkz_beta(self,bs,beta,jump,progressive_sieve = False,cost_model=1):
        (S0,avgG,dsvp0,dsvp0max,l0,cumulated_proba0) = deepcopy(bs)
        remaining_proba0 = 1. - cumulated_proba0

        d = len(l0)
        loop = 0

        #print(dsvp0,dsvp0max,S0)
    
        l1= simulate_pnjBKZ(l0, beta, jump, 1)
        proba = 1.
        proba *= self.chisquared_table[beta].cum_distribution_function(
                2**(2 * l1[d-beta]))
        
        #avgb += beta * remaining_proba0 * proba

        G1, _ = bkz_cost(d,beta,jump,cost_model=cost_model)
        avgG = log(2**avgG + ((2**G1) * remaining_proba0 * proba))

        cumulated_proba1 = cumulated_proba0 + remaining_proba0 * proba
        remaining_proba1 = 1. - cumulated_proba1
        (_,_,dsvp1,dsvp1max) = d_svp_prediction(l1,beta,d,cumulated_proba1,progressive_sieve = progressive_sieve,cost_model=cost_model)
        
        while dsvp0 - dsvp1 >= 1:
            loop += 1
            l0 = [_ for _ in l1]
            dsvp0 = dsvp1
            dsvp0max = dsvp1max
            cumulated_proba0 = cumulated_proba1
            remaining_proba0 = remaining_proba1
            S0.append((beta,jump))

            
            bs0 = (S0,avgG,ceil(dsvp0),dsvp0max,l0,cumulated_proba0)
            self.BS_add(bs0)
        
            l1 = simulate_pnjBKZ(l0, beta, jump, 1)
    

            #print("cost of BS_add: %3.2f s, cost of simulator: %3.2f s"%(T_BS_add,T_l1))

            proba = 1.
            proba *= self.chisquared_table[beta].cum_distribution_function(
                    2**(2 * l1[d-beta]))
            #avgb += beta * remaining_proba0 * proba

            G1, _ = bkz_cost(d,beta,jump,cost_model=cost_model)
            avgG = log(2**avgG + ((2**G1) * remaining_proba0 * proba))

            #print(avgb,beta,remaining_proba0,proba)

            cumulated_proba1 = cumulated_proba0 + remaining_proba0 * proba
            remaining_proba1 = remaining_proba0 * (1. - proba)
            (_,_,dsvp1,dsvp1_max) = d_svp_prediction(l1,beta,d,cumulated_proba1,progressive_sieve = progressive_sieve,cost_model=cost_model)
        

    #LWE estimation: progressive pnj-BKZs + Pump (EnumBS)
    def EnumBS_estimation(self,d, logvol, l, verbose=False, progressive_sieve = False, J = 1, gap=1, J_gap = 1, cost_model = 1):
        """
        Computes the beta value for given dimension and volumes
        It is assumed that the instance has been normalized and sphericized, 
        i.e. that the covariance matrices of the secret is the identity
        :d: integer
        :vol: float
        """
        
        Gmin, Bmin, avgGmin, dsvpmin = float("inf"), float("inf"), d, d

        # Keep increasing beta to be sure to catch the second intersection
        
        remaining_proba = 1.
        average_beta = 0.
        cumulated_proba = 0.

        #delta = compute_delta(2)
        #l = [log(bkzgsa_gso_len(logvol, i, d, delta=delta)) / log(2) for i in range(d)]
        (_,_,dsvp,dsvpmax) = d_svp_prediction(l,0,d,cumulated_proba,progressive_sieve = progressive_sieve,cost_model=cost_model)
        self.BS = [([],0,ceil(dsvp),dsvpmax,l,cumulated_proba)]

        k = 0
        
        beta_start = 50
        J = ((J-1)//J_gap)*J_gap + 1
        while( k < len(self.BS)):
            bs = self.BS[k]
            S = bs[0]
            if S ==[]:
                beta_start = 50
            else:
                beta_start = max([_[0] for _ in S])
            

            for beta in range(beta_start+1,min(MAX_DIM,d),gap):
                if verbose:
                    print('\r index: %4d, beta: %4d, goal index: %4d' %(k+1,beta,len(self.BS)),end='')
                for j in range(J,0,-J_gap):
                    self.max_tour_for_pnjbkz_beta(bs,beta,j, progressive_sieve = progressive_sieve,cost_model=cost_model)
                #print(beta_start,k,len(self.BS),S)
                

            k += 1
        print()
        
            
        
        for bs in self.BS:
            (S,avgG,dsvp,dsvpmax,l,cumulated_proba) = bs

            #G1,B1 = bkz_cost(d,average_beta,cost_model=cost_model)
            G_sieve,B_sieve = pump_cost(d,dsvp,cost_model=cost_model)

            G = log(2**avgG + 2**G_sieve)
            if S== []:
                B = B_sieve
            else:
                B = max(pump_cost(d,max([_[0] for _ in S]),cost_model=cost_model)[1], B_sieve)
            if(G!= float("inf") and B!= float("inf")):
                if(G < Gmin and B < Bmin):
                    Gmin, Bmin, G1min, dsvpmin, dsvpmaxmin, Smin = G, B, avgG, dsvp, dsvpmax, S

            #if verbose:
                #os.system('clear')
                #print("β= %d,  dsvp= %3.2f,  pr=%.2e,  cum-pr=%.2e,  rem-pr=%.2e,  β_min=%d,  dsvp_min=%3d,  G=%3.2f gate,  B=%3.2f bit"%
                #        (beta, dsvp, proba, cumulated_proba, remaining_proba, avgbetamin, dsvpmin, Gmin, Bmin), 
                #        end="\r" if cumulated_proba < 1e-4 else "\n")
                    #sys.stdout.flush()
        print(S)
        return S, G1min, dsvpmin, dsvpmaxmin, Gmin, Bmin, None

    @parallel(multiprocessing.cpu_count())
    def max_tour_for_pnjbkz_beta_in_parallel(self,bs,beta,jump,progressive_sieve = False,cost_model=1):
        #print("\r beta = %4d" %(beta), end='')
        tmpBS = []
        
        (S0,avgG,dsvp0,dsvp0max,l0,cumulated_proba0) = deepcopy(bs)
        remaining_proba0 = 1. - cumulated_proba0

        d = len(l0)
        loop = 0
    
        l1= simulate_pnjBKZ(l0, beta, jump, 1)
        proba = 1.
        proba *= chisquared_table[beta].cum_distribution_function(2**(2 * l1[d-beta]))
        #avgb += beta * remaining_proba0 * proba
        G1, _ = bkz_cost(d,beta,jump,cost_model=cost_model)
        avgG = log(2**avgG + ((2**G1) * remaining_proba0 * proba))

        
        cumulated_proba1 = cumulated_proba0 + remaining_proba0 * proba
        remaining_proba1 = remaining_proba0 * (1. - proba)
        (_,_,dsvp1,dsvp1max) = d_svp_prediction(l1,beta,d,cumulated_proba1,progressive_sieve = progressive_sieve,cost_model=cost_model)
        
        while dsvp0 - dsvp1 >= 1:
            loop += 1
            l0 = [_ for _ in l1]
            dsvp0 = dsvp1
            dsvp0max = dsvp1max
            cumulated_proba0 = cumulated_proba1
            remaining_proba0 = remaining_proba1
            S0.append((beta,jump))

            bs0=(S0,avgG,ceil(dsvp0),dsvp0max,l0,cumulated_proba0)

            T0_BS_add = time.time()
            tmpBS.append(bs0)
            T_BS_add = time.time()-T0_BS_add

            T0_l1 = time.time()
            l1 = simulate_pnjBKZ(l0, beta, jump, 1)
            T_l1 = time.time() - T0_l1

            #print("cost of BS_add: %3.2f s, cost of simulator: %3.2f s"%(T_BS_add,T_l1))

            proba = 1.
            proba *= chisquared_table[beta].cum_distribution_function(
                    2**(2 * l1[d-beta]))
            #avgb += beta * remaining_proba0 * proba

            G1, _ = bkz_cost(d,beta,jump,cost_model=cost_model)
            avgG = log(2**avgG + ((2**G1) * remaining_proba0 * proba))

            #print(avgb,beta,remaining_proba0,proba)

            cumulated_proba1 = cumulated_proba0 + remaining_proba0 * proba
            remaining_proba1 = remaining_proba0 * (1. - proba)
            (_,_,dsvp1,dsvp1max) = d_svp_prediction(l1,beta,d,cumulated_proba1,progressive_sieve = progressive_sieve,cost_model=cost_model)
        
        return tmpBS


    #LWE estimation: progressive pnj-BKZs + Pump (EnumBS)
    def EnumBS_estimation_in_parallel(self,d, logvol, l, verbose=False, progressive_sieve = False, J = 1, gap=1, J_gap = 1, cost_model = 1):
        """
        Computes the beta value for given dimension and volumes
        It is assumed that the instance has been normalized and sphericized, 
        i.e. that the covariance matrices of the secret is the identity
        :d: integer
        :vol: float
        """
        
        Gmin, Bmin, avgGmin, dsvpmin = float("inf"), float("inf"), d, d

        # Keep increasing beta to be sure to catch the second intersection
        
        remaining_proba = 1.
        average_beta = 0.
        cumulated_proba = 0.

        #delta = compute_delta(2)
        #l = [log(bkzgsa_gso_len(logvol, i, d, delta=delta)) / log(2) for i in range(d)]
        (_,_,dsvp,dsvpmax) = d_svp_prediction(l,0,d,cumulated_proba,progressive_sieve = progressive_sieve,cost_model=cost_model)
        self.BS = [([],0,ceil(dsvp),dsvpmax,l,cumulated_proba)]

        k = 0
        J = ((J-1)//J_gap)*J_gap + 1
               
        while( k < len(self.BS)):
            bs = self.BS[k]
            S = bs[0]
            if S ==[]:
                beta_start = 50
            else:
                beta_start = max([_[0] for _ in S])
            

            T0_sim = time.time()
            
           
            
            result = self.max_tour_for_pnjbkz_beta_in_parallel([(bs,beta,j,progressive_sieve,cost_model) for beta in range(beta_start+1,min(MAX_DIM,d), gap) for j in range(J,0,-J_gap)])

            
            tmpBS = []
            for _ in result:
                tmpBS.extend(_[1])


            
            #for beta in range(beta_start+1,min(100,d)):
            #    if verbose:
            #        print('\r index: %4d, beta: %4d, goal index: %4d' %(k+1,beta,len(self.BS)),end='')
            #    self.max_tour_for_pnjbkz_beta_in_parallel(bs,beta,1,cost_model=cost_model)
                #print(beta_start,k,len(self.BS),S)
            T_sim = time.time() - T0_sim
            T0_bs = time.time()
            for bs in tmpBS:
                self.BS_add(bs)

            T_bs = time.time() - T0_bs
            

            if verbose:
                print('\r index: %4d, beta from %4d to %4d, goal index: %10d, cost for simulation: %3.2f s, cost for BS_add: %3.2f s' %(k+1,beta_start+1,min(MAX_DIM,d),len(self.BS),T_sim,T_bs),end='')
            
            k += 1
        print()
        
            
        
        for bs in self.BS:
            (S,avgG,dsvp, dsvpmax, l,cumulated_proba) = bs

            #G1,B1 = bkz_cost(d,average_beta,cost_model=cost_model)
            G_sieve,B_sieve = pump_cost(d,dsvp,cost_model=cost_model)

            G = log(2**avgG + 2**G_sieve)
            if S== []:
                B = B_sieve
            else:
                B = max(pump_cost(d,max([_[0] for _ in S]),cost_model=cost_model)[1], B_sieve)
            if(G!= float("inf") and B!= float("inf")):
                if(G < Gmin and B < Bmin):
                    Gmin, Bmin, G1min,  dsvpmin, dsvpmaxmin, Smin = G, B, avgG, dsvp, dsvpmax, S

            #if verbose:
                #os.system('clear')
                #print("β= %d,  dsvp= %3.2f,  pr=%.2e,  cum-pr=%.2e,  rem-pr=%.2e,  β_min=%d,  dsvp_min=%3d,  G=%3.2f gate,  B=%3.2f bit"%
                #        (beta, dsvp, proba, cumulated_proba, remaining_proba, avgbetamin, dsvpmin, Gmin, Bmin), 
                #        end="\r" if cumulated_proba < 1e-4 else "\n")
                    #sys.stdout.flush()
        print("Min Strategy generated by EnumBS is", Smin)
        print("dsvp = %d" %dsvpmin)
        return Smin, G1min, dsvpmin, dsvpmaxmin, Gmin, Bmin, None
