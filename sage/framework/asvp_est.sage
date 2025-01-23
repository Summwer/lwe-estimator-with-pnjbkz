load("../framework/utils.sage")
load("../framework/d_svp_prediction.sage")
load("../framework/cost.sage")
load("../framework/simulator/pnjbkz_simulator.sage")

def gamma_in_colattice_pump(log_rr, blocksizes):
    gh = gaussian_heuristic(log_rr)
    l = 0
    r = 0
    ap_norm = 0.
    for beta in blocksizes:
        r+= beta
        ap_norm += gaussian_heuristic(log_rr[l:r])
        l = r
    return sqrt(ap_norm/gh)

def gamma_in_bkz(log_rr, blocksize):
    gh = gaussian_heuristic(log_rr)
    return sqrt(gaussian_heuristic(log_rr[:blocksize])/gh)



#for lattice challenge: n = q

q = 122
dim = q*2
dvol = q*log(q)
print("Generate gs-lengths by GSA assumption.")
delta = compute_delta(2)
log_rr = [2*log(bkzgsa_gso_len(dvol, i, dim, delta=delta))  for i in range(dim)]
blocksizes = [65, 32, 55]
gcp = gamma_in_colattice_pump(log_rr, blocksizes)
gh = gaussian_heuristic(log_rr)
sv_norm = gcp*gh

print("reduced norm = %.1f" %sv_norm)

ccost = colattice_pump_cost(dim, blocksizes)
print("approximate factor in colattice-pump with ",blocksizes, end="")
print(" is %.3f" %gcp)
print("colattice-pump time cost: %.3f log2(gates), memory cost: %.3f log2(bits)" %(ccost[0],ccost[1]))

for blocksize in range(dim-1,10,-1):
    bkz_cost = theo_bkz_cost(dim, blocksize)
    if bkz_cost[0] <= ccost[0]:
        break
blocksize = 50
gcpbkz = gamma_in_bkz(log_rr, blocksize)
print("approximate factor in bkz-%d with " %blocksize , end="")
print(" is %.3f" %gcpbkz)
print("bkz time cost: %.3f log2(gates), memory cost: %.3f log2(bits)" %(bkz_cost[0],bkz_cost[1]))
