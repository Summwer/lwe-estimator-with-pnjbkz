#include "utils.h"

void initialize_from_LWE_instance(ZZ_mat<ZT> &A, ZZ_mat<ZT> &B, ZZ_mat<FT> &S,vector<Z_NR<ZT>> &b, vector<Z_NR<ZT>> &u, vector<double> &mu, int n, int q, int m, map<int,rational<int>> D_e, map<int,rational<int>> D_s, bool verbosity=true);