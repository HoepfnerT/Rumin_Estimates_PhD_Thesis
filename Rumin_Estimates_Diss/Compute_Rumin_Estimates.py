##### Start Compute_Rumin_Estimates.py #####
import sympy as sp
from Lie_Complex import Lie_Complex
from Weight_Handler import Weight_Handler
# Compute the estimates on the Novikov-Shubin invariants accordingly to Rumin's method; Applying Hodge duality if flag is set True
# Returns a dict {k: (low_k, high_k)} for such k where estimate   low_k <= alpha_k <= high_k   follows from Rumin's approach
def Compute_Rumin_Estimates(G, hodge_duality = True):
    WEIGHT_HANDLER = Weight_Handler(G)
    NSinvars_estimates = {0: (G.get_growth_rate(), G.get_growth_rate()) }
    HG = Lie_Complex(G)
    for k in range(1, G.dim):
        k_weights = HG.get_homology_group_basis(k).get_weights(WEIGHT_HANDLER.initial_restraints)
        pure_weight_restraints = WEIGHT_HANDLER.find_pure_weight_restraits(k_weights)
        if 0 in pure_weight_restraints.values(): continue
        k_weights   = set([ w.subs(pure_weight_restraints) for w in WEIGHT_HANDLER.initial_homology_weights[k] ])
        k1_weights  = set([ w.subs(pure_weight_restraints) for w in WEIGHT_HANDLER.initial_homology_weights[k+1] ])
        N           = sp.sympify(WEIGHT_HANDLER.initial_homology_weights[G.dim][0].subs(pure_weight_restraints))
        diffs       = [ sp.sympify(w-v) for w in k1_weights for v in k_weights ]
        # check if still variables remain, replace by 1 
        free_vars = set().union(*[list(v.free_symbols) for v in diffs])
        if len(free_vars) > 0: 
            replace_dict = {fv: 1 for fv in free_vars}
            diffs = [v.subs(replace_dict) for v in diffs]
            N = N.subs(replace_dict)
            pure_weight_restraints = {k: v.subs(replace_dict) for k,v in pure_weight_restraints.items()}
        Nmax, Nmin  = max(diffs), min(diffs)
        if not Nmin > 0: Nmin = min([1/sp.sympify(w) for w in pure_weight_restraints.values()])
        NSinvars_estimates[k] = (N / Nmax, N / Nmin)
    if hodge_duality:
        for k in range(G.dim//2):
            if (k in NSinvars_estimates) and not (G.dim-k-1 in NSinvars_estimates):    
                NSinvars_estimates[G.dim-k-1] = NSinvars_estimates[k]
            elif not (k in NSinvars_estimates) and (G.dim-k-1 in NSinvars_estimates):
                NSinvars_estimates[k] = NSinvars_estimates[G.dim-k-1]
            elif (k in NSinvars_estimates) and (G.dim-k-1 in NSinvars_estimates):
                n_min = max(NSinvars_estimates[k][0], NSinvars_estimates[G.dim-k-1][0])
                n_max = min(NSinvars_estimates[k][1], NSinvars_estimates[G.dim-k-1][1])
                NSinvars_estimates[k] = (n_min, n_max)
                NSinvars_estimates[G.dim-k-1] = (n_min, n_max)
    return NSinvars_estimates
##### End Compute_Rumin_Estimates.py #####