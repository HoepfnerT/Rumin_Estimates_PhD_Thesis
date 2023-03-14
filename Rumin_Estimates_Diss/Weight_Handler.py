##### Start Weight_Handler.py #####
import sympy as sp
from sympy.solvers.solveset import linsolve
from Differential_Forms import Span_Differential_Forms
from Lie_Complex import Lie_Complex

class Weight_Handler():
    """For a given nilpotent Lie algebra G, 
    handle the computations and estimates concerning weights as used by Rumin."""

    def __init__(self, G):
        self.G = G
        self.initial_restraints = self.setup_weight_restraints()
        self.initial_homology_weights = self.setup_homology_weights()
    
    # Return a dict {i: w(X_i) for i = 1,..,G.dim} of possible weights on G satisfying w(X_i) + w(X_j) = w(X_k) if c_ij^k != 0.
    def setup_weight_restraints(self):
        weights     = [ sp.symbols("w" + str(i)) for i in range(self.G.dim) ]
        restraints  = [ weights[i] + weights[j] - weights[k] 
                        for k in self.G.structure_constants 
                        for i,j,c in self.G.structure_constants[k] 
                        if not c == 0 ] # setup linear system of equations based on structure constants
        if len(restraints) == 0:  return {sp.symbols("w" + str(i)): weights[i] for i in range(self.G.dim)}
        S = linsolve(restraints, weights[::-1])
        restrained_weights = list(list(S)[0])[::-1] # solution of linear system with minimal independent variables
        return {sp.symbols("w" + str(i)): wi for i, wi in enumerate(restrained_weights)}

    # Get the subspace of a linear space of a given pure weight.
    def get_pure_weight_part( self, span: Span_Differential_Forms, pure_weight ):
        out = []
        for form in span:
            weight = form.get_weights(self.initial_restraints)
            if len(weight) > 1: raise ValueError("Basisvector of mixed weights passed to Weight_Handler.get_pure_weight_part().")
            if list(weight)[0] == pure_weight: out.append(form)
        return Span_Differential_Forms(out)

    # Compute the weights appearing in the homology algebras of G under the initial weight restrictions
    def setup_homology_weights(self):
        HG = Lie_Complex(self.G)
        return {k: HG.get_homology_group_basis(k).get_weights(self.initial_restraints) for k in range(self.G.dim+1)}

    # Find restraints on weights to make a set of weights equal
    def find_pure_weight_restraits(self, weights_list: list):
        free_vars = set().union(*[sp.sympify(w).free_symbols for w in weights_list ])
        LGS = [weights_list[i] - weights_list[i+1] for i in range(len(weights_list)-1)]
        S = linsolve(LGS, list(free_vars)[::-1])
        if len(S) == 0: return self.initial_restraints
        else:           weight_restraints = {list(free_vars)[i]: list(S)[0][::-1][i] for i in range(len(free_vars))}
        if 0 in weight_restraints.values():     
            return {k: 0 for k,v in self.initial_restraints.items()}
        free_vars = set().union(*[expr.free_symbols for expr in weight_restraints.values() ])
        if len(free_vars) > 0:
            w = list(free_vars)[0]
            weight_restraints = {k: v.subs({w:1}) for k,v in weight_restraints.items()}
        return {k: v.subs(weight_restraints) for k,v in self.initial_restraints.items()}
##### End Weight_Handler.py #####