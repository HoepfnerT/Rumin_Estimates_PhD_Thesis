##### Start Lie_Complex.py #####
from Differential_Forms import Basis_Differential_Form, Differential_Form, Span_Differential_Forms, sort_sign
import itertools
import sympy as sp


class Lie_Complex():
    """
    Handle the Lie Complex and its homology for a given nilpotent Lie algebra G
    """
    def __init__(self, G):
        self.G = G
        self.chain_basis_list       = [ list(itertools.combinations(list(range(self.G.dim)), k)) for k in range(self.G.dim+1) ]
        self.differential_list      = self.setup_differential_list()
        self.homology_basis_list    = self.setup_homology_basis_list()
        self.envelope = [ direction for k in range(self.G.dim) for form in self.get_homology_group_basis(k) for direction in form.to_dict() ]


    # compute entries of matrix representing the differential
    def get_differential_entry(self, v: tuple, w: tuple):
        for k in v: 
            for i, j in itertools.combinations(w, 2):
                if set(list(w) + [k] ) == set(list(v) + [i,j]):
                    sgn = sort_sign(list(w) + [k]) * sort_sign(list(w) + [i,j])
                    return -sgn*self.G.get_structure_constant_ijk(i,j,k) 
        return 0

    def setup_differential_list(self):
        differential_list  = [ sp.Matrix([[0] for i in range(self.G.dim)]) ] # d^0: C^0(G) -> C^1(G)
        differential_list += [ sp.Matrix([[ self.get_differential_entry(v,w)    
                                              for v in self.chain_basis_list[m]   ]   
                                              for w in self.chain_basis_list[m+1] ]) 
                               for m in range(1, self.G.dim-1) ]             # d^m: C^m(G) -> C^(m+1)(G), n>m>0
        differential_list += [ sp.Matrix([[0 for i in range(self.G.dim)]]) ] # d^{n-1}: C^n(G) -> 0
        differential_list += [ sp.Matrix([[0]]) ]                            # d^n: C^n(G) -> 0
        return differential_list

    def setup_homology_basis_list(self):
        homology_groups = [ sp.Matrix([[1]]) ]  # H^0(G)
        for k in range(1,self.G.dim):                                                        
            rg = sp.Matrix(self.differential_list[k-1]).columnspace()    # orthogonal basis of image
            ker = sp.Matrix(self.differential_list[k]).nullspace()       # orthogonal basis of kernel 
            joined = rg + ker 
            H = joined[0]
            for m in joined[1:]: H = H.row_join(m)
            Q, _ = H.QRdecomposition()
            H = Q[:, len(rg):].T
            homology_groups.append(H)             # H^k(G)
        homology_groups.append(sp.Matrix([[1]]))  # H^n(G)
        return homology_groups
    
    def get_differential(self, k: int):
        if not k in range(self.G.dim+1): return 0
        return self.differential_list[k]

    # Return a list of span vectors of H^k(G) as a subspace in C^k(G).
    def get_homology_group_basis(self, k: int):
        if not k in range(self.G.dim+1): return 0
        H = self.homology_basis_list[k]
        basis = Span_Differential_Forms([
                    Differential_Form([ 
                        Basis_Differential_Form(c, bv) 
                        for c,bv in zip(H.row(j), self.chain_basis_list[k]) ])  
                    for j in range(H.shape[0]) ])
        basis.simplify()
        return basis

    def __str__(self): return f"Lie Complex of {self.G}."
    def __repr__(self):
        repr = f"For Lie Algebra {self.G}:\n" 
        for k in range(self.G.dim+1): 
            repr += f"\nDimension {k}:\nH^{k}(G) = {self.get_homology_group_basis(k)} \nd_{k} = {self.get_differential(k)}"
        return repr
##### End Lie_Complex.py #####
