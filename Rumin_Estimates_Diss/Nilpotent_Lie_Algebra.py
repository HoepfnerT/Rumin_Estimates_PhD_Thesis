##### Start Nilpotent_Lie_Algebra.py #####
class Nilpotent_Lie_Algebra():
    """Save the dimension and the structure constants of a nilpotent Lie Algebra."""
    def __init__(self, name, dim, structure_constants):
        self.name, self.dim, self.structure_constants = name, dim, structure_constants

    # Return c_ij^k
    def get_structure_constant_ijk(self, i,j,k):
        if not k in self.structure_constants:   return 0
        for i2,j2,c in self.structure_constants[k]: 
            if [i,j] == [i2, j2]:               return c
        return 0
    
    # Compute the lower central series 
    def get_lower_central_series(self):
        lcs = []
        g = list(range(self.dim)) # g0
        while len(g) > 0:
            lcs.append(g)
            g = [ k for k in g 
                  if (k in self.structure_constants)
                  and any([(i in g or j in g) and (not c == 0) 
                            for i,j,c in self.structure_constants[k]]) ] # g_n+1 <-- g_n
        return lcs

    # Read off the growth rate from the lower central series    
    def get_growth_rate(self):   return sum([len(g) for g in self.get_lower_central_series()])

    def __str__(self): return self.name
    def __repr__(self): return f"Nilpotent Lie Algebra {self.name} of dimension {self.dim} with structure constants {self.structure_constants}."
##### End Nilpotent_Lie_Algebra.py #####