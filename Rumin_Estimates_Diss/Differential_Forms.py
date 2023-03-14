##### Start Differential_Forms.py #####
import sympy as sp
# Sort list, return sign of permutation needed to sort
def sort_sign(w: tuple):
    swap_counter = 0
    for i in range(len(w)-1):
        for j in range(i,len(w)):
            if w[i] > w[j]: swap_counter += 1
    if swap_counter % 2 == 0: return 1
    return -1

class Basis_Differential_Form():
    """Describe a differential form of type   coefficient * dx_I   by its coefficient and the ordered tuple I."""
    def __init__(self, coefficient, direction: tuple):
        self.coefficient = sort_sign(direction)*sp.sympify(coefficient)
        self.direction   = tuple(sorted(direction))

    def simplify(self):
        try:    self.coefficient = sp.simplify(self.coefficient)
        except: pass
        return  self

    # Compute weight based on a weight directionary for 1-forms
    def get_weight(self, weights_dict: dict):
        if self == 0: raise ValueError("The 0-form has no well-defined weight.")
        try:    return sum([weights_dict[sp.symbols("w" + str(i))] for i in self.direction])
        except KeyError: 
            raise KeyError(f"{self.direction} not contained in the weights dictionary {weights}")

    ## Arithmetics
    def __add__(self, other):
        if self.direction == other.direction: 
            return Basis_Differential_Form(self.coefficient + other.coefficient, self.direction)
        return NotImplemented
    def __neg__(self):          return Basis_Differential_Form(-self.coefficient, self.direction)
    def __sub__(self, other):   return self + (-other)
    def __bool__(self):         return not(self.coefficient == 0)
    def __mul__(self, other):  
        if isinstance(other, (int, float)):
            return Basis_Differential_Form(other*self.coefficient, self.direction)
        if isinstance(other, Basis_Differential_Form):
            if not len(set(self.direction + other.direction)) == len(self.direction) + len(other.direction): 
                return 0
            return Basis_Differential_Form( self.coefficient * other.coefficient * sort_sign(self.direction + other.direction), 
                                            tuple(sorted(self.direction + other.direction)) )
    def __str__(self):          return f"({self.coefficient}) {self.direction}"
    def __repr__(self):         
        if len(self.direction) == 0:    return f"{self.coefficient}"
        if self.coefficient == 1:       return "\\vartheta_{" + "".join([str(i) for i in self.direction]) + "}"
        return f"({self.coefficient})" + "\\vartheta_{" + "".join([str(i) for i in self.direction]) + "}"

class Differential_Form():
    """Describe a differential form given as a sum of basis differential forms by the list of these basis forms"""
    def __init__(self, summands: list):
        self.summands = summands
    # Simplify a differential form by adding basis forms of same directon
    # Set flag rescale to true if rescaling is allowed
    def simplify(self, rescale = False):
        self_dict = {}
        for summand in self.summands:
            if summand.direction in self_dict:
                self_dict[summand.direction].append(summand.coefficient)
            else:
                self_dict[summand.direction] = [summand.coefficient]
        self_dict = {direction: sum( self_dict[direction] ) for direction in self_dict }
        self.summands = [ Basis_Differential_Form(self_dict[direction], direction).simplify()
                          for direction in self_dict ]
        self.summands = [ summand for summand in self if summand ] # remove 0 entries
        if rescale: 
            min_coeff = min([ abs(summand.coefficient) for summand in self ])
            for summand in self.summands: summand.coefficient /= min_coeff
        self.summands.sort(key = lambda x: x.direction)
        return self

    def to_dict(self):          return {summand.direction: summand.coefficient for summand in self}

    def get_weights(self, weights_dict: dict):
        return list(set([summand.get_weight(weights_dict) for summand in self]))

    ## Arithmetics
    def __add__(self, other):   return Differential_Form(self.summands + other.summands).simplify()
    def __neg__(self):          return Differential_Form([-summand for summand in self])
    def __sub__(self, other):   return self + (-other)
    def __bool__(self):         return not(len(self.summands) == 0)

    def __mul__(self, other): 
        if isinstance(other, (int, float)):
            if other == 0: return Differential_Form([Basis_Differential_Form(0, ())])
            return Differential_Form([bf1*other for bf1 in self.summands])
        elif isinstance(other, Differential_Form):
            return Differential_Form([bf1*bf2 for bf1 in self.summands for bf2 in other.summands if bf1*bf2]).simplify()

    def __iter__(self):         yield from self.summands

    def __str__(self):         
        if len(self.summands) == 0: return "0"
        return "  +  ".join([str(summand) for summand in self])
    def __repr__(self):
        if len(self.summands) == 0: return "0"
        return "  +  ".join([repr(summand) for summand in self])


class Span_Differential_Forms():
    """ Describe a linear subspace of differential forms given by its basis vectors by a list of these basis vectors """
    def __init__(self, basis: list):
        self.basis = basis
    
    def simplify(self):
        for form in self.basis: form.simplify(rescale = True)
        self.basis = [ form for form in self.basis if form ]

    def get_weights(self, weights_dict: dict):
        return list(set().union(*[form.get_weights(weights_dict) for form in self]))
        
    def __iter__(self):         yield from self.basis
    def __str__(self):          return "<   " + ",\n\t".join([str(form) for form in self.basis]) + "   >"
    def __repr__(self):         return "\\langle " + ",\n\t\t".join([repr(form) for form in self.basis]) + "  \\rangle"
##### End Differential_Forms.py #####