##### Start Print_Methods.py #####
from NLG_List import *
from Nilpotent_Lie_Algebra import Nilpotent_Lie_Algebra as NLG
from Compute_Rumin_Estimates import Compute_Rumin_Estimates
# Given a list of graded nilpotent Lie algebras and the maximal dimension 
# returns a table containing the NS-estimates from Rumin's method.
def make_NS_invars_table( NLG_list: list, MAX_DIM: int, COL_SEP = 14, hodge_duality = True, Latex = False ):
    outstr = ""
    if Latex:
        out_str = "Estimates on the Novikov-Shubin invariants of the deRham differentials.\n"
        out_str += "\\begin{tabular}{|c|" + "|c"*(MAX_DIM+1) + "|} \\hline \n"
        out_str += " $\\mathfrak{g}$ " + "".join([f" & $\\alpha_{k}$" for k in range(MAX_DIM+1)]) + " \\\\ \\hline \n"
    else:
        out_str = " " + "_"*((MAX_DIM+1)*(COL_SEP+3)-1) + "\n" + "|"
        out_str +=  "_"*COL_SEP + "__|_{}_|".format("_|_".join("_"*((COL_SEP-len(str(k))-2)//2) +"a_" 
                                                               + str(k) + "_"*((COL_SEP-len(str(k))-1)//2) for k in range(MAX_DIM)))
        out_str += "\n"
    for LG in NLG_list:
        G = NLG(*LG)
        print(f"Starting to estimate NS invars of {G}...")
        NSinvars_estimates = Compute_Rumin_Estimates(G, hodge_duality=hodge_duality)
        NS_str = []
        if Latex:   out_str += f"${G}$"
        for k in range(G.dim):
            if k not in NSinvars_estimates: 
                NS_str.append("-")
            elif NSinvars_estimates[k][0] == NSinvars_estimates[k][1]: 
                if Latex:   NS_str.append(f"{NSinvars_estimates[k][0]}")
                else: NS_str.append("{}".format(NSinvars_estimates[k][0]))
            else: 
                NS_str.append("[{}, {}]".format(*NSinvars_estimates[k]))
        if Latex: out_str += "".join([f" & ${a}$" for a in NS_str]) + " & $\infty^+$" + " & "*(6-G.dim) + "\\\\ \\hline \n"
        else: out_str += "| {} | {} | {}\n".format( " "*(COL_SEP-len(str(G))) + str(G), 
                        " | ".join( " "*((COL_SEP-len(s)+1)//2) +s + " "*((COL_SEP-len(s))//2) for s in NS_str),
                        " | ".join(" "*(COL_SEP) for i in range(MAX_DIM-len(NS_str)+1)))
        print(f"{G} completed.")
    if Latex: out_str += "\end{tabular}"
    else:     out_str += " " + "‾"*((MAX_DIM+1)*(COL_SEP+3)-1)
    return out_str
##### End Print_Methods.py #####