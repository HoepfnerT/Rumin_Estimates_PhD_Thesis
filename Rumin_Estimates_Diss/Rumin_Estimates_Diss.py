##### Start Rumin_Estimates_Diss.py #####
import Print_Methods
from NLG_List import * # Here the Lie algebras are defined

if __name__=="__main__":
    RE = Print_Methods.make_NS_invars_table(ALL_GNLG_LEQ6, MAX_DIM=6, Latex = False)
    print(RE)
    #
    # For some arbitrarely chosen nilpotent Lie algebras of dimension 7:
    #print(Print_Methods.make_NS_invars_table([H7, L7_37A, L7_257C], MAX_DIM=7))
##### End Rumin_Estimates_Diss.py #####


