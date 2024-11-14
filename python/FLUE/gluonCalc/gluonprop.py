from FLUE.compiled import flue_c as fc

def calc_mom_space_scalarD(U, mu_start, xi):
    return fc.calc_mom_space_scalarD_c(U, mu_start, xi)
