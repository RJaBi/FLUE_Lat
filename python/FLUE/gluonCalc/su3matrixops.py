from FLUE.compiled import flue_c as fc

def multiplymatmat(left, right):
    return fc.multiplymatmat_c(left, right)
