from FLUE.compiled import flue_c as fc


def get_qhat(coord, shape, a=1):
    return fc.get_qhat_c(coord, shape, a)
