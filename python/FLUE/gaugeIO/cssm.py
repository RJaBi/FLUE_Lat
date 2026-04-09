from FLUE.compiled import flue_c as fc

def readCSSM(filename, NS, NT):
    return fc.readgaugefield_cssm_c(filename, NS, NS, NS, NT)

