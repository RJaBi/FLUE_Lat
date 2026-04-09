from FLUE.compiled import flue_c as fc


def writeILDG(filename, U, NS, NT):
    fc.writegaugefield_ildg_c(filename, U, NS, NS, NS, NT)
    return

def readILDG(filename, NS, NT):
    return fc.readgaugefield_ildg_c(filename, NS, NS, NS, NT)

