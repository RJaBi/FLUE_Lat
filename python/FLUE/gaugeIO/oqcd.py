from FLUE.compiled import flue_c as fc

def writeOQCD(filename, U, NS, NT):
    fc.writegaugefield_openqcd_c(filename, U, NS, NS, NS, NT)
    return

def readOQCD(filename, NS, NT):
    return fc.readgaugefield_openqcd_c(filename, NS, NS, NS, NT)
