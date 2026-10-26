from FLUE.compiled import flue_c as fc

def writeNERSC(filename, U, NS, NT):
    fc.writegaugefield_nersc_c(filename, U, NS, NS, NS, NT)
    return

def readNERSC(filename, NS, NT):
    return fc.readgaugefield_nersc_c(filename, NS, NS, NS, NT)
