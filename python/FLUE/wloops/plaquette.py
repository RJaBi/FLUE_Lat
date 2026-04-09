from FLUE.compiled import flue_c as fc
import numpy as np

def plaq(U):
    """
    Average of all orientation of plaquettes
    """
    sumtrp, nplaq, time = fc.genplaquette_c(U, 1, 4, 4)
    return sumtrp / float(nplaq)

def splaq(U):
    """
    Average of space-space plaquettes
    """
    sumtrp, nplaq, time = fc.genplaquette_c(U, 2, 4, 4)
    return sumtrp / float(nplaq)

def tplaq(U):
    """
    Average of space-time plaquettes
    """
    sumtrp, nplaq, time = fc.genplaquette_c(U, 1, 1, 4)
    return sumtrp / float(nplaq)
