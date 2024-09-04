from FLUE.compiled import flue_c as fc

def q_average(NQ, q, D, q_slices, DOUT, QOUT, Qcount):
    return fc.q_average(NQ, q, D, q_slices, DOUT, QOUT, Qcount)

def cone_cut(NQ, radius, Q, D, D4, NT, NS, angleIN, xiIN, IRCutIN, IRRadiusIN):
    # Should be QOUT, Dout, D4out, qcount
    return fc.cone_cut(NQ, radius, Q, D, D4, NT, NS, angleIN, xiIN, IRCutIN, IRRadiusIN)


#def calculate_area(radius):
#    return fc.calculate_area_c(radius)


#def calculate_perimeter(radius):
#    return fc.calculate_perimeter_c(radius)
