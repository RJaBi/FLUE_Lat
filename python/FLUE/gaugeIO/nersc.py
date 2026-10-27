from FLUE.compiled import flue_c as fc
from datetime import datetime

def writeNERSC(filename, U, NS, NT):
    fc.writegaugefield_nersc_c(filename, U, NS, NS, NS, NT)
    return

def readNERSC(filename, NS, NT):
    return fc.readgaugefield_nersc_c(filename, NS, NS, NS, NT)

def writeNERSC_fullHeader(filename, U, NS, NT, 
                        sequence_number,
                        creator,
                        creator_hardware,
                        creation_date,
                        ):
    now = datetime.now()
    archive_date = now.strftime("%A %B %d %H:%M:%S %Y")
    # It gets NS, NT from U automatically
    # IT worked above because they are default arguments at the end
    fc.writegaugefield_nersc_fullheader_c(
        filename, U,
        '1.0',
        'IEEE64LITTLE',
        '4D_SU3_3X3',
        '1234567',
        sequence_number,
        creator,
        creator_hardware,
        creation_date,
        archive_date,
        '',
        '')
