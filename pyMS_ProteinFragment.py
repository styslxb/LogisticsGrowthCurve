import re
import random
import datetime
print """				pyMS v0.5
-----------------------------------------------------------------------------
This script was written to cut the protein sequence into small peptide
fragments and to calculate their molecular weight. For now it can calculate all
the fragment's Mw. derived form the source protein sequence. It needs to be
improved. It is free to use.Welcome! Any questions please email:
 styslxb@gmail.com.
-------------------------------------------------------------------------------"""
myProteinName = raw_input( 'Please input the protein name,max 8 chars,default is random chars:' )
myProteinSeq = raw_input( "Please input the protein sequence:" )
myCysteineForm = raw_input( 'Cysteine form: R for-SH, O for-S-S-, default is R[R/O?]:' )
myAtomAvgMw = {"H+":1.00739, "H":1.00794, "H2O":18.01524, }


def SeqCleanIllegalChar( PrSeq ):
    myRecogizePattern = re.compile( '[^ARNDCEQGHILKMFPSTWYV]', re.IGNORECASE )
    myOutput = myRecogizePattern.sub( '', PrSeq )
    return myOutput.upper()

def GenerateAminoAcidMwDict( CysteineForm , AtomMw ):
    CysteineForm = re.sub( '[^o|r]', '', CysteineForm, re.IGNORECASE )
    CysteineForm = CysteineForm.upper()
    CysteineForm = CysteineForm[0:1]
    myAminoAcidAvgMw = {'A':71.0788,
                        'R':156.1875,
                        'N':114.1038,
                        'D':115.0886,
                        'C':103.1388,
                        'E':129.1155,
                        'Q':128.1307,
                        'G':57.0519,
                        'H':137.1411,
                        'I':113.1594,
                        'L':113.1594,
                        'K':128.1741,
                        'M':131.1926,
                        'F':147.1766,
                        'P':97.1167,
                        'S':87.0782,
                        'T':101.1051,
                        'W':186.2132,
                        'Y':163.1760,
                        'V':99.1326,
                        'U':150.0388,
                        'O':237.3018}
    if CysteineForm == 'O':
        deltac = myAminoAcidAvgMw.get( 'C' ) - AtomMw.get( 'H' )
        myAminoAcidAvgMw.update( {'C':deltac} )
    return myAminoAcidAvgMw

def randomfilename():
    dictionary = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    pickedchars = random.sample( dictionary * 5, 8 )
    w = str( "" )
    for t in pickedchars:
        w = w + t
    return w

def GenerateFileName( ProteinName ):
    ProteinName = re.sub( '[\W]', '', ProteinName )
    d = datetime.datetime.now()
    part = d.strftime( '%Y%m%d%H%M%S' )
    ActualPrName = randomfilename()
    if ProteinName:
        ActualPrName = ProteinName[0:8]
    return 'SeqFrag_' + part + '_' + ActualPrName + '.csv'

def SeqFragmentsEm( PrSeq, AminoAcidMw, AtomMw ):
    f = open( GenerateFileName( myProteinName ), 'w' )
    FragmentCount = 0
    FragmentSize = 1
    SeqLength = len( PrSeq )
    while SeqLength > 0:
        if FragmentSize <= SeqLength:
            for TraversalTimes in range( 0, SeqLength ):
                print 'Processing the length of ' + str( FragmentSize ) + '-aa peptides.'
                for PossitionPointer in range( 0, SeqLength - TraversalTimes ):
                    FragmentCount += 1
                    FragmentStartPossition = PossitionPointer + 1
                    FragmentStopPossition = PossitionPointer + FragmentSize
                    FragmentSeq = PrSeq[PossitionPointer:FragmentStopPossition]
                    FragmentSeqMw = 0
                    for char in FragmentSeq:
                        FragmentSeqMw += AminoAcidMw.get( char )
                    FragmentSeqMw += AtomMw.get( 'H2O' )
#                     print FragmentSeqMw
#                     print FragmentCount, ",", PossitionPointer + 1, ',', FragmentSeq, ',', FragmentStopPossition, ',', FragmentSeqMw
                    f.write( str( FragmentCount ) + "," + str( FragmentStartPossition ) + ',' + FragmentSeq + ',' + str( FragmentStopPossition ) + ',' + str( FragmentSeqMw ) + "\n" )
                FragmentSize += 1
            SeqLength -= 1
        print "OK."
        break
    f.close()
SeqFragmentsEm( SeqCleanIllegalChar( myProteinSeq ), GenerateAminoAcidMwDict( myCysteineForm , myAtomAvgMw ), myAtomAvgMw )
