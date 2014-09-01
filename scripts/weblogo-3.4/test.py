
from corebio.matrix import Motif


fin = open('WM')
Motif.read_swissRegulon( fin, alphabet='ACGT' )
