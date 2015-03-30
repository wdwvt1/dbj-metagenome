
from numpy import array
from dbjm.dbj import Nodes2
from dbjm.analyze import traverse_within_radius


pairs = [('ACTGA', array([5,0,0,4])),
         ('CTGAA', array([0,5,3,0])),
         ('TGAAC', array([0,5,4,3])),
         ('GAACC', array([0,0,0,0])),
         ('GAACT', array([0,0,0,0])),
         ('GAACG', array([0,8,1,0])),
         ('TGAAT', array([0,0,0,0])),
         ('CTGAG', array([3,0,0,0])),
         ('TGAGA', array([0,0,0,0]))]


n = Nodes2()

for kmer, arr in pairs:
    n.add_from_file(kmer, arr)

traverse_within_radius('ACTGA', 3, n, 'AAAA', 4, 1, alphabet='ACTG')