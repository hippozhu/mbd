import math
import swalign
from collections import namedtuple
import numpy as np
from sklearn.neighbors import NearestNeighbors

def alnDistFunc(a1, a2):
  #return math.fabs(a1.ref_start-a2.ref_start)
  return math.fabs(a1[0]-a2[0])

ref = open('shfv.lower').read()
qry = 'uccuuaacc'
#qry_short = 'uccuuaacc'
qry_long = 'gcagacccuccuuaaccauguucug'
starting_pos = 10000
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring) 
Alignment = namedtuple('Alignment', ['score', 'ref_start', 'ref_end'])

def align(ref, qry, sw):
  step_size = 10
  seg_len = 20
  segs = [(x, min(x+seg_len, len(ref))) for x in xrange(starting_pos, len(ref), step_size)]
  aligns = ((start, sw.align(ref[start: end], qry)) for start, end in segs)
  return sorted(set([Alignment(align.score, start + align.r_pos, start + align.r_end) for start, align in aligns]), key = lambda x:x.score, reverse=True)

algns = align(ref, qry, sw)
nbrs = NearestNeighbors(metric = alnDistFunc, algorithm = 'brute').fit([[aln.ref_start] for aln in algns])
know_trs = np.loadtxt('known_trs', delimiter = ',', usecols=[2], dtype = np.int64)
distances, indices = nbrs.kneighbors(know_trs.reshape(32,1), 1)
allscores =  np.array([algn.score for algn in algns])
scores = np.array([algns[i].score for i in indices[:,0]])

'''
top_results = sorted(aln_results, key = lambda result: result[0], reverse=True)[:20]
aln_results_filtered = filter(lambda ar: ar[1]>10000, aln_results)
top_alns = sorted(aln_results_filtered, key  = lambda alf:alf[0], reverse=True)[:20]
filter(lambda x:x[1]>13000 and x[1]<14000, sorted(aln_results_filtered, key  = lambda alf:alf[0], reverse=True)[:100])
'''
