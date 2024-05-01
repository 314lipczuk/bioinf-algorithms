import sys
from matplotlib import pyplot
import os
import unittest
import sys

# exercise taken from MIT OCW 6.047 Computational Biology course, assignment 1, exercise 2
# along with makeDotplot fn

def find_matches(seq1, seq2, matching_Nmers=30):
  lookup = {}
  assert len(seq1) > matching_Nmers and len(seq2) > matching_Nmers
  for i in range(len(seq1) - matching_Nmers +1):
    lookup.setdefault(seq1[i:i+matching_Nmers], []).append(i)
  hits = []
  for i in range(len(seq2) - matching_Nmers + 1 ):
    hash = seq2[i:i+matching_Nmers]
    for hit in lookup.get(hash, []):
      hits.append((hit, i))
  return hits


# their plotting code, TODO: write my own
# def makeDotplot(filename, hits):
#     """generate a dotplot from a list of hits
#        filename may end in the following file extensions:
#          *.ps, *.png, *.jpg
#     """
#     x, y = zip(* hits)

#     slope1 = 1.0e6 / (825000 - 48000)
#     slope2 = 1.0e6 / (914000 - 141000)
#     offset1 = 0 - slope1*48000
#     offset2 = 0 - slope2*141000

#     hits2 = quality(hits)
#     print("%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits))))

#     # create plot
#     p = plotting.Gnuplot()
#     p.enableOutput(False)
#     p.plot(x, y, xlab="sequence 2", ylab="sequence 1")
#     p.plotfunc(lambda x: slope1 * x + offset1, 0, 1e6, 1e5)
#     p.plotfunc(lambda x: slope2 * x + offset2, 0, 1e6, 1e5)

#     # set plot labels
#     p.set(xmin=0, xmax=1e6, ymin=0, ymax=1e6)
#     p.set(main="dotplot (%d hits, %.5f%% hits on diagonal)" %
#           (len(hits), 100 * len(hits2) / float(len(hits))))
#     p.enableOutput(True)

#     # output plot
#     p.save(filename)

#     return p

def read_fasta(fn):
  s = ''
  with open(fn,'r') as f:
    s = ''.join(f.read().split('\n')[1:])
  return s

def quality(hits):
    """determines the quality of a list of hits"""

    slope1 = 1.0e6 / (825000 - 48000)
    slope2 = 1.0e6 / (914000 - 141000)
    offset1 = 0 - slope1*48000
    offset2 = 0 - slope2*141000

    goodhits = []

    for hit in hits:
        upper = slope1 * hit[0] + offset1
        lower = slope2 * hit[0] + offset2

        if lower < hit[1] < upper:
            goodhits.append(hit)

    return goodhits


def main():
  seq1 = read_fasta('../.static/ps1docs/human-hoxa-region.fa')
  seq2 = read_fasta('../.static/ps1docs/mouse-hoxa-region.fa')
  matches = find_matches(seq1, seq2)

def test():
  unittest.main()

class LocalAlignTester(unittest.TestCase):
  def test_findsbasic(self):
    s1 = ("a" * 30) + "dcb"
    s2 = ("a" * 30) + "bcd"
    matches = find_matches(s1, s2)
    eq = [(0,0)]
    self.assertEqual(len(eq), len(matches))
    for i, m in enumerate(matches):
      self.assertEqual(eq[i][0],m[0])
      self.assertEqual(eq[i][1],m[1])

  def test_findsame(self):
    s1 = ("a" * 30) + "cda"
    s2 = ("a" * 30) + "cda"
    matches = find_matches(s1, s2)
    eq = [(x,x) for x in range(4)]
    self.assertEqual(len(eq), len(matches))
    for i, m in enumerate(matches):
      self.assertEqual(eq[i][0],m[0])
      self.assertEqual(eq[i][1],m[1])

if __name__ == '__main__':
  test()