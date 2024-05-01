import sys
import os
from matplotlib import pyplot as plt
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

def plot_and_save(data, filename):
  x = list(map(lambda x:x[0], data))
  y = list(map(lambda x:x[1], data))
  fig = plt.figure()
  plot = fig.add_subplot()
  plot.scatter(x,y)
  fig.savefig(filename)


def main():
  seq1 = read_fasta('../.static/ps1docs/human-hoxa-region.fa')
  seq2 = read_fasta('../.static/ps1docs/mouse-hoxa-region.fa')
  matches = find_matches(seq1, seq2)
  plot_and_save(matches, 'figure.png')

if __name__ == '__main__': main()