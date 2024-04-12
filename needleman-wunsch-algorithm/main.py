from enum import Enum
import sys

# current version only works with simple gap costs f.ex G = -2 is the same option as querying BLAST with gap costs {existence:0, extension:2} - i might add support for more complex types in the future

S = [
	# A  G   C   T
	[1, -2, -2, -2], # A
	[-2, 1, -2, -2], # G
	[-2, -2, 1, -2], # C
	[-2, -2, -2, 1]  # T
	]
G = -2 # gap pen
base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }

class NWConfig:
  def from_matrix(self, m):
    self.matrix = m
    return self

  def from_rules(self, r, size=4):
    keys = r.keys()
    assert 'match' in keys and 'mismatch' in keys
    assert size > 0
    self.matrix = [[r['match'] if i==j else r['mismatch'] for i in range(size)] for j in range(size)]
    return self

  def with_penalty(self, p):
    self.penalty = p
    return self

  def lookup(self, a, b):
    return self.matrix[base_idx[a]][base_idx[b]]

  def fill(self, m, seq1, seq2, aggr_fun=max):
    for j in range(len(seq2)):
      for i in range(len(seq1)):
        self.fill_single(m, [j+1, i+1], seq1, seq2, aggr_fun=aggr_fun)

  def fill_single(self, m, pos, seq1, seq2, aggr_fun):
    j,i = pos
    assert i > 0 and j > 0
    l_j, l_i = [seq2[j-1], seq1[i-1]]
    score = aggr_fun([
      [(m[j-1][i-1][0] + self.lookup(l_j, l_i)), Direction.DIAG],
      [(m[j-1][i][0] + self.penalty), Direction.UP],
      [(m[j][i-1][0] + self.penalty), Direction.LEFT],
      ], key=lambda x: x[0])
    m[j][i] = score


class Direction(Enum):
  DONE = 0
  DIAG = 1
  UP = 2
  LEFT = 3

def make_initial_matrix(seq1, seq2, config):
  m = [[[0,None] for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
  for j in range(len(seq1)+1): m[0][j] = [config.penalty * j, Direction.LEFT]
  for i in range(len(seq2)+1): m[i][0] = [config.penalty * i, Direction.UP]
  m[0][0][1] = Direction.DONE
  return m



def solve(m, seq1, seq2):
  solution = []
  j, i = len(m)-1, len(m[0])-1
  score = m[j][i][0]
  while True:
    current_node = m[j][i]
    match current_node[1]:
      case Direction.DONE:
        break
      case Direction.DIAG:
        solution.append((seq1[i-1], seq2[j-1]))
        j -= 1
        i -= 1
      case Direction.LEFT:
        solution.append((seq1[i-1], ' '))
        i -= 1
      case Direction.UP:
        solution.append((' ', seq2[j-1]))
        j -= 1
  solution.reverse()
  l1,l2 = "",""
  for s in solution:
    l1 += s[0]
    l2 += s[1]
  return l1,l2, score

def print_matrix(m):
  for r in m:
   print(r,'\n') 

def read_fasta(f):
  return ''.join(f.split('\n')[1:])

def test():
  """testing procedure, utilizing human_HoxA13.fa and mouse_HoxA13.fa - compared with score from blast.ncbi.nlm.nih.gov"""
  import os
  if os.getcwd().split('/')[-1] != 'needleman-wunsch-algorithm':
    print("run test from within directory needleman-wunsch-algorithm, please")
    sys.exit(1)
  seq1 ,seq2 = "",""
  with open("./../static/test/needleman_wunsch/human_HoxA13.fa",'r') as f:
    seq2 = read_fasta(f.read())
  with open('./../static/test/needleman_wunsch/mouse_HoxA13.fa','r') as f:
    seq1 = read_fasta(f.read())
  config = NWConfig().from_matrix(S).with_penalty(G)
  R = make_initial_matrix(seq1, seq2, config)
  config.fill(R, seq1, seq2)
  _,_,score = solve(R, seq1, seq2)
  assert score == 840
  print("test 1 - alignment score for HoxA13 in human and mouse the same as from blast\npassed")

  R = make_initial_matrix(seq1, seq2, config)
  config = NWConfig().from_rules({'match':1, 'mismatch':-2}).with_penalty(G)
  config.fill(R, seq1, seq2)
  _,_,score2 = solve(R, seq1, seq2)
  assert score == score2
  print("test 2 - score computed using rule-based type of config matches the original\npassed")

  R = make_initial_matrix(seq1, seq1, config)
  config = NWConfig().from_rules({'match':0, 'mismatch':1}).with_penalty(-G)
  config.fill(R, seq1, seq1, aggr_fun=min)
  _,_,score3 = solve(R, seq1, seq1)
  print(score3)
  assert score3 == 0
  print("test 3 - score distance metric achieved by changing matrix and aggregation function\npassed")

def main():
  if len(sys.argv) == 2 and sys.argv[-1] == 'test':
    test()
    sys.exit(0)
  if len(sys.argv) -1 != 2:
    print("wrong arguments. Call this like this:\n./program <fasta_file1> <fasta_file2>")
    sys.exit(1)
  seq1,seq2 = "",""
  with open(sys.argv[1],'r') as f:
    seq2 = read_fasta(f.read())
  with open(sys.argv[2],'r') as f:
    seq1 = read_fasta(f.read())
  config = NWConfig().from_matrix(S).with_penalty(G)
  R = make_initial_matrix(seq1, seq2, config)
  config.fill(R, seq1, seq2)
  print(solve(R, seq1, seq2))
  
if __name__ == "__main__":
  main()

# how do you make a distance metric?
# well, in theory i could do it by substituting the aggregating function to min and the S matrix into diagonal zeros and all the rest positive ints
# what about gaps? change into positive numbers i guess
# aight, let's try that
