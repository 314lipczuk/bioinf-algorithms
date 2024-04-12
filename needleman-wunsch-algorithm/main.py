from enum import Enum
import sys

# current version only works with simple gap costs - G = -2 is the same option as querying BLAST with gap costs {existence:0, extension:2} - i might add support for more complex types in the future

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

class Direction(Enum):
  DONE = 0
  DIAG = 1
  UP = 2
  LEFT = 3

def make_initial_matrix(seq1, seq2):
  m = [[[0,None] for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
  for j in range(len(seq1)+1): m[0][j] = [G * j, Direction.LEFT]
  for i in range(len(seq2)+1): m[i][0] = [G * i, Direction.UP]
  m[0][0][1] = Direction.DONE
  return m


def fill_single(m, pos, seq1, seq2):
  j,i = pos
  assert i > 0 and j > 0
  l_j, l_i = [seq2[j-1], seq1[i-1]]
  score = max([
    [(m[j-1][i-1][0] + config.lookup(l_j, l_i)), Direction.DIAG],
    [(m[j-1][i][0] + G), Direction.UP],
    [(m[j][i-1][0] + G), Direction.LEFT],
    ], key=lambda x: x[0])
  m[j][i] = score

def fill(m, seq1, seq2):
  for j in range(len(seq2)):
    for i in range(len(seq1)):
      fill_single(m, [j+1, i+1], seq1, seq2)


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
  R = make_initial_matrix(seq1, seq2)
  global config
  config = NWConfig().from_matrix(S).with_penalty(R)
  fill(R, seq1, seq2)
  _,_,score = solve(R, seq1, seq2)
  assert score == 840
  print("test passed")

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
  R = make_initial_matrix(seq1, seq2)
  global config
  config = NWConfig().from_matrix(S).with_penalty(R)
  fill(R, seq1, seq2)
  print(solve(R, seq1, seq2))
  
if __name__ == "__main__":
  main()