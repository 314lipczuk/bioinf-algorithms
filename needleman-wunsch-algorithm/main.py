from enum import Enum


S = [
	# A  G   C   T
	[3, -1, -2, -2], # A
	[-1, 3, -2, -2], # G
	[-2, -2, 3, -1], # C
	[-2, -2, -1, 3]  # T
	]

S = [
	# A  G   C   T
	[1, -1, -1, -1], # A
	[-1, 1, -1, -1], # G
	[-1, -1, 1, -1], # C
	[-1, -1, -1, 1]  # T
	]

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
G = -4 # gap pen
G = -2 # gap pen

seq1, seq2 = "GCATGCG", "GATTACA"
seq1, seq2 = "AATCG", "AACG"

class Direction(Enum):
  DONE = 0
  DIAG = 1
  UP = 2
  LEFT = 3

def make_initial_matrix(seq1, seq2):
  m = [[[0,None] for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
  for j in range(len(seq1)+1): m[0][j] = [-10 * j, Direction.LEFT]
  for i in range(len(seq2)+1): m[i][0] = [-10 * i, Direction.UP]
  m[0][0][1] = Direction.DONE
  return m

R = make_initial_matrix(seq1, seq2)

#print(R)

def fill_single(m, pos):
  j,i = pos
  assert i > 0 and j > 0
  l_j, l_i = [seq2[j-1], seq1[i-1]]
  score = max([
    [(m[j-1][i-1][0] + S[base_idx[l_j]][base_idx[l_i]]), Direction.DIAG],
    [(m[j-1][i][0] + G), Direction.UP],
    [(m[j][i-1][0] + G), Direction.LEFT],
    ], key=lambda x: x[0])
  m[j][i] = score

def fill(m):
  for j in range(len(seq2)):
    for i in range(len(seq1)):
      fill_single(m, [j+1, i+1])


def solve(m):
  solution = []
  j, i = len(m)-1, len(m[0])-1
  #current_node = m[j][i]
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
  return l1,l2

def print_matrix(m):
  for r in m:
   print(r) 

fill(R)
print_matrix(R)
print(solve(R))