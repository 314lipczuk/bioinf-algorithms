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

def make_initial_matrix(seq1, seq2):
  m = [[0 for _ in range(len(seq1)+1)] for _ in range(len(seq2)+1)]
  for j in range(len(seq1)+1): m[0][j] = -10 * j
  for i in range(len(seq2)+1): m[i][0] = -10 * i
  return m

R = make_initial_matrix(seq1, seq2)

print(R)

def fill_single(m, pos):
  j,i = pos
  assert i > 0 and j > 0
  l_j, l_i = [seq2[j-1], seq1[i-1]]
  score = max([
    m[j-1][i-1] + S[base_idx[l_j]][base_idx[l_i]],
    m[j-1][i] + G,
    m[j][i-1] + G,
    ])
  m[j][i] = score

def fill(m):
  for j in range(len(seq2)):
    for i in range(len(seq1)):
      fill_single(m, [j+1, i+1])

fill(R)
def print_matrix(m):
  for r in m:
   print(r) 

print_matrix(R)
