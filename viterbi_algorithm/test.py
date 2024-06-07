import main

with open('toy_example') as f:
  global X
  X = [main.base_idx[i] for i in list(f.read())[:-1]]

print(main.viterbi(X))

