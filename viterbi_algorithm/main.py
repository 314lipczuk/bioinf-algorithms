import sys
from math import log

# DISCLAIMER: the setup code was taken from OCW 6.047/6.878/HST.507 Fall 2015, assignment 1, task 3
# implementation is my own

###############################################################################
# HMM PARAMETERS
# Conventions: + refers to High-GC, and - refers to Low-GC. When indexing
#  states, 0 is + and 1 is -.
###############################################################################

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
state_idx = { '+' : 0, '-' : 1 }

# initial distribution over states, i.e. probability of starting in state k
init_dist = [0.5,0.5]

# transition probabilities
tr = [
    #  to+   to-
    [ 0.99, 0.01 ], # from+
    [ 0.01, 0.99 ]  # from-
]


# emission probabilities
em = [
    #    A     G     C     T
    [ 0.20, 0.30, 0.30, 0.20], # +
    [ 0.30, 0.20, 0.20, 0.30]  # -
]

###############################################################################
# VITERBI ALGORITHM (you must complete)
# Note: The length of the sequences we are dealing with is large enough that it
#       is necessary to use log-probabilities for numerical stability. You will
#       need to adapt the formulae accordingly.
###############################################################################

def viterbi(X):
    """Returns the Viterbi path for the emission sequence X.
    X should be a list of integers, 0=A, 1=G, 2=C, 3=T.
    The returned Y is a list of integers, 0=High-GC, 1=Low-GC.
    """

    N = len(tr)
    L = len(X)
    assert len(em) == N

    V = [[0]*N for _ in range(L)]
    TB = [[0]*N for _ in range(L)]
    
    for i in range(0,L):
        Vprev = []
        if i == 0:
            Vprev = [log(pk0) for pk0 in init_dist]
        else:
            Vprev = V[i-1]

        for k in range(N):
            # Set V[i][k] to the appropriate value for the Viterbi matrix, based
            #  on Vprev (V[i-1]) and the model parameters.
            # based on transition probabilities for current round, and max of the results of the previous round
            # Vprev, params
            """
            e_l(i) * max_k(p_k(j,x-1) * p_kl)           

            e_l - prob to observe element i in state l
            p_k - probability of the most probable path ending at position x-1 in state k with element j
            p_pl - probability of the transition from state l to state k
            """
            

            c2 = max( [] )
            V[i][k] = 

            # Set TB[i][k] to the selected previous state (0 or 1 corresponding
            #  to + or -)
            # k - transition
            # i - state
            # get max of possible 

    # perform traceback and return the predicted hidden state sequence
    Y = [-1 for i in range(L)]
    _, yL = max([ (V[L-1][k], k) for k in range(N)])
    Y[L-1] = yL
    for i in range(L-2,-1,-1):
        Y[i] = TB[i+1][Y[i+1]]
    return Y

###############################################################################
# ANNOTATION BENCHMARKING
###############################################################################

def basecomp(X,anno):
    counts = [[0]*4,[0]*4]
    for i in range(len(X)):
        counts[anno[i]][X[i]] += 1
    sum0 = sum(counts[0])
    sum1 = sum(counts[1])
    return [[float(x)/sum0 for x in counts[0]],[float(x)/sum1 for x in counts[1]]]

def region_lengths(anno):
    lengths = [[],[]]
    curlen=1
    for i in range(1,len(anno)):
        if anno[i] == anno[i-1]:
            curlen += 1
        else:
            lengths[anno[i-1]].append(curlen)
            curlen=1
    lengths[anno[len(anno)-1]].append(curlen)
    return lengths

def anno_accuracy(refanno,testanno):
    correct = 0
    assert len(refanno) == len(testanno)
    for i in range(1,len(refanno)):
        if refanno[i] == testanno[i]:
            correct += 1
    return float(correct)/len(refanno)

def print_basecomp(b):
    print( "A=%.2f%% G=%.2f%% C=%.2f%% T=%.2f%%" % (100*b[0],100*b[1],100*b[2],100*b[3]))

def print_annostats(X,anno,filename):
    lengths = region_lengths(anno)
    basecomps = basecomp(X,anno)

    print( "High-GC mean region length: ", sum(lengths[0])/len(lengths[0]))
    print( "High-GC base composition:",)
    print_basecomp(basecomps[0])
    print( "Low-GC mean region length: ", sum(lengths[1])/len(lengths[1]))
    print( "Low-GC base composition:",)
    print_basecomp(basecomps[1])

    print( "Saving High-GC length histogram to %s_highgc.png" % filename)
    p = plothist(lengths[0],low=0)
    p.save(filename+"_highgc.png")
    print( "Saving Low-GC length histogram to %s_lowgc.png" % filename)
    p = plothist(lengths[1],low=0)
    p.save(filename+"_lowgc.png")

###############################################################################
# MAIN
###############################################################################

def main():
    if len(sys.argv) < 2:
        print( "you must call program as: ./viterbi.py <datafile>")
        sys.exit(1)
    
    datafile = sys.argv[1]

    f = open(datafile)
    X = f.readline()
    refanno = f.readline()
    f.close()

    if X[len(X)-1] == '\n': X=X[0:len(X)-1]
    if refanno[len(refanno)-1] == '\n': refanno=refanno[0:len(refanno)-1]

    X=list(X)
    refanno=list(refanno)
    for i in range(len(X)):
        X[i] = base_idx[X[i]]
    for i in range(len(refanno)):
        refanno[i] = state_idx[refanno[i]]

    print( "Authoritative annotation statistics")
    print( "-----------------------------------")
    print_annostats(X,refanno,datafile+"_authoritative")
    print( "")
    vanno = viterbi(X)
    print( "Viterbi annotation statistics")
    print( "-----------------------------")
    print_annostats(X,vanno,datafile+"_viterbi")
    print( "")
    print( "Accuracy: %.2f%%" % (100*anno_accuracy(refanno,vanno)))

if __name__ == "__main__":
    main()