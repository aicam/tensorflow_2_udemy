
"""
Python code for:
    * local alignment
    * global alignment with affine gaps
    * RNA folding

Example usage:

    From the command line run "python" (the code has been tested with python
    version 2.6 and 2.4). You will get the ">>>" prompt from python. Import
    this code by issuing the following command:

    >>> from align423 import *

    You can then run the local alignment, global alignment with affine gap
    penalities and rna folding algorithms using commands as shown below
    (commands are shown on lines starting with >>>):

    >>> local_align("acgt", "cg", ScoreParam(gap=-5, match=10, mismatch=-5))
    Scoring: match = 10; mismatch = -5; gap_start = 0; gap_extend = -5
    A matrix =
              *     a     c     g     t
        *     0     0     0     0     0
        c     0     0    10     5     0
        g     0     0     5    20    15
    Optimal Score = 20
    Max location in matrix = (3, 2)


    >>> affine_align("acgt", "cg", ScoreParam(gap_start=-15, gap=-5, match=10, mismatch=-5))
    x = acgt & y = cg
    Scoring: match = 10; mismatch = -5; gap_start = -15; gap_extend = -5
    M matrix =
              *     a     c     g     t
        *     0  -inf  -inf  -inf  -inf
        c  -inf    -5   -10   -30   -35
        g  -inf   -25   -10     0   -35
    X matrix =
              *     a     c     g     t
        *     0  -inf  -inf  -inf  -inf
        c   -20   -40   -45   -50   -55
        g   -25   -25   -30   -50   -55
    Y matrix =
              *     a     c     g     t
        *     0   -20   -25   -30   -35
        c  -inf   -40   -25   -30   -35
        g  -inf   -45   -45   -30   -20
    Optimal = -20


    >>> OPT, Arrow = rnafold("gaaccccttt")
    OPT Matrix =
              g     a     a     c     c     c     c     t     t     t
        g     0     0     0     0     0     0     0     0     0     0
        a     0     0     0     0     0     0     0     0     0     0
        a     0     0     0     0     0     0     0     0     0     0
        c     0     0     0     0     0     0     0     0     0     0
        c     0     0     0     0     0     0     0     0     0     0
        c     1     0     0     0     0     0     0     0     0     0
        c     1     0     0     0     0     0     0     0     0     0
        t     1     1     1     0     0     0     0     0     0     0
        t     2     2     1     0     0     0     0     0     0     0
        t     2     2     1     0     0     0     0     0     0     0
    Arrow Matrix =
              g     a     a     c     c     c     c     t     t     t
        g     0     0     0     0     0     0     0     0     0     0
        a     0     0     0     0     0     0     0     0     0     0
        a     0     0     0     0     0     0     0     0     0     0
        c     0     0     0     0     0     0     0     0     0     0
        c     0     0     0     0     0     0     0     0     0     0
        c     0     0     0     0     0     0     0     0     0     0
        c     0    -1     0     0     0     0     0     0     0     0
        t     2     2     2     0     0     0     0     0     0     0
        t     1     1     2    -1     0     0     0     0     0     0
        t     1     1     2    -1    -1     0     0     0     0     0
    >>> rna_backtrace(Arrow)
    Matched Pairs = (1,9), (2,8)
"""

#=============================================================
# Printing and utility functions
#=============================================================

Infinity = float('inf')


def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in xrange(sizex)]


def print_matrix(x, y, A):
    """Print the matrix with the (0,0) entry in the top left
    corner. Will label the rows by the sequence and add in the
    0-row if appropriate."""

    # decide whether there is a 0th row/column
    if len(x) == len(A):
        print "%5s" % (" "),
    else:
        print "%5s %5s" % (" ","*"),
        y = "*" + y

    # print the top row
    for c in x:
        print "%5s" % (c),
    print

    for j in xrange(len(A[0])):
        print "%5s" % (y[j]),
        for i in xrange(len(A)):
            print "%5.0f" % (A[i][j]),
        print


def is_complement(a, b):
    """Return True if character a is complmentary to character b"""
    assert len(a) == len(b) == 1
    return (a.upper(), b.upper()) in [
        ("A", "T"), ("T", "A"),
        ("C", "G"), ("G", "C"),
        ("A", "U"), ("U", "A")
    ]


#=============================================================
# Alignment Parameters
#=============================================================

class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    def __init__(self, match, mismatch, gap, gap_start=0):
        self.gap_start = gap_start
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

    def matchchar(self, a,b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a==b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (
                self.match, self.mismatch, self.gap_start, self.gap
        )

#=============================================================
# Sequence Alignment
#=============================================================

def local_align(x, y, score=ScoreParam(10, -5, -7)):
    """Do a local alignment between x and y with the given scoring parameters.
    We assume we are MAXIMIZING."""

    # create a zero-filled matrix
    A = make_matrix(len(x) + 1, len(y) + 1)

    best = 0
    optloc = (0,0)

    # fill in A in the right order
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

            # the local alignment recurrance rule:
            A[i][j] = max(
               A[i][j-1] + score.gap,
               A[i-1][j] + score.gap,
               A[i-1][j-1] + score.matchchar(x[i-1], y[j-1]),
               0
            )

            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)

    print "Scoring:", str(score)
    print "A matrix ="
    print_matrix(x, y, A)
    print "Optimal Score =", best
    print "Max location in matrix =", optloc
    # return the opt score and the best location
    return best, optloc


def affine_align(x, y, score=ScoreParam(10, -2, -7, -15)):
    """Global alignment with affine penalties. We assume we are maximizing."""
    M = make_matrix(len(x) + 1, len(y) + 1)
    X = make_matrix(len(x) + 1, len(y) + 1)
    Y = make_matrix(len(x) + 1, len(y) + 1)

    for i in xrange(1, len(x)+1):
        M[i][0] = -Infinity
        X[i][0] = -Infinity
        Y[i][0] = score.gap_start + i * score.gap

    for i in xrange(1, len(y)+1):
        M[0][i] = -Infinity
        X[0][i] = score.gap_start + i * score.gap
        Y[0][i] = -Infinity

    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

            M[i][j] = score.matchchar(x[i-1], y[j-1]) + max(
                    M[i-1][j-1],
                    X[i-1][j-1],
                    Y[i-1][j-1]
            )

            X[i][j] = max(
                    score.gap_start + score.gap + M[i][j-1],
                    score.gap + X[i][j-1],
                    score.gap_start + score.gap + Y[i][j-1]
            )

            Y[i][j] = max(
                    score.gap_start + score.gap + M[i-1][j],
                    score.gap_start + score.gap + X[i-1][j],
                    score.gap + Y[i-1][j]
            )

    opt = max(M[len(x)][len(y)], X[len(x)][len(y)], Y[len(x)][len(y)])

    print "x = %s & y = %s" % (x,y)
    print "Scoring:", str(score)
    print "M matrix ="
    print_matrix(x,y,M)
    print "X matrix ="
    print_matrix(x,y,X)
    print "Y matrix ="
    print_matrix(x,y,Y)
    print "Optimal =", opt

    return opt


#=============================================================
# RNA Folding
#=============================================================

def rnafold(rna):
    """Compute the dynamic programming matrix for the RNA folding
    algorithm for the given sequence rna."""
    n = len(rna)
    OPT = make_matrix(n, n)
    Arrows = make_matrix(n, n)

    for k in xrange(5, n):     # interval length
        for i in xrange(n-k):  # interval start
            j = i + k          # interval end

            # start with the values assuming j is not paired
            best_t = OPT[i][j-1]
            arrow = -1

            # search for the t that gives the best score
            for t in xrange(i, j-4):
                # only allow pairing between A-U and G-C
                if is_complement(rna[t], rna[j]):
                    if t > i:
                        val = 1 + OPT[i][t-1] + OPT[t+1][j-1]
                    # handle the case when we move past the main diagonal
                    else:
                        val = 1 + OPT[t+1][j-1]
                    if val >= best_t:
                       best_t = val
                       arrow = t

            OPT[i][j] = best_t
            Arrows[i][j] = arrow

    print "OPT Matrix ="
    print_matrix(rna, rna, OPT)
    print "Arrow Matrix ="
    print_matrix(rna, rna, Arrows)

    return OPT, Arrows


def rna_backtrace(Arrows):
    """Trace the RNA folding matrix (returned from rnafold) backward to find the
    actual pairs that are bonded."""
    Pairs = []  # holds the pairs in the optimal solution
    Stack = [(0, len(Arrows) - 1)]  # tracks where we have visited so far

    # while there are more items to visit
    while len(Stack) > 0:

        # take next cell off of list
        i, j = Stack.pop()

        # if cell is base case, skip it
        if j - i <= 4: continue

        # Arrow = -1 means we didn't match j
        if Arrows[i][j] == -1:
            Stack.append((i, j - 1))
        else:
            t = Arrows[i][j]
            assert j-4 > t >= i
            Pairs.append((t, j))  # save that j matched with t

            # add the two daughter problems
            if t > i: Stack.append((i, t - 1))
            Stack.append((t + 1, j - 1))
    print "Matched Pairs =",
    print ", ".join("(%d,%d)" % ij for ij in Pairs)
    return Pairs