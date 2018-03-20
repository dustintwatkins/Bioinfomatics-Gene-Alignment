#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


import time
import sys

def memoize(f):
    memo = {}
    def helper(*args):
        key = tuple(args)
        try:
            return memo[key]
        except:
            memo[key] = f(*args)
            return memo[key]
    return helper

def align_genes_helper(string1, string2, dist):
    str1 = ['-'] + list(string1)
    str2 = ['-'] + list(string2)
    path = {}
    def recurse(function, i, j):
        if abs(i-j) > dist:
            return (float("inf"), path)
        if i == 0:
            return (j*5, path)
        if j == 0:
            return (i*5, path)
        diagonal = function(function, i - 1, j - 1)[0]
        if str1[i] != str2[j]:
            diagonal += 1
        if str1[i] == str2[j]:
            diagonal += -3
        indel_left = function(function, i - 1, j)[0] + 5
        indel_up = function(function, i, j - 1)[0] + 5
        min_score = min(diagonal, indel_left, indel_up)
        if min_score == diagonal:
            path[(i,j)] = (i-1, j-1)
        if min_score == indel_left:
            path[(i,j)] = (i-1, j)
        if min_score == indel_up:
            path[(i,j)] = (i, j-1)
        return (min_score, path)
    memo_recurse = memoize(recurse)
    return memo_recurse(memo_recurse, len(string1), len(string2))

def align_genes(string1, string2):
    return align_genes_helper(string1, string2, float("inf"))

def banded_edit(string1, string2):
    return align_genes_helper(string1, string2, 3)

class GeneSequencing:
    def __init__( self ):
        pass

    def align_all( self, sequences, banded, align_length ):
        #print(banded)
        #print(align_length)
        results = []

        sys.setrecursionlimit(10**6)
        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                f = align_genes
                if banded:
                    f = banded_edit


                a = sequences[i][:align_length]
                b = sequences[j][:align_length]
                res = f(a, b)

                if res[0] == float("inf"):
                    jresults.append({'align_cost': 'No alignment possible',
                        'seqi_first100': 'None',
                        'seqj_first100': 'None'})
                    continue

                path = res[1]

                a = ['-'] + list(a)
                b = ['-'] + list(b)

                x = len(a) - 1
                y = len(b) - 1
                seqi_100 = []
                seqj_100 = []

                while x != 0 and y != 0:
                    next = path[(x,y)]
                    if x - 1 == next[0] and y == next[1]:
                        seqi_100.append("-")
                        seqj_100.append(b[y])
                        x = x - 1
                        continue
                    elif x == next[0] and y -1 == next[1]:
                        y = y - 1
                    else:
                        seqi_100.append(a[x])
                        seqj_100.append(b[y])
                        x = x - 1
                        y = y - 1

                s = {'align_cost': res[0],
                    'seqi_first100': ''.join(seqi_100[len(seqi_100) - 100:][::-1]),
                    'seqj_first100': ''.join(seqj_100[len(seqj_100) - 100:][::-1])}
                jresults.append(s)
            results.append(jresults)
        return results
