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

                #s = {'align_cost': f(sequences[i], f(sequences[j])),
                #     'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                #         len(sequences[i]), align_length, ',BANDED' if banded else ''),
                #     'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                #         len(sequences[j]), align_length, ',BANDED' if banded else '')}
                s = {'align_cost': f(sequences[i], sequences[j]),
                    'seqi_first100': 'stuff',
                    'seqj_first100': 'more_stuff'}
                jresults.append(s)
            results.append(jresults)
        return results
