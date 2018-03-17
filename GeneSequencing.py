#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time

def SUB = 1
def INDEL = 5
def MATCH = -3

def check_banded(row, col):
    bandwidth = 7

    if row == col:
        return True

    top = col
    btm = row

    if row > col:
        top = row
        btm = col

    while bandwidth > 0:
        top = top - 1
        btm = btm + 1
        if top <= btm:
            return True
        bandwidth = bandwidth - 1

    return False


def unrestriced(gene1_str, gene2_str):
    zero_list = [0]
    lst_one = zero_list + list(gene1_str)
    lst_two = zero_list + list(gene2_str)
    matrix = [[None for x in range(len(lst_one))] for y in range(len(lst_two))

    for i in range(len(lst_one)):
        for j in rand(len(lst_two)):






class GeneSequencing:
    def __init__( self ):
        pass

    def align_all( self, sequences, banded, align_length ):
        #print(banded)
        #print(align_length)
        results = []
        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                s = {'align_cost':i+j,
                     'seqi_first100':'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                         len(sequences[i]), align_length, ',BANDED' if banded else ''),
                     'seqj_first100':'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                         len(sequences[j]), align_length, ',BANDED' if banded else '')}
                jresults.append(s)
            results.append(jresults)
        return results
