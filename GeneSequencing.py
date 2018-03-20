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

#Time Complexity:
#Space Complexity: O(n)^2
def memoize(f):
	memo = {}																	#Dictionary memo
	def helper(*args):															#Takes in a list of arguments and applys them to function
    	key = tuple(args)														#Set the dictionary key
    	try:																	#If a key exists, return the value associated w/ it
        	return memo[key]
    	except:																	#If dne, the value returned from the function parameter is returned
        	memo[key] = f(*args)
        	return memo[key]
	return helper

#Time Complexity:
#Space Complexity: O(n)^2
def align_genes_helper(string1, string2, dist):
	str1 = ['-'] + list(string1)												#Create list from string one and make index 0 '-'
	str2 = ['-'] + list(string2)												#Create list from string two and make index 0 '-'
	path = {}																	#Create a dictionary for the back trace, so we can get the path

	#Space Complexity:
	#Time Complexity:
	def recurse(function, i, j):

		if abs(i-j) > dist:														#If banded, we do not calculate when distance > 7 (3 from diagonal)
        	return (float("inf"), path)

    	if i == 0:																#Set row/columns to multiples of 5
        	return (j*5, path)

    	if j == 0:																#Set row/columns to multiples of 5
        	return (i*5, path)

    	diagonal = function(function, i - 1, j - 1)[0]							#Get diagonal value

    	if str1[i] != str2[j]:													#If string chars do NOT match
        	diagonal += 1														# +1 score

    	if str1[i] == str2[j]:													#If string chars DO match
        	diagonal += -3														# -3 score

    	indel_left = function(function, i - 1, j)[0] + 5						#Get value of left and +5 score
    	indel_up = function(function, i, j - 1)[0] + 5							#Get value of up and +5 score

		min_score = min(diagonal, indel_left, indel_up)							#Get the min of diagonal, up and left

    	if min_score == diagonal:												#If diagonal is min_score
        	path[(i,j)] = (i-1, j-1)											#Set path indexes mapped to diagonal indexes
    	if min_score == indel_left:												#If left is min_score
        	path[(i,j)] = (i-1, j)												#Set path indexex mapped to left indexes
    	if min_score == indel_up:												#If up is min_score
        	path[(i,j)] = (i, j-1)												#Set path indexes mapped to up indexes

		return (min_score, path)												#Return the min score and dictionary
	memo_recurse = memoize(recurse)												#Store the values for quick look up, instead of recalculating
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

		#Space: O(n)
    	results = []

    	sys.setrecursionlimit(10**6)											#This was because my VM would stack overflow

		#Time: O(a*b) => O(nm)
		#Space: O(n)^2
		for i in range(len(sequences)):											#Iterate through one of the sequences
			#Space: O(n)
			jresults = []

        	for j in range(len(sequences)):										#Iterate through the other sequence
            	f = align_genes													#Set to align_genes function
            	if banded:														#If banded align
                	f = banded_edit												#Set distance to 3

				#Space: O(n)
            	a = sequences[i][:align_length]									#Up to the desired align length
				#Space: O(n)
            	b = sequences[j][:align_length]									#Up to the desired align length

				#Time:O(a*b)
				#Space: O(n)^2
            	res = f(a, b)													#Pass in lists to be aligned...

            	if res[0] == float("inf"):										#For banded case where no values is in btm right of 2d array
                	jresults.append({'align_cost': 'No alignment possible',
                    				 'seqi_first100': 'None',
                    				 'seqj_first100': 'None'})
                	continue

				#Space: O(1)
            	path = res[1]													#Set the path from res

				#Space: O(n)
            	a = ['-'] + list(a)												#Reset list a with '-' at beginning
				#Space: O(n)
            	b = ['-'] + list(b)												#Reset list b with '-' at beginning

            	x = len(a) - 1													#Start our backtrace
            	y = len(b) - 1

            	seqi_100 = []													#Alignment string for first sequence
            	seqj_100 = []													#Alignment string for second sequence

				#Time Complexity: O(x+y) => O(n)
            	while x != 0 and y != 0:										#Get alignment string by backtracing through path

					next = path[(x,y)]

					#Time: O(1)
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

				#Space: O(1)
            	s = {'align_cost': res[0],
                	'seqi_first100': ''.join(seqi_100[len(seqi_100) - 100:][::-1]),	#Convert list to string, cut to first 100 chars, and then reverse
                	'seqj_first100': ''.join(seqj_100[len(seqj_100) - 100:][::-1])} #Convert list to string, cut to first 100 chars, and then reverse
            	jresults.append(s)
        	results.append(jresults)
    	return results
