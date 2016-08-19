#!/bin/python

import sys  #check argparse
import numpy as np

lecture = open(sys.argv[1], 'r')
k = int(sys.argv[2])
output = open(sys.argv[3], 'w')



class DNAHashTable:
	'DNA Hash Table, reads a sequence and finds k-mers present'
	def __init__(self k):
		self.k=k
		self.countArr = np.zeros((4**k,), dtype=np.int)
		self.observedKmers = []
	def calculaFuerza(self,dt):
		
