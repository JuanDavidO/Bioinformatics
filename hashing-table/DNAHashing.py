#!/bin/python

import sys  #check argparse ### from sys import argv
import numpy as np

def main():
	"""						
	inputF = open(sys.argv[1], 'r') ### script, inputF, k, outputF = argv
	k = int(sys.argv[2])
	outputF = open(sys.argv[3], 'w')
	"""
	inputF = open(raw_input("Ingrese nombre de archivo de entrada: "),'r')
	k = int(raw_input("Ingrese valor de k: "))
	outputF = open(raw_input("Ingrese nombre de archivo de salida: "),'w')


	lecture=inputF.read()
	#lecture="CACGCGTCGTCGATCGTTA\nCACGCGTCGTCGATCGTTA\nCACGCGTCGTCGATCGTTAC"
	if lecture[-1]=="\n":
		last=len(lecture)-k-1
	else:
		last=len(lecture)-k
	
	DNA=DNAHashTable(k)
	
	for i in xrange(last+1):   ### range instead of xrange if porting to Python 3
		kmer=lecture[i:i+k]
		if "\n" in kmer:  ### skip (index+1)%61==0 if we can assume the same text format
			if not kmer.startswith("\n"):
				skip = kmer.find("\n")
				kmer = kmer[:skip]+kmer[skip+1:]+lecture[i+k]
				DNA.updateCounts(kmer)
		else:
			DNA.updateCounts(kmer)
			
	inputF.close()
	
	for observed in getattr(DNA, 'observedKmers'):
		outputF.write("{}\t{}\n".format(observed, DNA.getCount(observed)))
	outputF.close()
	
	
	
class DNAHashTable:
	'DNA Hash Table, reads a sequence and finds k-mers present'
	def __init__(self, k):
		self.k=k
		self.countArr = np.zeros((4**k,), dtype=np.int)
		self.observedKmers = []
	def getHash(self,kmer):   ### string substitute, safe-substitute - dictionary-like object
		if len(kmer)==self.k:
			nuc="ACGT"
			i=self.k-1
			HashVal=0
			try:
				for letter in kmer:
					HashVal+=nuc.index(letter)*(4**i)
					i-=1
			except:
       				print("Oops!",sys.exc_info()[0],"occured.")
				return -1
			return HashVal
		else:
			return -1
	def updateCounts(self,kmer):
		updateIndex=self.getHash(kmer)
		if(updateIndex >= 0):
			if self.countArr[updateIndex]==0:
				self.countArr[updateIndex]=1
				self.observedKmers.append(kmer)
			else:
				self.countArr[updateIndex]+=1
		else:
			print "Unknown k-mer"
	def getCount(self, kmer):
		return self.countArr[self.getHash(kmer)]
			
main()
