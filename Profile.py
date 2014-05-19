from Sequence import *
from Score import *
from Cluster import *
from math import log

class Profile:
	def __init__(self, sequences):
		self.cluster = Cluster(sequences)
		self.score = Score()
		pssm = [x[:] for x in [[0]*len(self.score.indexes)]*len(self.cluster[0])]
		

if __name__ == '__main__':
	sequences = Sequence.load("SH2-domain/msaresults-clustal.fasta")
	clustal = Profile(sequences)