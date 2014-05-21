from Sequence import *
from Score import *
from Cluster import *
from math import log, log10, sqrt

class Profile:
	def __init__(self, sequences):
		self.cluster = Cluster(sequences)
		self.score = Score.load("scoring-matrices/blosum62.txt")
		acides = self.score.indexes[:-4]
		Nseq = len(self.cluster)
		pssm = [x[:] for x in [[0]*len(acides)]*len(self.cluster[0])]

		p = [ 8.28, 5.53, 4.05, 5.45, 1.36, 3.94, 6.76, 7.09, 2.27, 5.99, 9.67, 5.85, 2.43, 3.86, 4.68, 6.5, 5.32, 1.07, 2.91, 6.87 ]
		for column in range(0, len(self.cluster[0])):
			for acide in range(0, len(acides)):
				alpha = len(self.cluster)
				beta = 1
				pa = p[acide]/100
				q = (alpha*self.cluster.getFrequencyInColumn(column, self.score.indexes[acide])+beta*pa)/(alpha+beta)
				pssm[column][acide] = round((log10(q/pa)), 3)
		print(acides)
		for line in pssm:
			print(line)

if __name__ == '__main__':
	# sequences = Sequence.load("SH2-domain/msaresults-clustal.fasta")
	sequences = Sequence.load("SH2-domain/msaresults-muscle.fasta")
	# sequences = [ "TGVEAENLLL", "PRAKAEEMLS", "GRKDAERQLL" ]
	clustal = Profile(sequences)