from Sequence import *
from Score import *
from Cluster import *
from math import log

class Profile:
	def __init__(self, sequences):
		self.cluster = Cluster(sequences)
		self.score = Score.load("scoring-matrices/blosum62.txt")
		pssm = [x[:] for x in [[0]*len(self.score.indexes)]*len(self.cluster[0])]

		# acides = []
		# for acide in self.score.indexes:
		# 	acides.append(self.cluster.getFrequencyOf(acide))
		for column in range(0, len(self.cluster[0])):
			for acide in range(0, len(self.score.indexes)):
				# alpha = len(self.cluster)-1
				# beta = 10
				# p = acides[acide]
				# q = (alpha*self.cluster.getFrequencyInColumn(column, self.score.indexes[acide])+beta*p)/(alpha+beta)
				# q = self.cluster.getFrequencyInColumn(column, self.score.indexes[acide])
				pssm[column][acide] = 0
				for b in self.score.indexes:
					f = self.cluster.getFrequencyInColumn(column, b)
					s = self.score[self.score.indexes[acide], b]
					N = len(self.cluster)
					pssm[column][acide] += (log(1-f)/log(1/(N+1)))*s
				# pssm[column][acide] = round(pssm[column][acide], 2)
		# for line in sequences:
		# 	print(line)
		print(self.score.indexes)
		for line in pssm:
			print(line)

if __name__ == '__main__':
	sequences = Sequence.load("SH2-domain/msaresults-clustal.fasta")
	# sequences = [ "TGVEAENLLL", "PRAKAEESLS", "GRKDAERQLL" ]
	clustal = Profile(sequences)