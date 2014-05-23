from Sequence import *
from Score import *
from Cluster import *
from math import log, log10, sqrt, exp

class Profile:
	def __init__(self, sequences):
		self.cluster = Cluster(sequences)
		self.score = Score.load("scoring-matrices/blosum62.txt")
		acides = self.score.indexes[:-4]
		# Nseq = len(self.cluster)
		pssm = [x[:] for x in [[0]*len(acides)]*len(self.cluster[0])]

		p = [ 8.28, 5.53, 4.05, 5.45, 1.36, 3.94, 6.76, 7.09, 2.27, 5.99, 9.67, 5.85, 2.43, 3.86, 4.68, 6.5, 5.32, 1.07, 2.91, 6.87 ]

		# g = [x[:] for x in [[0]*len(self.cluster[0])]*len(acides)]
		# for column in range(0, len(self.cluster[0])):
		# 	for acide in range(0, len(acides)):
		# 		for b in range(0, len(acides)):
		# 			fua = self.cluster.getFrequencyInColumn(column, self.score.indexes[acide])
		# 			fub = self.cluster.getFrequencyInColumn(column, self.score.indexes[b])
		# 			qab = fua*fub/(len(self.cluster)*(len(self.cluster)-1)*len(self.cluster[0])/2)
		# 			g[acide][column] += fub*qab/p[b]
		# for line in g:
		# 	print(sum(line))

		# err_rel = 0
		alpha = len(self.cluster)
		beta = sqrt(len(self.cluster))
		for column in range(0, len(self.cluster[0])):
			for acide in range(0, len(acides)):
				pa = p[acide]/100
				gua = 0
				for b in range(0, len(acides)):
					# pb = p[b]/100
					fub = self.cluster.getFrequencyInColumn(column, self.score.indexes[b])
					# qab = pa*pb*exp(0.8*self.score.matrix[acide][b])
					gua += fub*exp(self.score.matrix[acide][b])
				gua = gua * pa
				fua = self.cluster.getFrequencyInColumn(column, self.score.indexes[acide])
				q = (alpha*fua+beta*gua)/(alpha+beta)
				pssm[column][acide] = round(log10(q/pa), 3)
				# err_rel += abs(log10(q/pa)-pssm[column][acide])
		print(acides)
		for line in pssm:
			print(line)
		# print(err_rel/(len(acides)*len(self.cluster[0])))


if __name__ == '__main__':
	# sequences = Sequence.load("SH2-domain/msaresults-clustal.fasta")
	sequences = Sequence.load("SH2-domain/msaresults-muscle.fasta")
	# sequences = [ "TGVEAENLLL", "PRAKAEEMLS", "GRKDAERQLL" ]
	clustal = Profile(sequences)