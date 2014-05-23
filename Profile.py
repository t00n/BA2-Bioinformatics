from Sequence import *
from Score import *
from Cluster import *
from math import log, log10, sqrt, exp

class Profile:
	def __init__(self, sequences):
		self.cluster = Cluster(sequences)
		self.score = Score.load("scoring-matrices/blosum62.txt")
		acides = self.score.indexes[:-4]
		Nseq = len(self.cluster)
		Lseq = len(self.cluster[0])
		p = {'T': 5.34, 'C': 1.37, 'H': 2.27, 'I': 5.95, 'K': 5.83, 'W': 1.09, 'E': 6.74, 'F': 3.86, 'P': 4.71, 'S': 6.57, 'V': 6.87, 'L': 9.66, 'R': 5.53, 'N': 4.05, 'Q': 3.93, 'G': 7.08, 'M': 2.41, 'D': 5.46, 'A': 8.26, 'Y': 2.92}
		pssm = [x[:] for x in [[0]*len(acides)]*Lseq]

		# err_rel = 0
		alpha = Nseq
		beta = sqrt(Nseq)
		for column in range(0, Lseq):
			for acide in acides:
				pa = p[acide]/100
				gua = 0
				for b in acides:
					fub = self.cluster.getFrequencyInColumn(column, b)
					gua += fub*exp(self.score[acide, b])
				gua = gua * pa
				fua = self.cluster.getFrequencyInColumn(column, acide)
				q = (alpha*fua+beta*gua)/(alpha+beta)
				pssm[column][acides.index(acide)] = round(log10(q/pa), 3)
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