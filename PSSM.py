from Sequence import *
from Score import *
from Cluster import *
from math import log10, sqrt, exp
import json
import argparse

class PSSM:
	def __init__(self, sequences, scoring_matrix, a, b):
		cluster = Cluster(sequences)
		score = Score.load(scoring_matrix)
		acides = score.indexes[:-4]
		Nseq = len(cluster)
		Lseq = len(cluster[0])
		p = {'T': 5.34, 'C': 1.37, 'H': 2.27, 'I': 5.95, 'K': 5.83, 'W': 1.09, 'E': 6.74, 'F': 3.86, 'P': 4.71, 'S': 6.57, 'V': 6.87, 'L': 9.66, 'R': 5.53, 'N': 4.05, 'Q': 3.93, 'G': 7.08, 'M': 2.41, 'D': 5.46, 'A': 8.26, 'Y': 2.92}
		self.pssm = [x[:] for x in [[0]*len(acides)]*Lseq]

		alpha = a or Nseq-1
		beta = b or sqrt(Nseq)
		for column in range(0, Lseq):
			for acide in acides:
				pa = p[acide]/100
				gua = 0
				for b in acides:
					fub = cluster.getFrequencyInColumn(column, b)
					gua += fub*exp(score[acide, b])
				gua = gua * pa
				fua = cluster.getFrequencyInColumn(column, acide)
				q = (alpha*fua+beta*gua)/(alpha+beta)
				self.pssm[column][acides.index(acide)] = round(log10(q/pa), 2)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("alignment_file", type=str, help="Fasta file containing a multiple alignment from which to create a PSSM")
	parser.add_argument("-s", "--score", type=str, help="File containing a scoring matrix to compute score")
	parser.add_argument("-a", "--alpha", type=float, help="Alpha value to compute score")
	parser.add_argument("-b", "--beta", type=float, help="Beta value to compute score")
	parser.add_argument("-j", "--json", type=str, help="Will output PSSM in the specified JSON file")
	args = parser.parse_args()

	sequences = Sequence.load(args.alignment_file)
	clustal = PSSM(sequences, args.score or "scoring-matrices/blosum62.txt", args.alpha, args.beta)

	if (args.json):
		f = open(args.json, "w")
		json.dump(clustal.pssm, f)
		f.close()
	else:
		i = 0
		for line in clustal.pssm:
			print(i, line)
			i+=1
