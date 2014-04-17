from Sequence import Sequence
from Score import Score
from math import log

class Cluster(list):
	def __init__(self, *args):
		list.__init__(self, *args)

	def distance(self, other, C):
		result = 0
		for seq1 in self:
			for seq2 in other:
				if (seq1.identity(seq2) >= C):
					result += 1
		return 1-(result/(len(self)*len(other)))

	def getProbabilityRatioOf(self, column, proteinA, proteinB):
		cptA, cptB = 0, 0
		for seq in self:
			if (seq[column] == proteinA):
				cptA += 1
			if (seq[column] == proteinB):
				cptB += 1
		return (cptA/(len(self)), cptB/(len(self)))

	def getProbabilityOf(self, acide):
		ret = 0
		for seq in self:
			for a in seq:
				if (a == acide):
					ret += 1
		return ret/(len(self)*len(self[0]))


class Blosum(Score):
	def __init__(self, threshold):
		Score.__init__(self)
		self.threshold = threshold
		self.matrix = [x[:] for x in [[0]*len(self.indexes)]*len(self.indexes)]
		self.nbOfBlocks = 0

	def __setitem__(self, acides, value): # acide is a tuple (letter from seq A, letter from seq B)
		i = self.indexes.index(acides[0]) # get index of letter acide[0]
		j = self.indexes.index(acides[1]) # get index of letter acide[1]
		self.matrix[i][j] = value

	def addBlocks(self, blocks):
		for block in blocks:
			self.addBlock(block)

	def addBlock(self, block):
		self.nbOfBlocks += 1
		self.clusters = []
		for seq in block:
			self.clusters.append(Cluster([Sequence(seq)]))

		self.buildClusters()
		self.computeRandomModel()
		self.computeFrequencies()

	def buildClusters(self):
		finished = False
		while (not finished):
			minVal = 1
			finished = True
			case = None

			for i in range(0, len(self.clusters)):
				for j in range(i+1, len(self.clusters)):
					value = self.clusters[i].distance(self.clusters[j], self.threshold)
					if (value != 1 and value <= minVal):
						minVal = value
						case = (i, j)
						if(value == 0):
							finished = False

			self.clusters[case[0]].extend(self.clusters[case[1]])
			del self.clusters[case[1]]

	def computeRandomModel(self):
		self.randomModel = []
		for i in range(0, len(self.indexes)):
			cpt = 0
			for cluster in self.clusters:
				cpt += cluster.getProbabilityOf(self.indexes[i])
			self.randomModel.append(cpt/len(self.clusters))


	def computeFrequencies(self):
		maxNbOfPairs = len(self.clusters)*(len(self.clusters)-1)*len(self.clusters[0][0])/2
		for a in range(0, len(self.indexes)):
			proteinA = self.indexes[a]
			for b in range(a, len(self.indexes)):
				proteinB = self.indexes[b]
				frequencyAB = 0
				frequencyAA = 0
				for column in range(0, len(self.clusters[0][0])):
					pairs = []
					for cluster in self.clusters:
						pairs.append(cluster.getProbabilityRatioOf(column, proteinA, proteinB))

					for i in range(0, len(pairs)):
						for j in range(0, len(pairs)):
							if (i != j):
								frequencyAB += pairs[i][0]*pairs[j][1]
						for j in range(i+1, len(pairs)):
								frequencyAA += pairs[i][0]*pairs[j][0]

				if (proteinA != proteinB):
					self[proteinA, proteinB] += log((frequencyAB/maxNbOfPairs)/(self.randomModel[a]*self.randomModel[b]), 2) if frequencyAB/maxNbOfPairs != 0 else 0
				else:
					self[proteinA, proteinA] += log((frequencyAA/maxNbOfPairs)/(self.randomModel[a]*self.randomModel[a]), 2) if frequencyAA/maxNbOfPairs != 0 else 0

	def computeScores(self):
		for a in range(0, len(self.indexes)):
			proteinA = self.indexes[a]
			for b in range(a, len(self.indexes)):
				proteinB = self.indexes[b]
				self[proteinA, proteinB] = round(self[proteinA, proteinB]/self.nbOfBlocks)
				# self[proteinB, proteinA] = self[proteinA, proteinB]

if __name__ == '__main__':
	# sequences = [[ "TECRQ", "SSCRN", "SECEN", "ATCRN", "SDCEQ", "ASCKN", "ATCKQ" ]]
	sequences = []
	sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109A"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109B"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109C"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109D"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109E"))

	blosum = Blosum(50)
	blosum.addBlocks(sequences)
	blosum.computeScores()
	print(blosum)