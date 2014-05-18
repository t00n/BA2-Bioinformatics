from Sequence import Sequence
from Score import Score
from Cluster import Cluster
from math import log

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
		self.computeEvolutionaryModel()

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
				cpt += cluster.getFrequencyOf(self.indexes[i])
			self.randomModel.append(cpt/len(self.clusters))


	def computeEvolutionaryModel(self):
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
						pairs.append((cluster.getFrequencyInColumn(column, proteinA), cluster.getFrequencyInColumn(column, proteinB)))

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
				self[proteinB, proteinA] = self[proteinA, proteinB]

if __name__ == '__main__':
	sequences = []
	sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109A"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109B"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109C"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109D"))
	# sequences.append(Sequence.loadFromBlocks("blocks/TKC PR00109E"))

	blosum = Blosum(40)
	blosum.addBlocks(sequences)
	blosum.computeScores()
	print(blosum)