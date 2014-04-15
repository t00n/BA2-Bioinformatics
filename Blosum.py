from Sequence import Sequence
from Score import Score

class Cluster(list):
	def __init__(self, *args):
		list.__init__(self, *args)

	def distance(self, other, C):
		result = 0
		for seq1 in self:
			for seq2 in other:
				result += seq1.identity(seq2) >= C
		return 1-(result/(len(self)*len(other)))

	def getPairNumber(self, proteinA, proteinB):
		cptA, cptB = 0, 0
		for i in range (0, len(self[0])):
			for seq in self:
				cptA += seq[i] == proteinA
				cptB += seq[i] == proteinB
		return (cptA/(len(self)), cptB/(len(self)))

class Blosum(Score):
	def __init__(self, sequences, threshold):
		Score.__init__(self)
		self.threshold = threshold
		self.clusters = []
		for seq in sequences:
			self.clusters.append(Cluster([Sequence(seq)]))

		self.makeClusters()
		self.compute()

	def __setitem__(self, acides, value): # acide is a tuple (letter from seq A, letter from seq B)
		i = self.indexes.index(acides[0]) # get index of letter acide[0]
		j = self.indexes.index(acides[1]) # get index of letter acide[1]
		self.matrix[i][j] = value

	def makeClusters(self):
		finished = False
		while (not finished):
			minVal = 1
			finished = True
			case = None
			distances = [x[:] for x in [[0]*len(self.clusters)]*len(self.clusters)]

			for i in range(0, len(self.clusters)):
				for j in range(i+1, len(self.clusters)):
					value = self.clusters[i].distance(self.clusters[j], self.threshold)
					distances[i][j] = value
					if (value <= minVal and value != 1):
						minVal = value
						case = (i, j)
					if(value == 0):
						finished = False

			if (case):
				for seq in self.clusters[case[1]]:
					self.clusters[case[0]].append(seq)

				del self.clusters[case[1]]

	def compute(self):
		self.matrix = [x[:] for x in [[float("-inf")]*len(self.indexes)]*len(self.indexes)]
		for proteinA in self.indexes:
			for proteinB in self.indexes:
				pairs = []
				for cluster in self.clusters:
					pairs.append(cluster.getPairNumber(proteinA, proteinB))
				frequencyAB = 0
				frequencyAA = 0
				frequencyBB = 0
				for i in range(0, len(pairs)):
					for j in range(0, len(pairs)):
						if (i != j):
							frequencyAB += pairs[i][0]*pairs[j][1]
					for j in range(i+1, len(pairs)):
							frequencyAA += pairs[i][0]*pairs[j][0]
							frequencyBB += pairs[i][1]*pairs[j][1]

				self[proteinA, proteinB] = round(frequencyAB/(len(self.clusters)*len(self.clusters[0][0])), 3)
				self[proteinA, proteinA] = round(frequencyAA/(len(self.clusters)*len(self.clusters[0][0])), 3)
				self[proteinB, proteinB] = round(frequencyBB/(len(self.clusters)*len(self.clusters[0][0])), 3)

if __name__ == '__main__':
	sequences = [ "TECRQ", "SSCRN", "SECEN", "ATCRN", "SDCEQ", "ASCKN", "ATCKQ" ]
	# sequences = Sequence.loadFromBlocks("blocks/TKC PR00109A")

	blosum = Blosum(sequences, 50)
	# print(blosum)