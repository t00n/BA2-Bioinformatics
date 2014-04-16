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

	def getPairRatio(self, column, proteinA, proteinB):
		cptA, cptB = 0, 0
		for seq in self:
			cptA += seq[column] == proteinA
			cptB += seq[column] == proteinB
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
		self.matrix = [x[:] for x in [[0]*len(self.indexes)]*len(self.indexes)]
		for a in range(0, len(self.indexes)):
			proteinA = self.indexes[a]
			for b in range(a, len(self.indexes)):
				proteinB = self.indexes[b]
				for column in range(0, len(self.clusters[0][0])):
					pairs = []
					frequencyAB = 0
					frequencyAA = 0
					for cluster in self.clusters:
						pairs.append(cluster.getPairRatio(column, proteinA, proteinB))

					for i in range(0, len(pairs)):
						for j in range(0, len(pairs)):
							if (i != j):
								frequencyAB += pairs[i][0]*pairs[j][1]
						for j in range(i+1, len(pairs)):
								frequencyAA += pairs[i][0]*pairs[j][0]

					if (proteinA != proteinB):
						self[proteinA, proteinB] += frequencyAB/(len(self.clusters)*len(self.clusters[0][0]))
					else:
						self[proteinA, proteinA] += frequencyAA/(len(self.clusters)*len(self.clusters[0][0]))

				self[proteinA, proteinB] = round(self[proteinA, proteinB], 3)
				# self[proteinB, proteinA] = self[proteinA, proteinB]

if __name__ == '__main__':
	# sequences = [ "TECRQ", "SSCRN", "SECEN", "ATCRN", "SDCEQ", "ASCKN", "ATCKQ" ]
	sequences = Sequence.loadFromBlocks("blocks/TKC PR00109A")

	blosum = Blosum(sequences, 50)
	print(blosum)