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

class Blosum(Score):
	def __init__(self, sequences, threshold):
		self.threshold = threshold
		self.clusters = []
		for seq in sequences:
			self.clusters.append(Cluster([Sequence(seq)]))

		self.makeClusters()
		for cluster in self.clusters:
			print(cluster)

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

if __name__ == '__main__':
	sequences = [ "TECRQ", "SSCRN", "SECEN", "ATCRN", "SDCEQ", "ASCKN", "ATCKQ" ]
	blosum = Blosum(sequences, 62)