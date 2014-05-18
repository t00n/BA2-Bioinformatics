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

	def getFrequencyInColumn(self, column, acide):
		cpt = 0
		for seq in self:
			if (seq[column] == acide):
				cpt += 1
		return cpt/len(self)

	def getFrequencyOf(self, acide):
		ret = 0
		for seq in self:
			for a in seq:
				if (a == acide):
					ret += 1
		return ret/(len(self)*len(self[0]))