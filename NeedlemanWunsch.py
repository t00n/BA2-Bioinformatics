class NeedlemanWunsch:
	def __init__(self, seqA, seqB, score, penalty):
		self.seqA = seqA
		self.seqB = seqB
		self.score = score
		self.penalty = penalty
		self.alignmentA = []
		self.alignmentB = []

	def __repr__(self):
		str = ""
		for i in range(0, len(self.alignmentA)):
			str += self.alignmentA[i] + '\n' + self.alignmentB[i] + '\n'
		return str

	def computeScores(self):
		self.matrix = [x[:] for x in [[0]*len(self.seqA)]*len(self.seqB)]
		for i in range(0, len(self.seqA)):
			self.matrix[0][i] = self.penalty*i
		for i in range(0, len(self.seqB)):
			self.matrix[i][0] = self.penalty*i
		for i in range(1, len(self.seqB)):
			for j in range(1, len(self.seqA)):
				match = self.matrix[i-1][j-1] + self.score[self.seqA[j], self.seqB[i]]
				delete = self.matrix[i-1][j] + self.penalty
				insert = self.matrix[i][j-1] + self.penalty
				self.matrix[i][j] = max(match, delete, insert)

	def findAlignments(self, i, j, alignmentA, alignmentB):
		if (i > 0 or j > 0):
			if (i > 0 and j > 0 and self.matrix[i][j] == self.matrix[i-1][j-1] + self.score[self.seqA[j], self.seqB[i]]):
				self.findAlignments(i-1, j-1, self.seqA[j] + alignmentA, self.seqB[i] + alignmentB)
			if (j > 0 and self.matrix[i][j] == self.matrix[i][j-1] + self.penalty):
				self.findAlignments(i, j-1, self.seqA[j] + alignmentA, "-" + alignmentB)
			if (i > 0 and self.matrix[i][j] == self.matrix[i-1][j] + self.penalty):
				self.findAlignments(i-1, j, "-" + alignmentA, self.seqB[i] + alignmentB)
		else:
			self.alignmentA.append(self.seqA[0] + alignmentA)
			self.alignmentB.append(self.seqB[0] + alignmentB)

	def align(self):
		self.computeScores()
		self.findAlignments(len(self.seqB) - 1, len(self.seqA) - 1, "", "")