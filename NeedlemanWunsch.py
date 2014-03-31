class NeedlemanWunsch:
	def __init__(self, seqA, seqB, score, gap_start, gap_extend):
		self.seqA = seqA
		self.seqB = seqB
		self.score = score
		self.gap_start = gap_start
		self.gap_extend = gap_extend
		self.alignmentA = []
		self.alignmentB = []

	def __repr__(self):
		str = ""
		for i in range(0, len(self.alignmentA)):
			str += self.alignmentA[i] + '\n' + self.alignmentB[i] + '\n' + '\n'
		return str

	def computeScores(self):
		self.matrix = [x[:] for x in [[0]*len(self.seqB)]*len(self.seqA)]
		self.GapA = [x[:] for x in [[float("-inf")]*len(self.seqB)]*len(self.seqA)]
		self.GapB = [x[:] for x in [[float("-inf")]*len(self.seqB)]*len(self.seqA)]
		self.result = [x[:] for x in [[int]*len(self.seqB)]*len(self.seqA)]

		for i in range(1, len(self.seqA)):
			self.matrix[i][0] = float("-inf")
			self.GapA[i][0] = self.gap_start + (i-1) * self.gap_extend
			self.GapB[i][1] = self.gap_start

		for j in range(1, len(self.seqB)):
			self.matrix[0][j] = float("-inf")
			self.GapA[1][j] = self.gap_start
			self.GapB[0][j] = self.gap_start + (j-1) * self.gap_extend

		for i in range(0, len(self.seqA)):
			for j in range(0, len(self.seqB)):
				self.matrix[i][j] = max(self.matrix[i-1][j-1] + self.score[self.seqA[i], self.seqB[j]],
					self.GapA[i-1][j-1] + self.score[self.seqA[i], self.seqB[j]],
					self.GapB[i-1][j-1] + self.score[self.seqA[i], self.seqB[j]])
				self.GapA[i][j] = max(self.matrix[i-1][j] + self.gap_start, self.GapA[i-1][j] + self.gap_extend)
				self.GapB[i][j] = max(self.matrix[i][j-1] + self.gap_start, self.GapB[i][j-1] + self.gap_extend)
				self.result[i][j] = max(self.matrix[i][j], self.GapA[i][j], self.GapB[i][j])

	def findAlignments(self, i, j, alignmentA, alignmentB):
		# print(str(i) + ":" + str(j) + ":" + alignmentA + ":" + alignmentB)
		if (i > 0 or j > 0):
			# diagonal : same score
			if (i > 0 and j > 0 and self.result[i][j] == self.matrix[i-1][j-1] + self.score[self.seqA[i], self.seqB[j]]):
				self.findAlignments(i-1, j-1, self.seqA[i] + alignmentA, self.seqB[j] + alignmentB)
			# diagonal : gap A
			if (i > 0 and j > 0 and self.result[i][j] == self.GapA[i-1][j-1] + self.score[self.seqA[i], self.seqB[j]]):
				self.findAlignments(i-1, j-1, self.seqA[i] + alignmentA, self.seqB[j] + alignmentB)
			# diagonal : gap B
			if (i > 0 and j > 0 and self.result[i][j] == self.GapB[i-1][j-1] + self.score[self.seqA[i], self.seqB[j]]):
				self.findAlignments(i-1, j-1, self.seqA[i] + alignmentA, self.seqB[j] + alignmentB)
			# left : gap start
			if (i > 0 and self.result[i][j] == self.matrix[i-1][j] + self.gap_start):
				self.findAlignments(i-1, j, self.seqA[i] + alignmentA, "-" + alignmentB)
			# left : gap extend
			if (i > 0 and self.result[i][j] == self.GapA[i-1][j] + self.gap_extend):
				self.findAlignments(i-1, j, self.seqA[i] + alignmentA, "-" + alignmentB)
			# right : gap start
			if (j > 0 and self.result[i][j] == self.matrix[i][j-1] + self.gap_start):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j] + alignmentB)
			# right : gap extend
			if (j > 0 and self.result[i][j] == self.GapB[i][j-1] + self.gap_extend):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j] + alignmentB)
		else:
			self.alignmentA.append(self.seqA[0] + alignmentA)
			self.alignmentB.append(self.seqB[0] + alignmentB)

	def align(self):
		self.computeScores()
		self.findAlignments(len(self.seqA) - 1, len(self.seqB) - 1, "", "")