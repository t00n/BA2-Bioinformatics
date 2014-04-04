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
		ret = ""
		for i in range(0, len(self.alignmentA)):
			ret += self.alignmentA[i] + '\n' + self.alignmentB[i] + '\n' + '\n'
		ret += "Global score : " + str(self.S[len(self.seqA)][len(self.seqB)])
		return ret

	def computeScores(self):
		self.S = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.V = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.W = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.S[0][0] = 0
		self.V[0][0] = 0
		self.W[0][0] = 0

		for i in range(1, len(self.seqA)+1):
			self.S[i][0] = - self.gap_start - (i-1) * self.gap_extend
			self.W[i][0] = self.S[i][0]
			self.V[i][0] = self.S[i][0]

		for j in range(1, len(self.seqB)+1):
			self.S[0][j] = - self.gap_start - (j-1) * self.gap_extend
			self.V[0][j] = self.S[0][j]
			self.W[0][j] = self.S[0][j]

		# for line in self.S:
		# 	print(line)
		# for line in self.V:
		# 	print(line)
		# for line in self.W:
		# 	print(line)

		for i in range(1, len(self.seqA)+1):
			for j in range(1, len(self.seqB)+1):
				self.V[i][j] = max(
					self.S[i-1][j] - self.gap_start - self.gap_extend, 	# before = alignment and now = gap
					self.V[i-1][j] - self.gap_extend)	# before = gap and now = gap

				self.W[i][j] = max(
					self.S[i][j-1] - self.gap_start - self.gap_extend, 	# before = alignment and now = gap
					self.W[i][j-1] - self.gap_extend)	# before = gap and now = gap

				self.S[i][j] = max(
					self.S[i-1][j-1] + self.score[self.seqA[i-1], self.seqB[j-1]],	# alignment = diagonal
					self.V[i][j],	# gap in Sequence A = top
					self.W[i][j])	# gap in Sequence B = left

		for line in self.S:
			print(line)
		for line in self.V:
			print(line)
		for line in self.W:
			print(line)
			

	def findAlignments(self, i, j, alignmentA, alignmentB):
		print(str(i) + ":" + str(j) + ":" + alignmentA + ":" + alignmentB)
		print(str(self.S[i][j])+":"+str(self.W[i][j-1] - self.gap_extend)+":"+str(self.S[i][j-1] - self.gap_start - self.gap_extend))
		if (i > 0 or j > 0):
			if (i > 0 and j > 0 and self.S[i][j] == self.S[i-1][j-1] + self.score[self.seqA[i-1], self.seqB[j-1]]):
				self.findAlignments(i-1, j-1, self.seqA[i-1] + alignmentA, self.seqB[j-1] + alignmentB)
			
			if (i > 0 and self.S[i][j] == self.V[i][j]):
				self.findAlignments(i-1, j, self.seqA[i-1] + alignmentA, "-" + alignmentB)
			
			if (j > 0 and self.S[i][j] == self.W[i][j]):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j-1] + alignmentB)

		else:
			print(alignmentA)
			print(alignmentB)

	def align(self):
		self.computeScores()
		self.findAlignments(len(self.seqA), len(self.seqB), "", "")