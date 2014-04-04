class NeedlemanWunsch:
	def __init__(self, seqA, seqB, scoreMatrix, gap_start, gap_extend):
		self.seqA = seqA
		self.seqB = seqB
		self.scoreMatrix = scoreMatrix
		self.gap_start = gap_start
		self.gap_extend = gap_extend
		self.result = []

	def __repr__(self):
		ret = ""
		for i in range(0, len(self.result)):
			ret += "Result #" + str(i) + "\n"
			ret += self.result[i][0] + '\n'
			cpt = 0
			# lalign style notation ":" for identity, "." for similarity, " " for a gap
			for j in range(0, len(self.result[i][0])): 
				if (self.result[i][0][j] == self.result[i][1][j]):
					ret += ':'
					cpt += 1
				elif (self.result[i][0][j] == '-' or self.result[i][1][j] == '-'):
					ret += ' '
				else:
					ret += '.'
			ret += '\n' + self.result[i][1] + '\n'
			ret += str(round(100*cpt/len(self.result[i][0]), 1)) + "% identity\n" # % identity
		ret += "Global score : " + str(self.S[len(self.seqA)][len(self.seqB)]) # global score
		return ret

	def computeScores(self):
		# S is the "result" matrix, V holds the gap score for sequence A, W holds the gap score for sequence B
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

		# compute scores : gap score for sequence A then gap score for sequence B then max score (gaps or alignement)
		for i in range(1, len(self.seqA)+1):
			for j in range(1, len(self.seqB)+1):
				self.V[i][j] = max(
					self.S[i-1][j] - self.gap_start - self.gap_extend, 	# before = alignment and now = gap
					self.V[i-1][j] - self.gap_extend)	# before = gap and now = gap

				self.W[i][j] = max(
					self.S[i][j-1] - self.gap_start - self.gap_extend, 	# before = alignment and now = gap
					self.W[i][j-1] - self.gap_extend)	# before = gap and now = gap

				self.S[i][j] = max(
					self.S[i-1][j-1] + self.scoreMatrix[self.seqA[i-1], self.seqB[j-1]],	# alignment = diagonal
					self.V[i][j],	# gap in Sequence A = top
					self.W[i][j])	# gap in Sequence B = left
			
	# traceback the path we ran through
	def findAlignments(self, i, j, alignmentA, alignmentB):
		if (i > 0 or j > 0):
			if (i > 0 and j > 0 and self.S[i][j] == self.S[i-1][j-1] + self.scoreMatrix[self.seqA[i-1], self.seqB[j-1]]):
				self.findAlignments(i-1, j-1, self.seqA[i-1] + alignmentA, self.seqB[j-1] + alignmentB)
			
			if (i > 0 and self.S[i][j] == self.V[i][j]):
				self.findAlignments(i-1, j, self.seqA[i-1] + alignmentA, "-" + alignmentB)
			
			if (j > 0 and self.S[i][j] == self.W[i][j]):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j-1] + alignmentB)

		else:
			self.result.append([alignmentA, alignmentB])

	def align(self):
		self.computeScores()
		self.findAlignments(len(self.seqA), len(self.seqB), "", "")