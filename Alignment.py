class Alignment:
	LOCAL = 0
	GLOBAL = 1
	def __init__(self, seqA, seqB, scoreMatrix, gap_start, gap_extend, type = GLOBAL):
		self.seqA = seqA
		self.seqB = seqB
		self.scoreMatrix = scoreMatrix
		self.gap_start = gap_start
		self.gap_extend = gap_extend
		self.type = type
		self.result = []

		if (self.type == self.GLOBAL):
			start_value = float("-inf")
		else:
			start_value = 0

		# S is the "result" matrix, V holds the gap score for sequence A, W holds the gap score for sequence B
		self.S = [x[:] for x in [[start_value]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.V = [x[:] for x in [[start_value]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.W = [x[:] for x in [[start_value]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.S[0][0] = 0
		self.max = [0, 0]

		if (self.type == self.GLOBAL):
			for i in range(1, len(self.seqA)+1):
				self.S[i][0] = - self.gap_start - (i-1) * self.gap_extend

			for j in range(1, len(self.seqB)+1):
				self.S[0][j] = - self.gap_start - (j-1) * self.gap_extend

	def __repr__(self):
		ret = ""
		for i in range(0, len(self.result)):
			ret += "Result #" + str(i) + "\n" # result number
			ret += self.result[i][0] + '\n' # sequence A
			cpt = 0
			# lalign style notation ":" for identity, "." for similarity, " " for a gap
			for j in range(0, len(self.result[i][0])):
				# identiy
				if (self.result[i][0][j] == self.result[i][1][j]):
					ret += ':'
					cpt += 1
				# gap
				elif (self.result[i][0][j] == '-' or self.result[i][1][j] == '-'):
					ret += ' '
				# similarity
				else:
					ret += '.'
			ret += '\n' + self.result[i][1] + '\n' # sequence B
			ret += str(round(100*cpt/len(self.result[i][0]), 1)) + "% identity\n" # % identity
		ret += "Global score : " + str(self.S[self.max[0]][self.max[1]]) # global score
		return ret

	# fill 3 matrixes
	def computeScores(self):
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
					self.W[i][j],)	# gap in Sequence B = left

				if (self.type != self.GLOBAL):
					self.S[i][j] = max(self.S[i][j], 0)

				if (self.S[i][j] > self.S[self.max[0]][self.max[1]]):
					self.max = [i, j]
			
	# traceback the path we ran through
	def findAlignments(self, i, j, alignmentA, alignmentB):
		if (self.type == self.GLOBAL and (i > 0 or j > 0) or (self.type == self.LOCAL and self.S[i][j] > 0)):
			# alignement : diagonal
			if (i > 0 and j > 0 and self.S[i][j] == self.S[i-1][j-1] + self.scoreMatrix[self.seqA[i-1], self.seqB[j-1]]):
				self.findAlignments(i-1, j-1, self.seqA[i-1] + alignmentA, self.seqB[j-1] + alignmentB)
			# gap in sequence B : left
			elif (i > 0 and self.S[i][j] == self.V[i][j]):
				self.findAlignments(i-1, j, self.seqA[i-1] + alignmentA, "-" + alignmentB)
			# gap in sequence A : top
			elif (j > 0 and self.S[i][j] == self.W[i][j]):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j-1] + alignmentB)

		# end of backtracking : we are back in S[0][0]
		else:
			# if (self.type == self.LOCAL):
			# 	alignmentA = self.seqA[i-1] + alignmentA
			# 	alignmentB = self.seqB[j-1] + alignmentB
			self.result.append([alignmentA, alignmentB])

	def align(self):
		self.computeScores()
		if (self.type == self.GLOBAL):
			self.max = [len(self.seqA), len(self.seqB)]
		self.findAlignments(self.max[0], self.max[1], "", "")