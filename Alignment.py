class Alignment:
	LOCAL = "local"
	GLOBAL = "global"
	def __init__(self, seqA, seqB, scoreMatrix, gap_start, gap_extend, type = GLOBAL):
		self.seqA = seqA
		self.seqB = seqB
		self.scoreMatrix = scoreMatrix
		self.gap_start = gap_start
		self.gap_extend = gap_extend
		self.type = type

		# S is the "result" matrix, V holds the gap score for sequence B, W holds the gap score for sequence A
		self.S = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.V = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.W = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
		self.S[0][0] = 0

		for i in range(1, len(self.seqA)+1):
			self.S[i][0] = - self.gap_start - (i-1) * self.gap_extend
			self.W[i][0] = self.S[i][0]

		for j in range(1, len(self.seqB)+1):
			self.S[0][j] = - self.gap_start - (j-1) * self.gap_extend
			self.V[0][j] = self.S[0][j]

	def __repr__(self):																		# this method prints :
		ret = ""
		for i in range(0, len(self.result)):											# result number
			ret += self.result[i][0] + '\n' 												# sequence A
			identity = 0
			similarity = 0
			gap = 0
			# lalign style notation ":" for identity, "." for similarity, " " for a gap
			for j in range(0, len(self.result[i][0])):										
				if (self.result[i][0][j] == self.result[i][1][j]):							# identiy
					ret += ':'
					identity += 1
					similarity += 1															
				elif (self.result[i][0][j] == '-' or self.result[i][1][j] == '-'):			# gap
					ret += ' '
					gap += 1																
				else:																		# similarity
					ret += '.'
					similarity += 1
			ret += '\n' + self.result[i][1] + '\n' 											# sequence B
			ret += str(round(100*identity/len(self.result[i][0]), 1)) + "% identity\n" 		# % identity
			ret += str(round(100*similarity/len(self.result[i][0]), 1)) + "% similarity\n" 	# % similarity
			ret += str(round(100*gap/len(self.result[i][0]), 1)) + "% gap\n"				# % gap
			ret += "Length : " + str(len(self.result[i][0])) + "\n"
		ret += "Global score : " + str(self.S[self.max[0]][self.max[1]]) 					# global score
		return ret

	# fill 3 matrixes
	def computeScores(self):
		self.max = [len(self.seqA), len(self.seqB)]
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

				if (self.type == self.LOCAL):
					self.S[i][j] = max(self.S[i][j], 0)
					if (self.S[i][j] > self.S[self.max[0]][self.max[1]]):
						self.max = [i, j]
			
	def fromDiag(self, i, j):
		return i > 0 and j > 0 and self.S[i][j] == self.S[i-1][j-1] + self.scoreMatrix[self.seqA[i-1], self.seqB[j-1]]

	def fromTop(self, i, j):
		return j > 0 and self.S[i][j] == self.W[i][j]

	def fromLeft(self, i, j):
		return i > 0 and self.S[i][j] == self.V[i][j]

	# traceback the path we ran through
	def findAlignments(self, i, j, alignmentA, alignmentB, score):
		if ((self.type == self.GLOBAL and (i > 0 or j > 0)) or (self.type == self.LOCAL and self.S[i][j] > 0)):
			# alignement : diagonal
			if (self.fromDiag(i, j)):
				self.findAlignments(i-1, j-1, self.seqA[i-1] + alignmentA, self.seqB[j-1] + alignmentB, score + self.S[i][j])
			# gap in sequence B : left
			if (self.fromLeft(i, j)):
				self.findAlignments(i-1, j, self.seqA[i-1] + alignmentA, "-" + alignmentB, score + self.S[i][j])
			# gap in sequence A : top
			if (self.fromTop(i, j)):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j-1] + alignmentB, score + self.S[i][j])

		# end of backtracking : we are back in S[0][0]
		else:
			self.result.append((alignmentA, alignmentB))

	def align(self):
		self.computeScores()
		self.result = []
		self.maxScore = float("-inf")
		self.findAlignments(self.max[0], self.max[1], "", "", 0)