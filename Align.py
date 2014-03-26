class Sequence:
	def __init__(self, seq):
		self.seq = seq
		
	def load(filename):
		f = open(filename)
		sequences = []
		seq = ""
		for line in f:
			if line[0].isalpha():
				seq += line[:-1]
			else:
				if (seq != ""):
					sequences.append(Sequence(seq))
				seq = ""
		if (seq != ""):
			sequences.append(Sequence(seq))
		f.close()
		return sequences
	
	def __getitem__(self, index):
		return self.seq[index]

	def __repr__(self):
		return self.seq

	def __len__(self):
		return len(self.seq)

	def alignEnd(self, other, score, penalty, matrix, i, j, alignmentA, alignmentB):
		while (i > 0 or j > 0):
			if (i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + score[self[j], other[i]]):
				alignmentA = self[j] + alignmentA
				alignmentB = other[i] + alignmentB
				i -= 1
				j -= 1
			elif (j > 0 and matrix[i][j] == matrix[i][j-1] + penalty):
				alignmentA = self[j] + alignmentA
				alignmentB = "-" + alignmentB
				j -= 1
			elif (i > 0 and matrix[i][j] == matrix[i-1][j] + penalty):
				alignmentA = "-" + alignmentA
				alignmentB = other[i] + alignmentB
				i -= 1
		print(self[0] + alignmentA)
		print(other[0] + alignmentB)

	def align(self, other, score, penalty):
		matrix = [x[:] for x in [[0]*len(self)]*len(other)]
		for i in range(0, len(self)):
			matrix[0][i] = penalty*i
		for i in range(0, len(other)):
			matrix[i][0] = penalty*i
		for i in range(1, len(other)):
			for j in range(1, len(self)):
				match = matrix[i-1][j-1] + score[self[j], other[i]]
				delete = matrix[i-1][j] + penalty
				insert = matrix[i][j-1] + penalty
				matrix[i][j] = max(match, delete, insert)
		alignmentA = ""
		alignmentB = ""
		i = len(other) - 1
		j = len(self) - 1
		self.alignEnd(other, score, penalty, matrix, i, j, alignmentA, alignmentB)

class Score:
	def __init__(self, matrice, indexes):
		self.matrice = matrice
		self.indexes = indexes
	
	def load(filename):
		f = open(filename)
		matrice = []
		for line in f:
			if line[0] == '#':
				pass
			elif line[0] == ' ':
				indexes = line.split()
			else:
				l = line[1:]
				matrice.append(list(map(int, l.split())))
		f.close()
		return Score(matrice, indexes)
		
	def __getitem__(self, acide):
		j = self.indexes.index(acide[0]) # get index of letter acide[0]
		i = self.indexes.index(acide[1]) # get index of letter acide[1]
		return self.matrice[i][j]

	def __repr__(self):
		ret = ""
		for char in self.indexes:
			ret += ' ' + char + ' '
		ret += '\n'
		for line in self.matrice:
			for char in line:
				if (char >= 0):
					ret += ' '
				ret += str(char) + ' '
			ret += '\n'
		return ret


pdz = Sequence.load("PDZ-sequences.fasta")
score = Score.load("blosum80.txt")
# maguk = Sequence.load("maguk-sequences.fasta")

# for seq in sequences:
# 	print(seq)

# for seq in maguk:
# 	print(seq)

# print(score)


pdz[0].align(pdz[1], score, -4)