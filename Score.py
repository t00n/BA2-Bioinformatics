class Score:
	INDEXES = [ "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*" ]
	def __init__(self, matrix = [[]], indexes = INDEXES):
		self.matrix = matrix
		self.indexes = indexes
	
	def load(filename):
		f = open(filename)
		matrix = []
		indexes = Score.INDEXES
		for line in f:
			if line[0] == '#':
				pass
			elif line[0] == ' ' or line[0] == '\t':
				indexes = line.split()
			elif line[0].isalpha():
				l = line[1:]
				matrix.append(list(map(int, l.split())))
		f.close()
		assert(len(indexes) == 24)
		assert(len(matrix) == len(indexes) - 1)
		return Score(matrix, indexes)
		
	def __getitem__(self, acide): # acide is a tuple (letter from seq A, letter from seq B)
		i = self.indexes.index(acide[0]) # get index of letter acide[0]
		j = self.indexes.index(acide[1]) # get index of letter acide[1]
		return self.matrix[i][j]

	def __repr__(self):
		ret = "  "
		for char in self.indexes:
			ret += ' ' + char + ' '
		ret += '\n'
		i = 0
		for line in self.matrix:
			ret += self.indexes[i] + ' '
			i+=1
			for char in line:
				if (char >= 0): # for a good alignment
					ret += ' '
				ret += str(char) + ' '
			ret += '\n'
		return ret

if __name__ == '__main__':
	test = Score.load("scoring-matrices/blosum62.txt")
	assert(test["A","A"] == 4)
	assert(test["Q","H"] == 0)
	assert(test["X","X"] == -1)
	assert(test["X","P"] == -2)
	print(test)