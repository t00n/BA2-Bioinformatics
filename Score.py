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
		
	def __getitem__(self, acide): # acide is a tuple (letter from seq A, letter from seq B)
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