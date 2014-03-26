class Sequence:
	def __init__(self, seq):
		self.seq = seq
		
	@staticmethod
	def load(filename):
		f = open(filename)
		sequences = []
		seq = ""
		for line in f:
			if line[0].isalpha():
				seq += line[:-1]
			else:
				sequences.append(Sequence(seq))
				seq = ""
		if (seq != ""):
			sequences.append(Sequence(seq))
		f.close()
		return sequences
	
	def __repr__(self):
		return self.seq
	
	def align(other, score):
		pass
	
class Score:
	def __init__(self, matrice, acides):
		self.matrice = matrice
		self.acides = acides
	
	@staticmethod
	def load(filename):
		f = open(filename)
		matrice = []
		for line in f:
			if line[0] == '#':
				pass
			elif line[0] == ' ':
				acides = line.split()
			else:
				l = line[1:]
				matrice.append(map(int, l.split()))
		f.close()
		return Score(matrice, acides)
		
	def __repr__(self):
		ret = ""
		for char in self.acides:
			ret += ' ' + char + ' '
		ret += '\n'
		for line in self.matrice:
			for char in line:
				if (char >= 0):
					ret += ' '
				ret += str(char) + ' '
			ret += '\n'
		return ret


sequences = Sequence.load("PDZ-sequences.fasta")
for seq in sequences:
	print(seq)

score = Score.load("blosum80.txt")
print(score)