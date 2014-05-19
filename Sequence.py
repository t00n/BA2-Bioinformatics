import re

class Sequence(str):
	def __new__(cls, *args, **kw):
		return str.__new__(cls, *args, **kw)
		
	def load(filename):
		f = open(filename)
		sequences = []
		seq = ""
		for line in f:
			if line[0] != '>':
				seq += line[:-1]
			else:
				if (seq != ""):
					sequences.append(Sequence(seq))
				seq = ""
		if (seq != ""):
			sequences.append(Sequence(seq))
		f.close()
		return sequences

	def loadFromBlocks(filename):
		f = open(filename)
		sequences = []
		for line in f:
			sequences.append(Sequence(line[26:].split()[0]))
		f.close()
		return sequences

	def identity(self, other):
		if (len(self) == len(other)):
			identity = 0
			for i in range(0, len(self)):
				if (self[i] == other[i]):
					identity += 1
			return identity*100/len(self)
		return 0

	def distance(self, other):
		return 100-self.identity(other)

if __name__ == '__main__':
	seqA = Sequence("ATCKQ")
	seqB = Sequence("ATCRN")
	seqC = Sequence("ASCKN")
	seqD = Sequence("ABABABC")
	assert(seqA.identity(seqB) == 60)
	assert(seqA.identity(seqC) == 60)

	assert(seqA.identity(seqD) == 0)