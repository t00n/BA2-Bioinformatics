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

	def identity(self, other):
		assert(len(self) == len(other))
		identity = 0
		for i in range(0, len(self)):
			identity += self[i] == other[i]
		return identity*100/len(self)

	def distance(self, other):
		return 100-self.identity(other)

if __name__ == '__main__':
	seqA = Sequence("ATCKQ")
	seqB = Sequence("ATCRN")
	seqC = Sequence("ASCKN")
	seqD = Sequence("ABABABC")
	assert(seqA.identity(seqB) == 60)
	assert(seqA.identity(seqC) == 60)

	try:
		seqA.identity(seqD)
	except Exception as e:
		assert(type(e) is AssertionError)

