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
		for i in range(0, min(len(self), len(other))):
			if (self[i] == other[i]):
				identity += 1
		return identity*100/max(len(self), len(other))

if __name__ == '__main__':
	seqA = Sequence("ABCD")
	seqB = Sequence("ABAB")
	seqC = Sequence("ABABAC")
	seqD = Sequence("ABCD")
	assert(seqA.identity(seqB) == 50)
	assert(seqA.identity(seqD) == 100)

	try:
		seqA.identity(seqC)
	except Exception as e:
		assert(type(e) is AssertionError)

