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
