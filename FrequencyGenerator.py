from Score import Score
from math import log10

class FrequencyGenerator:
	def printIdSR(self):
		print("   ", end="")
		for i in range(0, len(self.acides)):
			print(self.acides[i] + "     ", end="")
		print()
		for i in range(0, len(self.structures)):
			print(self.structures[i]+"  ", end="")
			for j in range(0, len(self.acides)):
				print('{:05}'.format(self.IdSR(i, j)) + " ", end="")
			print()

	def loadFromFile(self, filename):
		file = open(filename, "r")
		score = Score()
		self.acides = score.indexes[:-4]
		self.structures = "HETC"
		self.FSR = [x[:] for x in [[0]*len(self.structures)]*len(self.acides)]
		self.FSR2 = [
			[
				[
					[
						0 for i in range(len(self.acides))
					] 
					for j in range(16)
				]
				for k in range(len(self.structures))
			] 
			for l in range(len(self.acides))
		]

		for line in file:
			if (line[0] == ">"):
				seq = ""
				struct = ""
			elif (not seq):
				seq = line[:-1]
			elif (not struct):
				struct = line[:-1]
				for aa in range(len(seq)):
					self.FSR[self.acides.index(seq[aa])][self.structures.index(struct[aa])] += 1
					for m in range(-8, 8):
						bb = aa+m
						if (bb >= 0 and bb < len(seq) and bb != aa):
							# print(seq[aa], struct[aa], m, seq[bb])
							self.FSR2[self.acides.index(seq[aa])] \
									 [self.structures.index(struct[aa])] \
									 [m] \
									 [self.acides.index(seq[bb])] += 1
		file.close()

	def IdSR(self, S , R):
		fsr = self.FSR[R][S] # number of structure j for acide i
		fnsr = sum(self.FSR[R])-fsr # number of any structure for acide minus number of structure j for that acid i
		fs = sum([x[S] for x in self.FSR]) # number of structure j
		fns = sum(sum(x) for x in self.FSR) - fs # number of any structure minus number of structure j
		return round(log10(fsr/fnsr) + log10(fns/fs), 2)

	def IdSjRjmRj(self, S, Rj, m, Rjm):
		aa = self.acides.index(Rj)
		s = self.structures.index(S)
		bb = self.acides.index(Rjm)
		fsrr = self.FSR2[aa][s][m][bb]
		fnsrr = sum([x[:][m][bb] for x in self.FSR2[aa]]) - fsrr
		fsr = self.FSR[aa][s]
		fnsr = sum(self.FSR[aa]) - fsr
		return round(log10(fsrr/fnsrr) + log10(fnsr/fsr), 2)

if __name__ == '__main__':
	generator = FrequencyGenerator()
	generator.loadFromFile("dataset/summary.txt")
	# generator.printIdSR()
	# print(generator.FSR2)
	print(generator.IdSjRjmRj("H","A",1,"R"))
