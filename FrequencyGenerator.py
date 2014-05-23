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

		for line in file:
			if (line[0] == ">"):
				seq = ""
				struct = ""
			elif (not seq):
				seq = line[:-1]
			elif (not struct):
				struct = line[:-1]
				for aa in range(0, len(seq)):
					self.FSR[self.acides.index(seq[aa])][self.structures.index(struct[aa])] += 1
		file.close()

	def IdSR(self, S , R):
		fsr = self.FSR[R][S] # number of structure j for acide i
		fnsr = sum(self.FSR[R])-fsr # number of any structure for acide minus number of structure j for that acid i
		fs = sum([x[S] for x in self.FSR]) # number of structure j
		fns = sum(sum(x) for x in self.FSR) - fs # number of any structure minus number of structure j
		return round(log10(fsr/fnsr) + log10(fns/fs), 2)

if __name__ == '__main__':
	generator = FrequencyGenerator()
	generator.loadFromFile("dataset/summary.txt")
	generator.printIdSR()