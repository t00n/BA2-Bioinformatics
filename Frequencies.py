from Score import Score
from math import log10
import json

class Frequencies:
	def __init__(self):
		score = Score()
		self.acides = score.indexes[:-4]
		self.structures = "HETC"

	def computeFromFile(self, filename):
		file = open(filename, "r")
		self.FSR = [x[:] for x in [[0]*len(self.structures)]*len(self.acides)]
		self.FSR2 = [
			[
				[
					[
						0 for i in range(len(self.acides))
					] 
					for j in range(17)
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
					for m in range(-8, 9):
						bb = aa+m
						if (bb >= 0 and bb < len(seq) and m != 0):
							# print(seq[aa], struct[aa], m, seq[bb])
							self.FSR2[self.acides.index(seq[aa])] \
									 [self.structures.index(struct[aa])] \
									 [m] \
									 [self.acides.index(seq[bb])] += 1
		file.close()

	def loadFromCache(self):
		f = open("dataset/FSR.txt", "r")
		self.FSR = json.load(f)
		f.close()
		f = open("dataset/FSR2.txt", "r")
		self.FSR2 = json.load(f)
		f.close()

	def saveToCache(self):
		f = open("dataset/FSR.txt", "w")
		json.dump(self.FSR, f)
		f.close()
		f = open("dataset/FSR2.txt", "w")
		json.dump(self.FSR2, f)
		f.close()

	def self_info(self, S , R):
		s = self.structures.index(S)
		r = self.acides.index(R)
		fsr = self.FSR[r][s] # number of structure j for acide i
		fnsr = sum(self.FSR[r])-fsr # number of any structure for acide minus number of structure j for that acid i
		fs = sum([x[s] for x in self.FSR]) # number of structure j
		fns = sum(sum(x) for x in self.FSR) - fs # number of any structure minus number of structure j
		return round(log10(fsr/fnsr) + log10(fns/fs), 2)

	def local_info(self, S, Rj, m, Rjm):
		aa = self.acides.index(Rj)
		s = self.structures.index(S)
		bb = self.acides.index(Rjm)
		fsrr = self.FSR2[aa][s][m][bb]
		fnsrr = sum([x[:][m][bb] for x in self.FSR2[aa]]) - fsrr
		fsr = self.FSR[aa][s]
		fnsr = sum(self.FSR[aa]) - fsr
		return round(log10(fsrr/fnsrr) + log10(fnsr/fsr), 2)

if __name__ == '__main__':
	freq = Frequencies()
	freq.loadFromCache()
	freq.computeFromFile("dataset/summary.txt")
	freq.saveToCache()
