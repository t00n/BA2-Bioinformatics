from math import log10, sqrt
import json
import argparse
import os
import sys
from Score import Score
from DSSP import DSSPData

class Predictor:
	CACHE_SFREQ = "sfreq"
	CACHE_PFREQ = "pfreq"
	def __init__(self):
		score = Score()
		self.acides = score.indexes
		self.structures = "HETC"

	def computeFromFile(self, filename):
		file = open(filename, "r")
		self.sfreq = [x[:] for x in [[0]*len(self.structures)]*len(self.acides)]
		self.pfreq = [
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
					self.sfreq[self.acides.index(seq[aa])][self.structures.index(struct[aa])] += 1
					for m in range(-8, 9):
						bb = aa+m
						if (bb >= 0 and bb < len(seq) and m != 0):
							self.pfreq[self.acides.index(seq[aa])] \
									 [self.structures.index(struct[aa])] \
									 [m] \
									 [self.acides.index(seq[bb])] += 1
		file.close()

	def loadFromCache(self, cache_dir):
		with open(cache_dir + self.CACHE_SFREQ, "r") as f:
			self.sfreq = json.load(f)
		with open(cache_dir + self.CACHE_PFREQ, "r") as f:
			self.pfreq = json.load(f)

	def saveToCache(self, cache_dir):
		with open(cache_dir + self.CACHE_SFREQ, "w") as f:
			json.dump(self.sfreq, f)
		with open(cache_dir + self.CACHE_PFREQ, "w") as f:
			json.dump(self.pfreq, f)

	def self_info(self, S , R):
		s = self.structures.index(S)
		r = self.acides.index(R)
		fsr = self.sfreq[r][s] # number of structure s for acide r
		fnsr = sum(self.sfreq[r])-fsr # number of any structure for acide r minus number of structure j for that acid r
		fs = sum([self.sfreq[i][s] for i in range(len(self.acides))]) # number of structure j
		fns = sum(sum(x) for x in self.sfreq) - fs # number of any structure minus number of structure j
		return log10(fsr/fnsr) + log10(fns/fs)

	def local_info(self, S, Rj, m, Rjm):
		aa = self.acides.index(Rj)
		s = self.structures.index(S)
		bb = self.acides.index(Rjm)
		fsrr = self.pfreq[aa][s][m][bb]
		fnsrr = sum([self.pfreq[aa][i][m][bb] for i in range(len(self.structures))]) - fsrr
		fsr = self.sfreq[aa][s]
		fnsr = sum(self.sfreq[aa]) - fsr
		return log10(fsrr/fnsrr) + log10(fnsr/fsr)

	def dir_info(self, S, R, j):
		ret = self.self_info(S, R[j])
		for m in range(-8, 9):
			aa = j+m
			if (aa >= 0 and aa < len(R) and m != 0):
				ret += self.self_info(S, R[aa])/abs(m)
		return ret

	def pair_info(self, S, R, j):
		ret = self.self_info(S, R[j])
		for m in range(-8, 9):
			aa = j+m
			if (aa >= 0 and aa < len(R) and m != 0):
				ret += self.local_info(S, R[j], m, R[aa])
		return ret

	def predict(self, sequence):
		prediction = [x for x in range(len(sequence))]
		for aa in range(len(sequence)):
			maX = float("-inf")
			beg = max(0, aa-8)
			end = min(len(sequence), aa+9)
			for S in "HETC":
				tmp = self.pair_info(S, sequence[beg:end], aa-beg) \
					+ self.dir_info(S, sequence[beg:end], aa-beg )
				if (tmp > maX):
					prediction[aa] = S
					maX = tmp
		return "".join(prediction)

def convertDSSP(f, d, cache, overwrite = False):
	if (not os.path.exists(cache) or overwrite):
		print("Parsing sequences from " + f + " and " + d, end="... ")
		sys.stdout.flush()
		dssp = DSSPData()
		dssp.loadMany(f, d)
		print("Done")
		print("Saving sequences to cache...", end=" ")
		dssp.saveResult(cache)
		print("Done")


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("info_test", help="File listing the sequences to predict")
	parser.add_argument("directory_test", help="Directory containing DSSP files of the sequences to predict")
	parser.add_argument("-c", "--compute", help="Do not use cache (compute frequencies from files)", action="store_true")
	parser.add_argument("-i", "--info", type=str, help="File listing the sequences to load")
	parser.add_argument("-d", "--directory", type=str, help="Directory containing DSSP files of the sequences to load")
	args = parser.parse_args()
	predictor = Predictor()

	FILE_INFO = args.info or "dataset/CATH_info.txt"
	DIR_DSSP = args.directory or "dataset/dssp/"
	FILE_INFO_TEST = args.info_test
	DIR_DSSP_TEST = args.directory_test
	CACHE_DIR = ".predictor/"
	CACHE_SEQ = CACHE_DIR + os.path.basename(FILE_INFO)
	CACHE_SEQ_TEST = CACHE_DIR + os.path.basename(FILE_INFO_TEST)

	if (args.compute or not os.path.exists(CACHE_DIR)):
		if (not os.path.exists(CACHE_DIR)):
			print("Creating cache dir...", end=" ")
			os.mkdir(CACHE_DIR)
			print("Done")
		convertDSSP(FILE_INFO, DIR_DSSP, CACHE_SEQ, True)
		print("Computing frequencies...", end=" ")
		sys.stdout.flush()
		predictor.computeFromFile(CACHE_SEQ)
		print("Done")
		print("Saving frequencies to cache...", end=" ")
		sys.stdout.flush()
		predictor.saveToCache(CACHE_DIR)
		print("Done")
	else:
		print("Loading frequencies from cache...", end=" ")
		sys.stdout.flush()
		predictor.loadFromCache(CACHE_DIR)
		print("Done")

	convertDSSP(FILE_INFO_TEST, DIR_DSSP_TEST, CACHE_SEQ_TEST)
	
	print("Predicting structure...", end=" ")

	f = open(CACHE_SEQ_TEST)
	structures = "HETC"
	results = []
	for line in f:
		if (line[0] == ">"):
			seq = ""
			struct = ""
		elif (not seq):
			seq = line[:-1]
		elif (not struct):
			struct = line[:-1]
			predicted = predictor.predict(seq)
			Q = dict.fromkeys(structures)
			for i in Q:
				Q[i] = dict.fromkeys(structures, 0)
			for i in range(len(predicted)):
				Q[predicted[i]][struct[i]] += 1
			results.append(Q)
	f.close()
	print("Done")
	
	average = 0
	for res in results:
		identity = sum(res[a][a] for a in structures)
		total = sum(res[a][b] for a in structures for b in structures)
		Q3 = identity*100/total
		average += Q3
		print("Predicted/Observed residues : ")
		print("   ", end="")
		[print(a, end="        ") for a in structures]
		print()
		for a in structures:
			print(a, end="  ")
			t = sum(res[a][b] for b in structures)
			for b in structures:
				q = round(res[a][b]*100/t, 2)
				print("{:04}".format(q), "%", end="  ")
			print()
		print("Q3 :", round(Q3, 2), "%")
	print("Q3 avg : ", average/len(results))
