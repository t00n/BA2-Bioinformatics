from Frequencies import Frequencies
from Sequence import Sequence
from math import sqrt

class Predictor:
	def __init__(self, freq):
		self.frequencies = freq

	def predict(self, sequence):
		prediction = [x for x in range(len(sequence))]
		for aa in range(len(sequence)):
			maX = float("-inf")
			beg = max(0, aa-8)
			end = min(len(sequence), aa+9)
			for S in "HETC":
				tmp = self.frequencies.self_info(S, sequence[aa]) \
					+ self.frequencies.pair_info(S, sequence[beg:end], aa-beg) \
					+ self.frequencies.dir_info(S, sequence[beg:end], aa-beg )
				if (tmp > maX):
					prediction[aa] = S
					maX = tmp
		return "".join(prediction)

if __name__ == '__main__':
	freq = Frequencies()
	freq.loadFromCache()
	predictor = Predictor(freq)
	# print(predictor.pair_info("H", "RRRRRRRARRRRRRRR", 8))
	
	f = open("dataset/test.txt")
	identity = 0
	H = 0
	E = 0
	T = 0
	C = 0
	TP = 0
	TN = 0
	FP = 0
	FN = 0
	total = 0
	for line in f:
		if (line[0] == ">"):
			seq = ""
			struct = ""
		elif (not seq):
			seq = line[:-1]
		elif (not struct):
			struct = line[:-1]
			predicted = predictor.predict(seq)
			local = 0
			for i in range(len(predicted)):
				if (predicted[i] == struct[i]):
					identity += 1
					local += 1
					if (predicted[i] == "H"):
						H += 1
					if (predicted[i] == "E"):
						E += 1
					if (predicted[i] == "T"):
						T += 1
					if (predicted[i] == "C"):
						C += 1
				total += 1
			print(local*100/len(predicted))
			# print(predicted)
			# print(struct)

	Q3 = identity*100/total
	print(Q3)
	# print(H*100/total)
	# print(E*100/total)
	# print(T*100/total)
	# print(C*100/total)
	# print((TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
	f.close()