from Frequencies import Frequencies
from Sequence import Sequence

class Predictor:
	def __init__(self, freq):
		self.frequencies = freq

	def pair_info(self, S, R, j):
		ret = self.frequencies.self_info(S, R[j])
		for m in range(len(R)):
			if (m != j):
				ret += self.frequencies.local_info(S, R[j], m-j, R[m-j])
		return ret

	def predict(self, sequence):
		prediction = [x for x in range(len(sequence))]
		for aa in range(len(sequence)):
			info = float("-inf")
			beg = max(0, aa-8)
			end = min(len(sequence), aa+9)
			for S in "HETC":
				tmp = self.frequencies.self_info(S, sequence[aa]) \
					+ self.pair_info(S, sequence[beg:end], aa-beg)
				if (tmp > info):
					prediction[aa] = S
					info = tmp
		return "".join(prediction)

if __name__ == '__main__':
	freq = Frequencies()
	freq.loadFromCache()
	predictor = Predictor(freq)
	# print(predictor.pair_info("H", "RRRRRRRARRRRRRRR", 8))
	
	f = open("dataset/test.txt")
	average = 0
	for line in f:
		if (line[0] == ">"):
			seq = ""
			struct = ""
		elif (not seq):
			seq = line[:-1]
		elif (not struct):
			struct = line[:-1]
			predicted = Sequence(predictor.predict(seq))
			real = Sequence(struct)
			identity = real.identity(predicted)
			print(identity)
			average += identity

	average /= 5
	print(average)
	f.close()