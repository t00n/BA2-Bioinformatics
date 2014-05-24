from Frequencies import Frequencies

class Predictor:
	def __init__(self, freq):
		self.frequencies = freq

	def pair_info(self, S, R):
		j = 8
		ret = self.frequencies.self_info(S, R[j])
		for m in range(-8, 9):
			if (m != 0):
				ret += self.frequencies.distance_info(S, R[j], m, R[m+j])
		return ret

if __name__ == '__main__':
	freq = Frequencies()
	freq.loadFromCache()
	predictor = Predictor(freq)
	# print(predictor.pair_info("H", "RRRRRRRRARRRRRRRR"))