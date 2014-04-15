from Sequence import Sequence

class Cluster:
	def __init__(self, sequences, identity):
		self.identity = identity
		self.sequences = []

		# self.matrix = [x[:] for x in [[False]*len(sequences)]*len(sequences)]
		# for i in range(0, len(sequences)):
		# 	for j in range(0, len(sequences)):
		# 		self.matrix[i][j] = sequences[i].identity(sequences[j]) >= identity

	def __getitem__(self, index):
		return self.sequences[index]

	def __repr__(self):
		ret = ""
		for seq in self.sequences:
			ret += seq.__repr__() + "\n"
		return ret

	def match(self, i):
		pass



	def build(self, sequences, i):
		if (sequences[i] not in self.sequences):
			self.sequences.append(sequences[i])
			match = False
			index = -1
			for j in range(i+1, len(sequences)):
				print(str(self.sequences)+":"+str(match)+":"+str(i)+":"+str(j)+":"+str(index))
				if (match):
					if (sequences[index].identity(sequences[j]) >= self.identity):
						self.build(sequences, j)
				elif (sequences[i].identity(sequences[j]) >= self.identity):
					match = True
					index = j
			self.build(sequences, index)


class Blosum:
	def __init__(self, sequences):
		pass


def distance(cluster1, cluster2):
	result = 0
	for i in range(0, len(cluster1)):
		for j in range(0, len(cluster2)):
			result += cluster1[i].distance(cluster2[j])
	return result/(len(cluster1)*len(cluster2))

def global_dissim(clusters):
	result = True
	for cluster in clusters:
		tmp_result = 0
		for seq1 in cluster:
			for seq2 in cluster:
				tmp_result = max(seq1.distance(seq2), tmp_result)

		result = tmp_result <= 50 and result
	return result

if __name__ == '__main__':
	sequences = [ Sequence("ATCKQ"), Sequence("ATCRN"), Sequence("ASCKN"), Sequence("SSCRN"), Sequence("SDCEQ"), Sequence("SECEN"), Sequence("TECRQ") ]
	clusters = []

	for seq in sequences:
		clusters.append([seq])

	while (True):
		minVal = 100
		case = None
		distances = [x[:] for x in [[0]*len(clusters)]*len(clusters)]

		for i in range(0, len(clusters)):
			for j in range(i+1, len(clusters)):
				value = distance(clusters[i], clusters[j])
				distances[i][j] = value
				if (value <= minVal):
					minVal = value
					case = (i, j)

		if (case):
			for seq in clusters[case[1]]:
				clusters[case[0]].append(seq)

			del clusters[case[1]]

		for line in distances:
			print(line)

		for line in clusters:
			print(line)

		print(global_dissim(clusters))
		print("----------------------------------------")
		if (global_dissim(clusters) == False):
			break


	for line in clusters:
		print(line)