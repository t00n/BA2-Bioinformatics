from Sequence import Sequence

class Cluster(list):
	def __init__(self, *args):
		list.__init__(self, *args)

	def distance(self, other, C):
		result = 0
		for seq1 in self:
			for seq2 in other:
				result += seq1.identity(seq2) >= C
		return 1-(result/(len(self)*len(other)))

if __name__ == '__main__':
	sequences = [ "TECRQ", "SSCRN", "SECEN", "ATCRN", "SDCEQ", "ASCKN", "ATCKQ" ]
	clusters = []

	for seq in sequences:
		clusters.append(Cluster([Sequence(seq)]))

	finished = False
	while (not finished):
		minVal = 1
		finished = True
		case = None
		distances = [x[:] for x in [[0]*len(clusters)]*len(clusters)]

		for i in range(0, len(clusters)):
			for j in range(i+1, len(clusters)):
				value = clusters[i].distance(clusters[j], 62)
				distances[i][j] = value
				if (value <= minVal and value != 1):
					minVal = value
					case = (i, j)
				if(value == 0):
					finished = False

		if (case):
			for seq in clusters[case[1]]:
				clusters[case[0]].append(seq)

			del clusters[case[1]]

	for cluster in clusters:
		print(cluster)