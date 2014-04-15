from Sequence import Sequence

def distance(cluster1, cluster2):
	result = 0
	for i in range(0, len(cluster1)):
		for j in range(0, len(cluster2)):
			result += cluster1[i].identity(cluster2[j]) >= 50 # or 38 or 60
	return 1-(result/(len(cluster1)*len(cluster2)))

if __name__ == '__main__':
	sequences = [ Sequence("TECRQ"), Sequence("SSCRN"), Sequence("SECEN"), Sequence("ATCRN"), Sequence("SDCEQ"), Sequence("ASCKN"), Sequence("ATCKQ") ]
	clusters = []

	for seq in sequences:
		clusters.append([seq])

	finished = False
	while (not finished):
		minVal = 1
		finished = True
		case = None
		distances = [x[:] for x in [[0]*len(clusters)]*len(clusters)]

		for line in clusters:
			print(line)
		for i in range(0, len(clusters)):
			for j in range(i+1, len(clusters)):
				value = distance(clusters[i], clusters[j])
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

		for line in distances:
			print(line)

		print("----------------------------------------")
		if (finished):
			break


	for line in clusters:
		print(line)