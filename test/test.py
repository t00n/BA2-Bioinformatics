import json

diff = 0
total = 0
with open("toon.json", "r") as t:
	f = open("clustal.json", "r")
	toon = json.load(t)
	flo = json.load(f)
	t.close()
	f.close()
	index = 0
	for aa in flo:
		for pos in range(len(flo)):
			diff += abs(flo[aa][pos]-toon[pos][index])
			total+=1
		index+=1
	t.close()
	f.close()

print(diff/total)