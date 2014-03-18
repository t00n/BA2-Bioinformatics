from Align import *

sequences = Sequence.load("maguk-sequences.fasta")
for seq in sequences:
	print(seq)

score = Score.load("blosum80.txt")
print(score)
