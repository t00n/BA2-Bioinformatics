from Score import *
from Sequence import *
from NeedlemanWunsch import *

pdz = Sequence.load("PDZ-sequences.fasta")
score = Score.load("blosum80.txt")
# maguk = Sequence.load("maguk-sequences.fasta")

# for seq in sequences:
# 	print(seq)

# for seq in maguk:
# 	print(seq)

# print(score)
needlemanwunsch = NeedlemanWunsch(pdz[0], pdz[1], score, -4)
needlemanwunsch.align()

print(needlemanwunsch)