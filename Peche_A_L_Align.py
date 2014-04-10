from Score import *
from Sequence import *
from Alignment import *

pdz = Sequence.load("PDZ-sequences.fasta")
blosum62 = Score.load("blosum62.txt")
pam120 = Score.load("pam120.txt")
maguk = Sequence.load("maguk-sequences.fasta")

# for seq in sequences:
# 	print(seq)

# for seq in maguk:
# 	print(seq)

# print(score)
needlemanwunsch = Alignment(
	# "ABCDEF",
	# "ABABABA",
	# pdz[0],
	# pdz[1],
	# "ISALIGNED",
	# "THISLINE",	
	# "CARS",
	# "CATS",
	maguk[0],
	maguk[1],
	blosum62, 10, 1, Alignment.GLOBAL)

needlemanwunsch.align()
	
print(needlemanwunsch)