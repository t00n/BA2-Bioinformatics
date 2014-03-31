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
needlemanwunsch = NeedlemanWunsch(
	# "EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAKP", 
	# "RIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYK", 
	# "ABCDEF",
	# "ABABABA",
	# pdz[0],
	# pdz[1],
	"CARTS",
	"CATS",
	score, -14, -4)
needlemanwunsch.align()

# print ("M")
# for line in needlemanwunsch.matrix:
# 	print(line)

# print("I")
# for line in needlemanwunsch.GapA:
# 	print(line)

# print("J")
# for line in needlemanwunsch.GapB:
# 	print(line)

# print("H")
# for line in needlemanwunsch.result:
# 	print(line)
for line in needlemanwunsch.result:
	print(line)
	
print(needlemanwunsch)