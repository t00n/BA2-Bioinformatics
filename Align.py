from Score import *
from Sequence import *
from NeedlemanWunsch import *

pdz = Sequence.load("PDZ-sequences.fasta")
score = Score.load("pam120.txt")
# maguk = Sequence.load("maguk-sequences.fasta")

# for seq in sequences:
# 	print(seq)

# for seq in maguk:
# 	print(seq)

# print(score)
needlemanwunsch = NeedlemanWunsch(
	"EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAKP", 
	"RIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYK", 
	# "ABCDEF",
	# "ABABABA",
	# pdz[0],
	# pdz[1],
	# "CATS",	
	# "CARS",
	score, 14, 4)

needlemanwunsch.align()
	
print(needlemanwunsch)