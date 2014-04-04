from Score import *
from Sequence import *
from NeedlemanWunsch import *

pdz = Sequence.load("PDZ-sequences.fasta")
score = Score.load("blosum62.txt")
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
	"THISLINE",	
	"ISALIGNED",
	score, 4, 1)

needlemanwunsch.align()
	
print(needlemanwunsch)