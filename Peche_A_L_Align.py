from Score import *
from Sequence import *
from Alignment import *

pdz = Sequence.load("PDZ-sequences.fasta")
blosum62 = Score.load("blosum62.txt")
# maguk = Sequence.load("maguk-sequences.fasta")

# for seq in sequences:
# 	print(seq)

# for seq in maguk:
# 	print(seq)

# print(score)
needlemanwunsch = Alignment(
	# "EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAKP", 
	# "RIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYK", 
	# "ABCDEF",
	# "ABABABA",
	pdz[0],
	pdz[1],
	# "ISALIGNED",	
	# "THISLINE",
	# "CARS",
	# "CATS",
	blosum62, 4, 4)

needlemanwunsch.align()
	
print(needlemanwunsch)