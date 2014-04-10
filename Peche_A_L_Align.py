from Score import *
from Sequence import *
from Alignment import *
import glob
import sys

#title
print ("This is Peche à l'Align, the protein sequence aligner !")

# alignement type
alignType = -1
while (not (alignType == Alignment.LOCAL or alignType == Alignment.GLOBAL)):
	alignType = str(input("Choose the alignment type (enter 'local' or 'global')\n>"))

# scoring matrix
substitutionMatrixes = []
for file in glob.glob("*.txt"):
	substitutionMatrixes.append(file)

matrix = ""
while (matrix not in substitutionMatrixes):
	for file in substitutionMatrixes:
		print(file)
	matrix = str(input("Which substitution matrix do you want to use (enter the file name) ?\n>"))

matrix = Score.load(matrix)

gap_open = -1
while(gap_open < 0):
	gap_open = int(input("Enter the opening gap penalty (>=0)\n>"))

gap_extend = -1
while(gap_extend < 0):
	gap_extend = int(input("Enter the extending gap penalty (>=0)\n>"))

# sequences
seqA = str(input("Enter the first sequence as a filename or a string embedded in quotes (ex: filename.fasta or \"THISLINE\")\n>"))

if (seqA[0] == "\""):
	seqA = seqA[1:-1]
else:
	sequences = Sequence.load(seqA)
	i = int(input("Which sequence in the file ? (0 -> " + str(len(sequences)-1) + ")\n>"))
	seqA = sequences[i]

seqB = str(input("Enter the first sequence as a filename or a string embedded in quotes (ex: filename.fasta or \"THISLINE\")\n>"))

if (seqB[0] == "\""):
	seqB = seqB[1:-1]
else:
	sequences = Sequence.load(seqB)
	i = int(input("Which sequence in the file ? (0 -> " + str(len(sequences)-1) + ")\n>"))
	seqB = sequences[i]


print("Please wait...")

needlemanwunsch = Alignment(
	seqA,
	seqB,
	# "ABCDEF",
	# "ABABABA",
	# pdz[0],
	# pdz[1],
	# "ISALIGNED",
	# "THISLINE",	
	# "CARS",
	# "CATS",
	# maguk[0],
	# maguk[1],
	matrix, gap_open, gap_extend, alignType)

needlemanwunsch.align()

print("Done.")

print(needlemanwunsch, file = sys.stderr)