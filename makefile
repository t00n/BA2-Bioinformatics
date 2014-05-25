


all: project2 project3 project4

project2:
	zip project2-blosum.zip Blosum.py Score.py Sequence.py Cluster.py blocks/*

project3:
	zip project3-pssm.zip PSSM.py Sequence.py Score.py Cluster.py SH2-domain/* PTB-domain/*

project4:
	zip project4-GORIII.zip Predictor.py DSSP.py Score.py .predictor/*

clean:
	rm *.zip