


all: project3 project4

project3:
	zip project3.zip PSSM.py Sequence.py Score.py Cluster.py SH2-domain/* "multicopper oxidase family"/*

project4:
	zip project4.zip Predictor.py DSSP.py Score.py

clean:
	rm *.zip