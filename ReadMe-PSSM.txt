* Auteur : Carpentier Antoine
* Année : BA2 sciences informatiques ULB 2013-2014
* Cours : Introduction à la bioinformatique
* Professeur : Tom Lenaerts
* Assistants : Catharina Olsen et Elisa Cilia

Projet 3 : Construction d'une PSSM (Position Specific Scoring Matrix):

Utilisation de PSSM.py : 

	python3 PSSM.py [-h] [-s SCORE] [-a ALPHA] [-b BETA] [-j JSON] alignment_file

	Argument obligatoire :
	  alignment_file        Chemin vers un fichier fasta contenant un alignement
	  						multiple à partir duquel on créera la PSSM.

	Arguments optionnels:
	  -h, --help            Montre le message d'aide.

	  -s SCORE, --score SCORE
	                        Chemin vers un fichier contenant une matrix de substitution
	                        utilisée pour le calcul des pseudocounts.
	                        La matrice utilisée par défaut est une BLOSUM62.

	  -a ALPHA, --alpha ALPHA
	                        Permet de forcer la valeur de la variable alpha lors
	                        du calcul des pseudocounts.
	                        La valeur par défaut est la longueur des séquences moins un.

	  -b BETA, --beta BETA  Permet de forcer la valeur de la variable beta lors
	                        du calcul des pseudocounts.
	                        La valeur par défaut est la racine carrée de la longueur des
	                        séquences.

	  -j JSON, --json JSON  Permet d'enregistrer la PSSM au format JSON dans un fichier.
	  						Par défaut, le programme l'affichera dans le terminal.