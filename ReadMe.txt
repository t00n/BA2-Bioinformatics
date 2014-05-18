* Auteur : Carpentier Antoine
* Année : BA2 sciences informatiques ULB 2013-2014
* Cours : Introduction à la bioinformatique
* Professeur : Tom Lenaerts
* Assistants : Catharina Olsen et Elisa Cilia

Projet 1 : Alignements locaux et globaux de séquences

A) Test prédéfini de Peche à l'Align : l'outil d'alignement de séquences

	Dans un terminal :

	1) se placer dans le dossier racine de l'archive "zip"

	2) rendre le script de test exécutable : 

		chmod +x test.sh

	3) pour tester l'alignement global (algorithe de Needleman-Wunsch), lancer

		./test.sh global

	Les paramètres utilisés sont : 
		matrice de substitution : blosum62
		pénalité d'ouverture de gap : 5
		pénalité d'extension de gap : 2

	Toutes les séquences du fichier "PDZ-sequences.fasta" sont testées avec les autres deux par deux

	4) pour tester l'alignement local (algorithme de Smith-Waterman), lancer

		./test.sh local

	Les paramètres utilisés sont : 
		matrice de substitution : pam120
		pénalité d'ouverture de gap : 5
		pénalité d'extension de gap : 1

	Les séquences testées sont les deux premières du fichier "maguk-sequences.fasta"

	5) pour les résultats, cf section B4) et B5)

B) Utilisation de Peche à l'Align : l'outil d'alignement de séquences

	Dans un terminal : 

	1) se placer dans le dossier racine de l'archive "zip"

	2) lancer

		python3 Peche_A_L_Align.py

	3) suivre les instructions à l'écran
		Sont requis dans l'ordre : 
			* le type de l'alignement (local ou global)
			* la matrice de substitution (à choisir parmi tous les fichiers du répertoire "scores/")
			* la pénalité d'ouverture d'un gap
			* la pénalité d'extension d'un gap
			* la séquence A sous forme d'un nom de fichier ou d'une chaine de caractères entourée de guillemets (exemple: PDK-sequences.fasta ou "THISLINE")
			* la séquence B (cf séquence A)

	4) les résultats affichés sont
		* les 2 séquences alignées
		* le pourcentage d'identité
		* le pourcentage de similarité
		* le pourcentage de gap
		* la longueur des chaines
		* le score global de l'alignement

	5) pour rediriger les résultats dans un fichier, lancer

		python3 Peche_A_L_Align.py 2> fichier