* Auteur : Carpentier Antoine
* Année : BA2 sciences informatiques ULB 2013-2014
* Cours : Introduction à la bioinformatique
* Professeur : Tom Lenaerts
* Assistants : Catharina Olsen et Elisa Cilia

Projet 4 : Prédiction de la structure secondaire d'une séquence

Utilisation de Predictor.py

	python3 Predictor.py [-h] [-c] [-i INFO] [-d DIRECTORY] info_test directory_test

	Arguments obligatoires:
	  info_test             Fichier CATH contenant la liste des séquences à prédire

	  directory_test        Répertoire contenant les fichiers DSSP des séquences à
	  						prédire

	Arguments optionnels:
	  -h, --help            Affiche l'aide

	  -c, --compute         Empeche l'utilisation du cache (.predictor/) et recalcule
	  						toutes les fréquences.

	  -i INFO, --info INFO  Fichier CATH contenant la liste des séquences à charger
	  						pour calculer les fréquences
	  						Par défaut : dataset/CATH_info.txt

	  -d DIRECTORY, --directory DIRECTORY
	                        Répertoire contenant les fichiers DSSP des séquences à
	                        charger pour calculer les fréquences.
	                        Par défaut : dataset/dssp/