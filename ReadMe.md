# Introduction à la bioinformatique

## Mini-projet 1 : Alignement de séquences

### Introduction

Plusieurs algorithmes existent dont le but est d'aligner deux ou plusieurs séquences de protéines et/ou de nucléotides.
Mon implémentation est basée sur l'algorithme original de Needleman et Wunsch pour l'alignement global et de Smith-Waterman pour l'alignement local qui servent à aligner deux séquences de protéines. Mes principaux ajouts sont : 
 1.  la récursivité
 2.  la pénalité affine

Pour commencer, je discuterai ces 2 points et ensuite j'aborderai le code proprement dit et j'expliquerai les classes et les méthodes les plus importantes, c'est-à-dire : 
 1.  les classes Score et Sequence
 2.  la méthode Alignment.init
 3.  la méthode Alignment.computeScores
 4.  la méthode Alignment.findAlignments

Je finirai par montrer des résultats obtenus avec mon programme et les comparerai à ceux de LALIGN

L'archive zip est disponible ici : {{:f20814:start:324440:pechealalign.zip|}}
### La récursivité

La récursivité est indispensable pour obtenir toutes les solutions optimales et pas uniquement la première. Le coût dû à la récursivité est largement remboursé par la facilité de la pile implicite et la clarté du code qui en découle. J'aurais pu utiliser une boucle mais il aurait alors fallu garder explicitement une trace des solutions en construction, ce qui aurait considérablement alourdi le code. L'inconvénient majeur est la mémoire utilisée lorsqu'on travaille sur de longue chaines (plus de 1000 protéines).
Le principe est le suivant : lors du "backtrace", on appelle la fonction récursivement sur chaque chemin qui est emprunté. On ajoute la solution à l'ensemble des solutions lorsqu'on arrive à la condition d'arrêt (en haut à gauche de la matrice dans le cas de l'algorithme global ou sur une case qui vaut 0 dans le cas de l'algorithme local).
Toutes les solutions ont le même score global puisqu'on part toujours du même endroit.

### La pénalité affine

Afin d'implémenter la pénalité affine, il est nécessaire d'utiliser deux autres matrices pour calculer le score de "gap", c'est à dire la probabilité qu'un "gap" ait été inséré. La différence avec une pénalité linéaire est que la probabilité qu'un "gap" soit inséré dépend de la probabilité qu'un "gap" ait été inséré juste avant ou non. Les deux matrices supplémentaires servent à calculer cette probabilité pour chacune des deux séquences.
Chaque case de ces matrices est calculée en prenant le score le plus haut parmi le score d'ouverture d'un nouveau "gap" si on avait un alignement précédemment et le score d'une extension de "gap" si on avait déjà un "gap" précédemment, c'est-à-dire la probabilité d'ouvrir un "gap" si on avait obtenu un alignement précédemment et celle d'étendre un "gap" si on avait obtenu un "gap" précédemment. En mathématique, cela donne, avec S la matrice des scores d'alignement, V la matrice des scores de "gap" dans la séquence A, I la pénalité d'ouverture d'un "gap" et E la pénalité d'extension d'un "gap" 

$$V_{i,j} = max(V_{i-1,j} - E, S_{i-1,j} - I - E)$$

Le score de la matrice des scores d'alignement vaut alors le maximum entre le score d'alignement, le score d'un "gap" dans la séquence A et le score d'un "gap" dans la séquence B (avec W la matrice des scores de gap dans la séquence B et M la matrice de substitution choisie pour calculer le score d'alignement)


$$ S_{i,j} = max(S_{i-1,j-1 + M_{i,j}, V_{i,j}, W{i,j}) $$
### Les classes Score et Séquence

Les classes Score et Sequence permettent de charger depuis un fichier respectivement des matrices de substitution et des séquences de protéines et de conserver leurs informations. Les deux méthodes **load** (cf les fichiers **Score.py** et **Sequence.py** dans l'archive) sont triviales et ne nécessitent pas d'explication.

La classe Score stocke la matrice de substitution et également la liste des acides aminés dans l'ordre. Ceci permet de renvoyer avec **getitem** (l'opérateur []) directement le score pour une pair d'acides aminés.
Dans le code suivant, acide est une pair d'acides aminés. Aux lignes 2 et 3, on récupère dans la liste d'acides aminés les index correspondant dans la matrice de substitution.

	:::python
	def __getitem__(self, acide):
		j = self.indexes.index(acide[0])
		i = self.indexes.index(acide[1])
		return self.matrice[i][j]

### La méthode Alignment.init

Le constructeur de Alignment initialise les 3 matrices nécessaires. Le code pertinent est 

	:::python
	self.S = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
	self.V = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
	self.W = [x[:] for x in [[float("-inf")]*(len(self.seqB)+1)]*(len(self.seqA)+1)]
	self.S[0][0] = 0
	for i in range(1, len(self.seqA)+1):
		self.S[i][0] = - self.gap_start - (i-1) * self.gap_extend
		self.W[i][0] = self.S[i][0]
	for j in range(1, len(self.seqB)+1):
		self.S[0][j] = - self.gap_start - (j-1) * self.gap_extend
		self.V[0][j] = self.S[0][j]


Dans les 3 premières lignes, on initialise toute la matrice à la plus petite valeur possible, étant donné qu'on cherchera toujours le score maximum par la suite. Cela permet de ne pas avoir de résultat faussé et d'initialiser de la même manière pour un alignement global ou local. Les matrices sont plus longues d'une case que les séquences pour permettre de calculer le score en début de séquence à partir de valeurs par défaut. 
Ces dernières sont déterminées à partir de la pénalité d'ouverture et d'extension de "gap" aux lignes 6, 7, 9 et 10. Initialiser de telle manière permet de prévoir un "gap" en tout début de séquence, de la même manière que le coin supérieur gauche est initialisé à 0 pour prévoir un alignement en début de séquence.
### La méthode Alignement.computeScores

Cette méthode remplit les 3 matrices à l'aide de formules récursives. Pour commencer (ligne 2), on spécifie l'endroit dans la matrice ou se trouvera le maximum par défaut, c'est-à-dire en bas à droite. Dans le cas d'un alignement local et puisque la matrice est initialisée à $ -\infty $, on mettra cette variable "max" à jour (ligne 10 et 11) à chaque tour de boucle de manière à savoir où démarrer la prochaine étape (le "backtrace").
Ensuite, pour chaque case des 3 matrices, on met à jour les scores. Les 2 matrices de scores de "gap" sont mises à jour en prenant le maximum entre le score d'alignement précédent diminué de la pénalité d'ouverture d'un "gap" (c'est-à-dire la probabilité d'avoir un alignement précédemment et un "gap" maintenant) et le score de "gap" précédent diminué de la pénalité d'extension d'un "gap" (c'est-à-dire la probabilité d'avoir un "gap" précédemment qui se poursuit maintenant).
Avec S la matrice de scores d'alignement, V la matrice de scores de "gap" de la séquence A, W la matrice de scores de "gap" de la séquence B, I la pénalité d'ouverture d'un "gap" et E la pénalité d'extension d'un "gap" : 

$$ V_{i,j} = max(S_{i-1, j} - I - E, V_{i-1, j} - E) $$

La matrice des scores d'alignement est mise à jour en prenant le maximum entre les scores de "gap" pour chaque séquence et le score d'alignement précédent augmenté du score provenant de la matrice de substitution pour les 2 acides aminés actuels (c'est-à-dire la probabilité d'avoir un alignement précédemment et un alignement maintenant).
Avec les mêmes notations et M la matrice de substitution : 

$$ S_{i,j} = max(V_{i,j}, W_{i,j}, S_{i-1,j-1} + M_{i, j}) $$

Enfin, aux lignes 8 et 9, si on fait un alignement local, on ramène tous les résultats négatifs à 0

	:::python
	def computeScores(self):
		self.max = [len(self.seqA), len(self.seqB)]
		for i in range(1, len(self.seqA)+1):
			for j in range(1, len(self.seqB)+1):
				self.V[i][j] = max(self.S[i-1][j] - self.gap_start - self.gap_extend, self.V[i-1][j] - self.gap_extend)
				self.W[i][j] = max(self.S[i][j-1] - self.gap_start - self.gap_extend, self.W[i][j-1] - self.gap_extend)
				self.S[i][j] = max(self.S[i-1][j-1] + self.scoreMatrix[self.seqA[i-1], self.seqB[j-1]], self.V[i][j], self.W[i][j],)
				if (self.type == self.LOCAL):
					self.S[i][j] = max(self.S[i][j], 0)
					if (self.S[i][j] > self.S[self.max[0]][self.max[1]]):
						self.max = [i, j]

### La méthode Alignment.findALignments

Cette méthode récursive parcoure en sens inverse tous les chemins ayant mené au score global optimal. Elle est appelée la première fois avec i et j valant les coordonnées de la case de la matrice des scores d'alignements ayant le score global optimal. Comme expliqué dans la précédente section, celle-ci est le coin en bas à droite dans le cas d'un alignement global et le la case ayant le score maximal dans le cas d'un alignement local.
La condition d'arrêt est dans le cas d'un alignement global, l'arrivée à l'autre bout de la matrice (en coordonnées (0,0)) et dans le cas d'un alignement local, l'arrivée à une case valant 0 (puisqu'après cette case, les alignements locaux ne seront pas optimaux). Lorsque la condition d'arrêt est atteinte on a trouvé une solution et on l'ajoute à l'ensemble des solutions.
Le chemin en sens inverse s'obtient en calculant les scores de "gap" et d'alignement pour savoir si on se déplace en diagonale, à gauche ou en haut.


	:::python
	def findAlignments(self, i, j, alignmentA, alignmentB):
		if ((self.type == self.GLOBAL and (i > 0 or j > 0)) or (self.type == self.LOCAL and self.S[i][j] > 0)):
			if (self.fromDiag(i, j)):
				self.findAlignments(i-1, j-1, self.seqA[i-1] + alignmentA, self.seqB[j-1] + alignmentB)
			if (self.fromLeft(i, j)):
				self.findAlignments(i-1, j, self.seqA[i-1] + alignmentA, "-" + alignmentB)
			if (self.fromTop(i, j)):
				self.findAlignments(i, j-1, "-" + alignmentA, self.seqB[j-1] + alignmentB)
		else:
			self.result.append((alignmentA, alignmentB))


### Résultats

Je présenterai tous mes résultats avec un alignement des deux premières séquences du fichier "PDZ-sequences.fasta" qui m'a été fourni, ceci pour des questions pratiques de longueur, les séquences du fichier maguk étant beaucoup trop longues pour rentrer sur cette page.

#### Alignement global

Les résultats que j'ai obtenu pour l'alignement global sont bons en ce qui concerne le pourcentage d'identité, la longueur et le score global de l'alignement mais LALIGN crée de plus gros "gap" et en crée moins, ce qui est théoriquement optimal en terme d'évolution. En effet, celle-ci a tendance à favoriser quelques grosses insertions ou changements d'acides aminés et non beaucoup de petites.

Avec une matrice de substitution PAM120, une pénalité d'ouverture de 4 et une pénalité d'extension de 1, mon algorithme donne le même résultat que LALIGN. Ceci s'explique par la petite différence entre la pénalité d'ouverture et d'extension et les petites valeurs de celles-ci. Mon algorithme n'est alors pas pénalisé dans ses choix de "gap".

	
	Peche à l'Align
	
	EIKLI-KGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVA--KP
	.: .: .:..::::.:.::   ..  :. .:....:..::.:...:.:..::.::.::.:.:....::.:..::::....:.. .:  : 
	RI-VIHRGSTGLGFNIVGG---ED--GE-GIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTI-IAQYK-
	40.0% identity
	86.7% similarity
	13.3% gap
	Length : 90
	Global score : 135

	
	LALIGN
	
	 using matrix file: /usr/molbio/share/fasta3/pam120.mat, gap open/ext: -4/-1
	  40.0% identity in 90 aa overlap;		 Global score: 135
	
	                10        20        30        40        50         
	unknow EIKLI-KGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSV
	        : .: .:. ::::.:.::   ..  :. .:... :..::.:. .: :. ::.::.::.:
	unknow RI-VIHRGSTGLGFNIVGG---ED--GE-GIFISFILAGGPADLSGELRKGDQILSVNGV
	                10           20           30        40        50   
	
	      60        70        80         
	unknow GLEDVMHEDAVAALKNTYDVVYLKVA--KP
	       .: .. ::.:. ::::. ..: . .:  : 
	unknow DLRNASHEQAAIALKNAGQTVTI-IAQYK-
	            60        70         80  


Par contre pour des pénalités plus grandes (14 et 4), l'alignement n'est plus le même bien qu'il s'en rapproche dans la taille du gap et son emplacement.

	
	Peche a l'Align
	
	EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAKP
	.:....:..::::.:.::      .....:....:..::.:...:.:..::.::.::.:.:....::.:..::::....:.......
	RIVIHRGSTGLGFNIVGG------EDGEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYK
	36.8% identity
	93.1% similarity
	6.9% gap
	Length : 87
	Global score : 94


	
	LALIGN
	
	 using matrix file: /usr/molbio/share/fasta3/pam120.mat, gap open/ext: -14/-4
	  36.8% identity in 87 aa overlap;		 Global score: 94
	
	               10        20        30        40        50        60
	unknow EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVG
	        : . .:. ::::.:.:: ...      .:... :..::.:. .: :. ::.::.::.:.
	unknow RIVIHRGSTGLGFNIVGGEDGE------GIFISFILAGGPADLSGELRKGDQILSVNGVD
	               10        20              30        40        50    
	
	               70        80       
	unknow LEDVMHEDAVAALKNTYDVVYLKVAKP
	       : .. ::.:. ::::. ..: . .   
	unknow LRNASHEQAAIALKNAGQTVTIIAQYK
	           60        70        80 


On peut voir dans ce dernier résultat que mon algorithme accorde plus d'importance aux alignements et celui de LALIGN aux "gaps". En effet dans mon résultat, juste avant le gap, il y a 3 identités et 1 similarité, la où il n'y que 3 similarités dans le résultat de LALIGN.

De la même manière avec la matrice de substitution BLOSUM62, et des pénalités de 4 et 1, les résultats sont les mêmes (voir ci-dessous) mais avec des pénalités de 14 et 4, les résultats sont encore plus différents car les matrices BLOSUM sont plus efficaces pour repérer des évolutions dans les séquences que des similarités, à cause de leur méthode de construction (non montré)

	
	Peche a l'Align
	
	EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAK-P
	.:....:..::::.:.::   ..  :.. :....:..::.:...:.:..::.::.::.:.:....::.:..::::....:.. .:. .
	RIVIHRGSTGLGFNIVGG---ED--GEG-IFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTI-IAQYK
	38.6% identity
	90.9% similarity
	9.1% gap
	Length : 88
	Global score : 142
	


	
	LALIGN
	
	 using matrix file: /usr/molbio/share/fasta3/blosum62.mat, gap open/ext: -4/-1
	  38.6% identity in 88 aa overlap;		 Global score: 142
	
	               10        20        30        40        50        60
	unknow EIKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVG
	       .: . .:  ::::.:.::   .   :.. :... :. :: :  .:.:. ::.::.::.: 
	unknow RIVIHRGSTGLGFNIVGG---ED--GEG-IFISFILAGGPADLSGELRKGDQILSVNGVD
	               10           20           30        40        50    
	
	               70        80        
	unknow LEDVMHEDAVAALKNTYDVVYLKVAK-P
	       :... ::.:. ::::. ..: . .:.  
	unknow LRNASHEQAAIALKNAGQTVTI-IAQYK
	           60        70         80 


#### Alignement local

Mon implémentation pour l'alignement local fonctionne bien avec des longues chaines de caractères (les "maguk" par exemple) mais souffre de quelques défauts avec les plus courtes. Ceci est du au même problème que dans l'alignement global, à savoir que mon algorithme favorise les similarités au dépend des "gaps". Dans les longues chaines de caractères, les score d'alignement et de gap deviennent très différents, l'ouverture d'un gap plus rare et l'extension plus fréquente à cause des scores très grands dans les matrices, ce qui "règle" les défauts de mon implémentation.
Je ne peux pas montrer ici les résultats fonctionnels avec "maguk" pour des raisons de place mais je vous renvoie au script de test (cf ReadMe.txt dans l'archive).
Je vous montre donc ici le résultat presque correct que j'obtiens avec les séquences "pdz".

	
	Peche a l'Align
	
	IKLIKGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVGLEDVMHEDAVAALKNTYDVVYLKVAK
	:....:..::::.:.::   ..  :.. :....:..::.:...:.:..::.::.::.:.:....::.:..::::....:.. .:.
	IVIHRGSTGLGFNIVGG---ED--GEG-IFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTI-IAQ
	40.0% identity
	91.8% similarity
	8.2% gap
	Length : 85
	Global score : 148
	


	
	LALIGN
	
	 using matrix file: /usr/molbio/share/fasta3/blosum62.mat (11/-4), gap-open/ext: -4/-1 E(limit)   0.05
	
	  38.4% identity in 86 aa overlap (3-86:1-79); score:  148 E(10000):     -1
	
	               10        20        30        40        50        60
	unknow KLI--KGPKGLGFSIAGGVGNQHIPGDNSIYVTKIIEGGAAHKDGRLQIGDKILAVNSVG
	       ...  .:  ::::.:.::   .   :.. :... :. :: :  .:.:. ::.::.::.: 
	unknow RIVIHRGSTGLGFNIVGG---ED--GEG-IFISFILAGGPADLSGELRKGDQILSVNGVD
	               10           20           30        40        50    
	
	               70        80      
	unknow LEDVMHEDAVAALKNTYDVVYLKVAK
	       :... ::.:. ::::. ..: . .:.
	unknow LRNASHEQAAIALKNAGQTVTI-IAQ
	           60        70          
	     


Comme on peut le voir dans le code, le problème se situe en fin de chaîne et est dû à l'absence d'un gap par rapport à LALIGN
## Mini-projet 2 : Création de matrices BLOSUM

{{:f20814:start:324440:project2-blosum.zip|}}

BLOSUM40-SH2

	
	   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * 
	A  2  0  1  1  0  1  1  0  0 -1 -1  1  1 -1  1  2  2  0 -2  0  0  0  0  0 
	R  0  2  1  1  0  1  1  0  0 -2 -1  1  1  0  0  0  1  1  0 -1  0  0  0  0 
	N  1  1  2  2  1  2  1  1  1 -2 -1  2  0 -2  1  1  1  0 -1 -1  0  0  0  0 
	D  1  1  2  2  1  2  2  0  1 -2  0  1  0 -2  1  1  1  0 -2 -1  0  0  0  0 
	C  0  0  1  1  0  1  0  0  0  0  0  0  1  1  0  1  1  0  0  0  0  0  0  0 
	Q  1  1  2  2  1  1  2  0  1 -1  0  2  1 -2  1  1  1  1 -1 -1  0  0  0  0 
	E  1  1  1  2  0  2  2  0  0 -2 -1  2  0 -2  0  0  0  1 -1 -1  0  0  0  0 
	G  0  0  1  0  0  0  0  3 -1 -1 -4  1 -1 -1  0  0  0 -1 -2 -1  0  0  0  0 
	H  0  0  1  1  0  1  0 -1  2 -3 -2  1  0  0 -1 -1  1  0  1 -1  0  0  0  0 
	I -1 -2 -2 -2  0 -1 -2 -1 -3  2  2  0  1 -1 -2 -2 -2  0 -3  3  0  0  0  0 
	L -1 -1 -1  0  0  0 -1 -4 -2  2  2 -1  3  0 -1 -1 -1  0 -1  0  0  0  0  0 
	K  1  1  2  1  0  2  2  1  1  0 -1  1  0 -2  2  1  2  0 -1  0  0  0  0  0 
	M  1  1  0  0  1  1  0 -1  0  1  3  0  1  0  0  0  0  0 -1  0  0  0  0  0 
	F -1  0 -2 -2  1 -2 -2 -1  0 -1  0 -2  0  3 -3 -2 -1  0  3 -2  0  0  0  0 
	P  1  0  1  1  0  1  0  0 -1 -2 -1  2  0 -3  2  0 -1  0 -2  1  0  0  0  0 
	S  2  0  1  1  1  1  0  0 -1 -2 -1  1  0 -2  0  2  2  1 -2 -1  0  0  0  0 
	T  2  1  1  1  1  1  0  0  1 -2 -1  2  0 -1 -1  2  2  0 -1  0  0  0  0  0 
	W  0  1  0  0  0  1  1 -1  0  0  0  0  0  0  0  1  0  1  0 -1  0  0  0  0 
	Y -2  0 -1 -2  0 -1 -1 -2  1 -3 -1 -1 -1  3 -2 -2 -1  0  3 -1  0  0  0  0 
	V  0 -1 -1 -1  0 -1 -1 -1 -1  3  0  0  0 -2  1 -1  0 -1 -1  2  0  0  0  0 
	B  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	Z  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 

	*  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 


BLOSUM40-TKC

	
	   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * 
	A  2  0 -1 -1  0 -1  1 -1 -3  0 -3  0  0 -1  1  2  1 -3 -1 -1  0  0  0  0 
	R  0  3  1  0 -2  2  0 -1  0 -3 -2  3 -2 -2  0  0 -1 -3  0 -1  0  0  0  0 
	N -1  1  2  1  0  1  0  2  1 -1 -3  2 -1 -2  0  1 -1 -1  0  0  0  0  0  0 
	D -1  0  1  2 -2  1  0  0 -1 -2 -3  0 -1 -4  0  0 -2 -2 -1 -1  0  0  0  0 
	C  0 -2  0 -2  2  0 -2 -1 -2  0 -1  0  0 -1  0  1  1 -1  1  0  0  0  0  0 
	Q -1  2  1  1  0  2  1  1  1 -1 -1  2  1 -3  1  1 -1  0 -1 -1  0  0  0  0 
	E  1  0  0  0 -2  1  3  0 -2 -2 -3  1  0 -1  0  1 -1 -1  0 -2  0  0  0  0 
	G -1 -1  2  0 -1  1  0  2  0 -3 -1  1 -2 -3 -1  0 -2 -1 -2 -2  0  0  0  0 
	H -3  0  1 -1 -2  1 -2  0  2 -4 -1  2 -2 -1 -1 -1  0  0  1 -1  0  0  0  0 
	I  0 -3 -1 -2  0 -1 -2 -3 -4  2  2 -3  0  0 -1 -4  0 -2 -3  3  0  0  0  0 
	L -3 -2 -3 -3 -1 -1 -3 -1 -1  2  2 -2  0  0 -2 -4 -2 -3  0  1  0  0  0  0 
	K  0  3  2  0  0  2  1  1  2 -3 -2  2  0 -1 -1 -1 -2  0  1 -2  0  0  0  0 
	M  0 -2 -1 -1  0  1  0 -2 -2  0  0  0  3  0 -3 -1  1 -3  0 -1  0  0  0  0 
	F -1 -2 -2 -4 -1 -3 -1 -3 -1  0  0 -1  0  3 -2 -1 -2 -2  2  0  0  0  0  0 
	P  1  0  0  0  0  1  0 -1 -1 -1 -2 -1 -3 -2  2  0 -1 -1 -2 -2  0  0  0  0 
	S  2  0  1  0  1  1  1  0 -1 -4 -4 -1 -1 -1  0  2  2 -2 -4 -2  0  0  0  0 
	T  1 -1 -1 -2  1 -1 -1 -2  0  0 -2 -2  1 -2 -1  2  3 -1  0 -1  0  0  0  0 
	W -3 -3 -1 -2 -1  0 -1 -1  0 -2 -3  0 -3 -2 -1 -2 -1  2  0 -3  0  0  0  0 
	Y -1  0  0 -1  1 -1  0 -2  1 -3  0  1  0  2 -2 -4  0  0  3 -1  0  0  0  0 
	V -1 -1  0 -1  0 -1 -2 -2 -1  3  1 -2 -1  0 -2 -2 -1 -3 -1  3  0  0  0  0 
	B  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	Z  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 

	*  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 


BLOSUM70-SH2

	
	   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * 
	A  1  0  1  2  0  1  1  0  0  0 -1  1  1 -1  1  1  2 -1 -2  0  0  0  0  0 
	R  0  2  1  1  1  1  1  0  0 -2  0  1  1  0  0  0  1  1  0 -1  0  0  0  0 
	N  1  1  1  2  0  2  1  1  0 -2 -1  2  0 -1  2  1  1  1 -1 -1  0  0  0  0 
	D  2  1  2  2  1  2  2  0  1 -2  0  1  0 -1  1  1  1  0 -2 -1  0  0  0  0 
	C  0  1  0  1  0  1  0  0  0  0  0  0  1  1  0  1  1  1  0  0  0  0  0  0 
	Q  1  1  2  2  1  1  2  0  0 -1  0  2  1 -1  1  1  1  1 -1 -1  0  0  0  0 
	E  1  1  1  2  0  2  2  0  0 -2 -1  2  0 -1  1  1  1  1 -1 -1  0  0  0  0 
	G  0  0  1  0  0  0  0  2 -1 -1 -3  1 -1 -1  1  1  0 -1 -2 -1  0  0  0  0 
	H  0  0  0  1  0  0  0 -1  2 -2 -1  1  1  0  0  0  1  1  0 -1  0  0  0  0 
	I  0 -2 -2 -2  0 -1 -2 -1 -2  2  2  0  1 -1 -1 -2 -2  0 -2  3  0  0  0  0 
	L -1  0 -1  0  0  0 -1 -3 -1  2  2  0  3 -1 -1 -1 -1  0 -1  0  0  0  0  0 
	K  1  1  2  1  0  2  2  1  1  0  0  1  0 -1  2  1  2  0 -1  0  0  0  0  0 
	M  1  1  0  0  1  1  0 -1  1  1  3  0  1  0  0  0  0  1 -1 -1  0  0  0  0 
	F -1  0 -1 -1  1 -1 -1 -1  0 -1 -1 -1  0  3 -1 -1 -1  0  3 -1  0  0  0  0 
	P  1  0  2  1  0  1  1  1  0 -1 -1  2  0 -1  2  1  0  0 -2  0  0  0  0  0 
	S  1  0  1  1  1  1  1  1  0 -2 -1  1  0 -1  1  2  2  1 -2 -1  0  0  0  0 
	T  2  1  1  1  1  1  1  0  1 -2 -1  2  0 -1  0  2  2  0 -1  0  0  0  0  0 
	W -1  1  1  0  1  1  1 -1  1  0  0  0  1  0  0  1  0  2  0 -1  0  0  0  0 
	Y -2  0 -1 -2  0 -1 -1 -2  0 -2 -1 -1 -1  3 -2 -2 -1  0  3 -1  0  0  0  0 
	V  0 -1 -1 -1  0 -1 -1 -1 -1  3  0  0 -1 -1  0 -1  0 -1 -1  2  0  0  0  0 
	B  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	Z  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 

	*  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 


BLOSUM70-TKC

	
	   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  * 
	A  2  0  0  0  0 -1  1  0 -3 -2 -2  0  0  0  1  2  1 -3  0 -1  0  0  0  0 
	R  0  3  1  0 -3  2  0  0  0 -2 -1  3 -1 -2  0 -1 -1 -2  0 -2  0  0  0  0 
	N  0  1  2  1 -1  1  1  1  1 -3 -1  2  0 -2  0  1  0 -1  0 -3  0  0  0  0 
	D  0  0  1  2 -2  1  0 -1 -1 -2 -2  0 -1 -4  0  0 -1 -2 -1 -1  0  0  0  0 
	C  0 -3 -1 -2  2 -1 -2 -1 -2  1 -2 -1  0 -1  0  1  1 -1  0  1  0  0  0  0 
	Q -1  2  1  1 -1  2  1  1  1 -1  0  2  1 -2  1  1  0  0  0  0  0  0  0  0 
	E  1  0  1  0 -2  1  3  0  0 -2 -2  1  0 -2  0  1  0 -1  0 -1  0  0  0  0 
	G  0  0  1 -1 -1  1  0  2  1 -2 -1  1 -1 -2 -1  1 -2 -1 -2 -1  0  0  0  0 
	H -3  0  1 -1 -2  1  0  1  2 -5 -1  2 -2 -1 -1 -1  0  1  1 -3  0  0  0  0 
	I -2 -2 -3 -2  1 -1 -2 -2 -5  2  2 -2  1  0 -1 -4 -1 -2 -3  3  0  0  0  0 
	L -2 -1 -1 -2 -2  0 -2 -1 -1  2  2 -1  1  0 -2 -2 -1 -2 -1  1  0  0  0  0 
	K  0  3  2  0 -1  2  1  1  2 -2 -1  2  1 -1  0 -1 -1  0  0 -1  0  0  0  0 
	M  0 -1  0 -1  0  1  0 -1 -2  1  1  1  3  1 -3  0  1 -2  0  0  0  0  0  0 
	F  0 -2 -2 -4 -1 -2 -2 -2 -1  0  0 -1  1  3 -1 -1 -1 -1  2  1  0  0  0  0 
	P  1  0  0  0  0  1  0 -1 -1 -1 -2  0 -3 -1  2  0 -1 -2 -2 -2  0  0  0  0 
	S  2 -1  1  0  1  1  1  1 -1 -4 -2 -1  0 -1  0  2  2 -2 -3 -1  0  0  0  0 
	T  1 -1  0 -1  1  0  0 -2  0 -1 -1 -1  1 -1 -1  2  3 -1  0 -1  0  0  0  0 
	W -3 -2 -1 -2 -1  0 -1 -1  1 -2 -2  0 -2 -1 -2 -2 -1  2  0 -2  0  0  0  0 
	Y  0  0  0 -1  0  0  0 -2  1 -3 -1  0  0  2 -2 -3  0  0  3 -1  0  0  0  0 
	V -1 -2 -3 -1  1  0 -1 -1 -3  3  1 -1  0  1 -2 -1 -1 -2 -1  2  0  0  0  0 
	B  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	Z  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
	X  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 

	*  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 

## Mini-projet 3 : Création d'un PSSM

### Introduction

Une PSSM (Position Specific Scoring Matrix) permet de représenter un motif dans un groupe de séquences (de protéines dans notre cas). Elles sont très utiles pour représenter une famille ou un domaine de protéines dans son ensemble.
Ce troisième projet consiste en l'implémentation d'un algorithme de création d'une PSSM et dans la vérification de sa pertinence biologique. Le principe est semblable à la construction d'une matrice de substitution mais en calculant le score pour chaque position de chaque acide aminé dans le groupe de séquence donné. On travaille avec un alignement multiple de séquences provenant du même domaine ou de la même famille dans le but de créer une PSSM spécifique à ce domaine ou cette famille.

Pour commencer, je vous décrirai mon algorithme et mon implémentation. Puis je vous présenterai deux exemples d'utilisation avec le domaine SH2 et le domaine PTB afin de prouver le bon fonctionnement de mon programme en introduisant d'abord brièvement ces deux domaines.

L'archive zip est disponible ici : {{:f20814:start:324440:project3-pssm.zip|}}
### Algorithme

Pour l'utilisation du programme PSSM, se reporter au fichier ReadMe-PSSM.txt dans l'archive zip.

L'algorithme que j'ai utilisé est assez similaire à la création d'une matrice de substitution mais au lieu d'obtenir un score d'alignement (c'est à dire la probabilité d'alignement) pour deux acides aminés, on calcule le score d'un acide aminé à une certaine position dans la séquence (c'est à dire la probabilité qu'un certain acide aminé apparaisse à une certaine position dans une séquence).
Puisque le but est de représenter un ensemble (famille ou domaine) de protéines, on commence par aligner cet ensemble de séquences. Le calcul du score est effectué de la manière suivante : 

$$m_{ua} = log(\frac{q_{ua}}{p_{a}}$$

où $$q_{ua}$$ est la probabilité d'apparition de l'acide aminé a dans la colonne u et $$p_{a}$$ est la probabilité d'apparition de l'acide aminé a à n'importe quelle position dans une séquence.
Pour plus de précision, la quantité $$p_{a}$$ est calculée à partir d'un nombre de séquences beaucoup plus grand (dans notre cas, la base de données UniProt, cf http://web.expasy.org/docs/relnotes/relstat.html).

Le calcul de $$q_{ua}$$ est un peu particulier puisque dans le cas ou un acide aminé n'apparait pas dans une colonne, il vaudrait 0 et donc le calcul du score serait mathématiquement invalide.
On utilise donc des pseudocounts, c'est à dire qu'on introduit des données non présentes dans les séquences mais expérimentalement correctes au moyen d'une matrice de substitution.
On ajoute au calcul de la fréquence de l'acide aminé a dans la colonne u $$f_{ua}$$ le calcul de la probabilité d'apparition de cet acide aminé à partir de la probabilité qu'il soit aligné à un acide aminé b (qui correspond au score $$s_{ab}$$ de la matrice de substitution). Pour obtenir une probabilité à partir du score $$s_{ab}$$ d'une matrice de substitution, on effectue le calcul suivant :

$$s_{ab}=log(\frac{q_{ab}}{p_{a}*p_{b}}) ==> \frac{q_{ab}}{p_{a}*p_{b}} = e^{s_{ab}}$$

On obtient alors la probabilité que le résidu b soit aligné en face du résidu a.
La probabilité rectifiée par les pseudocounts vaut donc la somme des probabilités d'alignement du résidu a multipliée par la probabilité d'apparition du résidu a :

$$g_{ua}=p_{a}*\sum_{b}f_{ub}*e^{s_{ab}}$$

Il nous reste alors à calculer la probabilité d'apparition de l'acide aminé a dans la colonne u $$q_{ua}$$ : 

$$q_{ua}=\frac{α*f_{ua}+β*g_{ua}}{α+β}$$

où α et β sont deux coefficients utilisés pour pondérer la fréquence observée dans les séquences $$f_{ua}$$ et l'information venant de la matrice de substitution. Ceci permet de donner plus d'importance aux pseudocounts quand le nombre de séquences n'est pas élevé et inversément.

Mon algorithme utilise par défaut une matrice de substitution BLOSUM62, une valeur de α égale à la longueur des séquences moins un et une valeur de β égale à la racine carrée de la longueur des séquences mais ces données sont paramètrables (se reporter au fichier ReadMe-PSSM.txt dans l'archive zip).


### Domaine SH2

Le premier domaine que j'ai analysé est le domaine SH2. Ce domaine est contenu principalement dans l'oncoprotéine Src dont il tire son nom Src Homologuous 2 mais également dans d'autre protéines chargées de l'échange de messages entre les cellules.

Le domaine SH2 permet à une protéine de s'attacher à d'autre protéines par un résidu tyrosine, c'est pourquoi on le trouve souvent dans des protéines impliquées dans la transduction d'un signal par les récepteurs tyrosines kinases qui sont liées entre autres choses à la conservation d'énergie sous forme d'ATP dans notre corps. Ces protéines, dont la famille Src fait partie, jouent un rôle dans beaucoup de fonctions cellulaire en permettant d'activer ou de désactiver des enzymes en leur attachant ou en leur détachant un ATP.
On le trouve également dans des protéines impliquées dans la transcription ou dans des enzymes commet les phosphatases dans lesquelles son but est de détecter des tyrosines.

Sa structure est formée de deux α-hélices à l'extérieur et de sept ß-feuillets anti-parallèles à l'intérieur. On peut le trouver dans de nombreux organismes, notamment l'être humain chez qui il est présent dans au moins 115 protéines identifiées à ce jour.

##### Analyse des PSSM

Les PSSM que j'ai obtenues pour le domaine SH2 se trouve dans le dossier "SH2-domain" dans l'archive zip.
En premier lieu, on remarque que les deux PSSM sont très différentes car une est basée sur l'alignement Clustal et l'autre sur l'alignement MUSCLE. J'ai commencé par chercher le domaine SH2 dans le weblogo en m'aidant du logo HMM trouvé sur PFAM (cf http://pfam.xfam.org/family/PF00017#tabview=tab4).

En observant le weblogo basé sur l'alignement MUSCLE, on remarque le domaine SH2 entre les colonnes 992 et 1278. La chaine WYHG*R*A*E*L au début et la chaine F*L*L*HY à la fin sont typiques de ce domaine. 
A cet endroit dans la PSSM, on observe une augmentation locale des scores. Ils passent d'une tranche de [-4;-2] avant 992 à une tranche de [-1;1] après la position 992 et jusqu'à la position 1313. On constate quand même des suites de 5 à 50 résidus où une diminution des scores se produit.
Cette augmentation montre qu'à cet endroit certains acides aminés ont une plus grande probabilité d'apparaitre, ce qui prouve la présence d'un motif récurrent parmi l'ensemble des séquences analysées. Ce motif est le domaine SH2. La présence de diminutions à certains endroit montre l'évolution du domaine et les insertions qui se sont produites.

Le domaine SH2 apparait clairement sur le weblogo Clustal comparé au weblogo MUSCLE. On distingue en entier le motif WYHG*I*R se trouvant au début et le motif LV*HY à la fin. Le domaine commence à la position 2947 et finit en position 3093. De la même manière qu'avec la PSSM MUSCLE, on aperçoit distinctement l'augmentation locale du score dans la PSSM Clustal entre les positions 2933 et 3070.

En comparant les deux PSSM et les deux weblogo, on peut se rendre compte des différences entre les algorithmes d'alignement multiples de MUSCLE et de Clustal. Alors que MUSCLE produit un alignement plus court et donc introduit moins de "gap", il semble moins efficace pour détecter les domaines puisque le domaine SH2 est plus dispersé dans le weblogo MUSCLE (il s'étend sur 286 positions) que dans le weblogo Clustal (il s'étend sur 137 positions). L'ordre de grandeur des acides aminés dans le weblogo Clustal est également plus proche du logo HMM que celui du weblogo MUSCLE.

### Domaine PTB

Le Phosphotyrosine-binding domain (PTB) est un domaine de protéine qui permet de se lier avec une phosphotyrosine, c'est à dire le phosphate présente sur un résidu tyrosine.

On le trouve dans la protéine tensine, présente dans la membrane cellulaire, qui permet de se lier à d'autre cellules pour créer une structure plus grande (tissu, organe...). La tensine humain possède également une fonction suppresseur de tumeurs en ralentissant la division cellulaire.

On le trouve également dans les récepteurs insuline ISR1 présents dans la membrance cellulaire qui transmettent les messages insuliques extracellulaires à l'intérieur de la cellule. Il a été prouvé que les souris déficientes en ISR1 développent une petite tendance diabétique et grandissent moins que les souris saines.

La structure du domaine PTB est constitué de 7 ß-feuillets compactes et se finit par une α-hélice au C-terminus.
## Mini-projet 4 : Prédiction de la structure secondaire d'une séquence

### Introduction

La fonction d'une protéine et les intéractions entre plusieurs protéines dépendent directement de leur structure en trois dimensions qui elle-même dépend des différentes structures de segments locaux dans les protéines, appelées structure secondaire. Celle-ci est déterminée par les liens hydrogènes formés entre deux acides aminés (ou plus) de la protéine qui donnent au segment local une certaine structure. Les structures les plus courantes sont les α-helices, les ß-feuillet, les ß-tours et les bobines. 

Un des objectifs les plus difficiles et les plus importants de la bioinformatique est la prédiction de ces structures secondaires à partir d'une séquence de protéines afin de déterminer les fonctions et les liaisons possibles d'une protéine nouvellement séquencée.

L'algorithme GOR proposé par Jean Garnier et al. en 1978 analyse la probabilité qu'un acide aminé dans une séquence se trouve dans une certaine structure en fonction des acides aminés dans son entourage proche. Cet algorithme a été plusieurs fois amélioré depuis sa création (la dernière version en date est GOR V). J'ai implémenté ici la version III.

Mon programme fonctionne en trois phases. Premièrement, il s'agit de parser un ensemble de fichier DSSP afin de récupérer un ensemble de séquences associées à leur structure. Deuxièment, le programme extrait des fréquences de cet ensemble de séquences. Pour terminer, le programme tente de prédire la structure d'un ensemble de séquences à partir de ces fréquences.

Je commencerai par expliquer ces trois phases et ensuite je montrerai les résultats obtenus avec les séquences de test.

L'archive zip est disponible ici : {{:f20814:start:324440:project4-goriii.zip|}}
### Parsing des fichiers DSSP

Le parsing des fichiers DSSP est trivial mais un point mérite d'être soulevé : afin d'optimiser l'utilisation de la mémoire, le programme parcoure d'abord le fichier CATH et récupère pour chaque nom de fichier l'ensemble des identifiants de chaines à récupérer. Ceci permet de n'ouvrir chacun des fichiers DSSP qu'une seule fois.

Le programme récupère alors les chaines dans les fichiers DSSP et les écrit dans un fichier (une chaine par ligne). Pour ce faire, il parcoure le fichier ligne par ligne en utilisant les indices de colonnes pour récupérer les données (les champs de données dans un fichier DSSP ont une taille fixe).


### Calcul des fréquences

Le programme parcoure le fichier créé lors de la phase précédente et à partir de chaque association séquence/structure, enregistre deux types de fréquences : la fréquence individuelle et la fréquence locale. Celles-ci seront cruciales lors de la troisième phase.

#### Fréquence individuelle

La fréquence individuelle $$sfreq_{RS}$$ est le nombre de fois où le résidu R est impliqué dans une structure S. Elle permet de déterminer le score individuel, c'est-à-dire la probabilité que le résidu R se trouve dans une structure S.

#### Fréquence locale

La fréquence locale $$pfreq_{RSmR_{m}}$$ est le nombre de fois où le résidu R est impliqué dans une structure S et situé à une distance m du résidu R_{m}. Cela permet de déterminer le score local, c'est-à-dire la probabilité que le résidu R se trouve dans une structure S en fonction des résidus de son entourage. Dans l'algorithme GORIII, les valeurs de m sont prises dans une fenêtre [-8,8] qui a été calculée à partir de la longueur moyenne des structures secondaires des protéines connues.

### Prédiction de la structure

La troisième étape est la plus importante car elle effectue la prédiction. Le principe de base est de parcourir la séquence dont on veut prédire la structure et pour chaque résidu de calculer l'information d'appartenance à chaque structure. La structure dont l'information est la plus haute pour ce résidu sera la structure prédite.
Mon algorithme additionne deux informations pour obtenir cette information d'appartenance : l'information directionnelle et l'information de pair. Elles ont pour point commun l'utilisation de la "window". Pour chaque résidu, l'information dépendra donc des résidus qui l'entourent.
Ces informations sont calculées à partir du score individuel d'un résidu et du score local d'un résidu, eux-mêmes calculés à partir des fréquences obtenues lors de la phase précédente. Pour plus de détails, je vous renvoie à l'article de Jean Garnier et al. 1996 (http://www.ulb.ac.be/di/map/tlenaert/Home_Tom_Lenaerts/INFO-F-208_files/1996%20Garnier.pdf)

Dans les formules suivantes, $$I(S_{j};R_{j})$$ est le score individuel d'un résidu R à l'emplacement j pour une structure S et $$I(S_{j};R_{j+m}|R_{j})$$ est le score local d'un résidu R à l'emplacement j pour une structure S en tant compte d'une window [-m;+m]

#### Information directionnelle

L'information directionnelle comprend le score individuel du résidu auquel est ajouté le score individuel des résidus dans la "window" pondéré par la distance à laquelle ils se trouvent du résidu.

$$I(S_{j};R_{1}..R_{m*2}) = I(S_{j};R_{j}) + \sum{m,m`<>`0}I(S_{m};R_{m})/abs(m)$$

#### Information de pair

L'information de pair comprend le score individuel du résidu auquel est ajouté le score local des résidus dans la "window".

$$I(S_{j};R_{1}..R_{m*2}) = I(S_{j};R_{j}) + \sum{m,m`<>`0}I(S_{j};R_{j+m}|R_{j})$$
### Test du programme

Voici les résultats d'une utilisation du programme avec les DSSP de test fournis. 
Pour commencer, on constate que la moyenne d'identité des structures prédites vaut presque 50%. Lors de test sur une plus grande quantité de séquences, j'ai obtenu un score de 55%. L'algorithme GORIII originel obtient une moyenne de 58.8% mais cet algorithme ne prend en compte que trois types de structures (les tours et les bobines ne sont pas différenciés). Comme il ne doit choisir qu'entre trois types de structures, le nombre de possibilités différentes est moindre et la certitude plus grande au moment du choix d'une structure.
Ensuite on peut remarquer que le programme prédit très bien les bobines (C) avec 66% de réussites en moyenne. La deuxième structure la mieux prédite est l'hélice (H) avec 53%. Viennent ensuite les feuillets (E) avec 37% puis les tours (T) avec 24% de prédictions correctes.

	
	Predicted/Observed residues : 
	   H        E        T        C        
	H  67.53 %  11.69 %  9.09 %  11.69 %  
	E  36.36 %  32.73 %  2.73 %  28.18 %  
	T  25.93 %  03.7 %  22.22 %  48.15 %  
	C  24.24 %  4.55 %  12.12 %  59.09 %  
	Q3 : 45.28 %
	Predicted/Observed residues : 
	   H        E        T        C        
	H  76.67 %  6.67 %  8.89 %  7.78 %  
	E  47.06 %  35.29 %  00.0 %  17.65 %  
	T  33.33 %  4.76 %  38.1 %  23.81 %  
	C  22.86 %  8.57 %  11.43 %  57.14 %  
	Q3 : 58.21 %
	Predicted/Observed residues : 
	   H        E        T        C        
	H  5.13 %  28.21 %  10.26 %  56.41 %  
	E  00.0 %  60.87 %  4.35 %  34.78 %  
	T  7.14 %  23.81 %  21.43 %  47.62 %  
	C  1.92 %  19.23 %  7.69 %  71.15 %  
	Q3 : 42.46 %
	Predicted/Observed residues : 
	   H        E        T        C        
	H  30.0 %  34.0 %  06.0 %  30.0 %  
	E  6.93 %  60.4 %  4.95 %  27.72 %  
	T  8.62 %  29.31 %  15.52 %  46.55 %  
	C  3.42 %  26.5 %  6.84 %  63.25 %  
	Q3 : 48.77 %
	Predicted/Observed residues : 
	   H        E        T        C        
	H  87.04 %  00.0 %  5.56 %  7.41 %  
	E  78.12 %  00.0 %  6.25 %  15.62 %  
	T  50.0 %  00.0 %  25.0 %  25.0 %  
	C  13.33 %  00.0 %  6.67 %  80.0 %  
	Q3 : 54.87 %
	Q3 avg :  49.91683830376933



