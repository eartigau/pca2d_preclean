# HD 10700 HARPS, la chaîne complète du 2026-10-15

Nettoyage tellurique par la méthode de la dérivée, puis pca2d sur les résidus,
puis LBL sur les deux séries. Rapport :
`outputs/HD10700_TC/0-2/HD10700_TC_0-2_4b99e3.pdf` (28 pages).

## Le résultat

Sur les **86 poses communes** aux deux séries, **aucune aberrante** (les quatre
nuits qui sortaient à plus de dix MAD dans le run précédent ont disparu) :

| série | RMS | RMS nocturne | erreur médiane |
|---|---|---|---|
| nettoyé (dérivée, modèle à vitesse zéro) | **2,17 m/s** | 2,04 | 0,38 |
| nettoyé + pca2d (0 étoile, 2 observateur) | **2,59 m/s** | 2,55 | 0,89 |

pca2d **dégrade** d'un facteur 1,19, et l'erreur médiane plus que double.

## Pourquoi, mesuré

Ce n'est **pas** le bruit des spectres : la correction n'ajoute que **1,9 %**
de bruit pixel à pixel (dispersion robuste de la différence seconde, médiane
sur 8 poses). Ce n'est **pas** le masque non plus : 17 376 raies pour la série
brute contre 17 395 pour la corrigée, soit 0,1 % d'écart. Et l'extension
d'erreur est intacte, le rapport err/√flux valant 1,331 des deux côtés.

**C'est que la correction laisse des trous.** Les spectres corrigés portent
**13,7 % de pixels NaN contre 0,07 %** dans les livrés, et LBL y perd le tiers
de ses raies :

| | raies valides | erreur par raie |
|---|---|---|
| livré | 11 134 | 145 m/s |
| corrigé | 7 860 | 340 m/s |

Ces trous ne sont pas aux bords des ordres (11 % des pertes dans 10 % des
pixels, donc uniforme) mais **dans le bleu**, là où le flux de cette étoile
est le plus faible :

| domaine | perdu | flux médian |
|---|---|---|
| 380-420 nm | **29,8 %** | 8 626 |
| 420-460 nm | 21,7 % | 27 025 |
| 460-500 nm | 11,1 % | 60 494 |
| 500-560 nm | 11,8 % | 95 490 |
| 560-620 nm | 4,0 % | 113 289 |
| 620-692 nm | **3,6 %** | 122 862 |

Les sept ordres les plus touchés sont tous sous 430 nm, à 30-37 % de perte.

## La piste pour la suite

Le seuil qui écarte ces échantillons fait son travail, mais il coûte à LBL
30 % de ses raies là où le domaine HARPS en porte beaucoup. Deux choses à
regarder :

1. **`correct.mask_share`**, à 0,50 pour HARPS : un échantillon doit être vivant
   dans la moitié des poses. Dans le bleu, où le SNR est bas, cette fraction
   est dure à atteindre.
2. **restreindre le domaine** de pca2d au rouge (au-delà de 460 nm), où la
   correction ne perd que 4 %, et laisser le bleu intact plutôt que troué :
   une raie que pca2d n'a pas touchée vaut mieux pour LBL qu'une raie à moitié
   effacée.

La comparaison doit aussi se faire **à masque commun** pour que les deux séries
mesurent les mêmes raies (`tools/mask_from_model.py` existe pour ça).

## Le balayage du nombre de composantes, refait, et ce qu'il ne dit pas

Refait le 2026-10-16 avec le seuil par pixel en SNR (les balayages précédents
tournaient avec le seuil de profondeur, qui vidait 13,7 % des pixels dont le
fond des raies : ils mesuraient quel `n_earth` reconstruit le mieux un spectre
troué).

| composantes | RMS | MAD |
|---|---|---|
| 0 (non corrigé) | 2,19 | 1,74 |
| 1 | 2,11 | 1,63 |
| 2 | 2,08 | 1,58 |
| 3 | 2,09 | 1,61 |
| 4 | 2,04 | **1,40** |
| 6 | 2,07 | 1,60 |

**Et aucun de ces écarts n'est significatif.** Le bootstrap sur les 86 poses
donne une incertitude de **0,20 à 0,24** sur chaque MAD, alors que l'écart
entre le meilleur (1,40) et le non-corrigé (1,74) vaut 0,34. Les intervalles à
68 % se recouvrent tous.

Le test apparié, qui est la bonne comparaison puisqu'il porte sur les mêmes
poses et annule le bruit commun, le confirme :

| composantes | gain de MAD | tirages donnant un gain |
|---|---|---|
| 1 | +0,105 | 55 % |
| 2 | +0,156 | 71 % |
| 3 | +0,125 | 79 % |
| 4 | +0,333 | **82 %** |
| 6 | +0,133 | 63 % |

82 % au mieux, c'est-à-dire pas une preuve. La non-monotonie de la série
(1,63 → 1,58 → 1,61 → 1,40 → 1,60) est elle-même la signature du bruit : un
vrai optimum se détacherait de ses voisins.

**Ce qu'il faudrait pour trancher** : plus de poses. À 86 poses et un MAD de
l'ordre de 1,5 m/s, l'incertitude est de 0,22, donc distinguer deux réglages
qui diffèrent de 0,2 m/s demande environ quatre fois plus d'époques. HD 69830
en a 713 dans l'archive et GJ 581 252, ce qui est une raison de plus de les
mesurer.

### Complément : le MAD et le RMS ne mesurent pas la même chose, et ici ils se contredisent

La section ci-dessus conclut, sur le MAD, qu'aucun écart n'est significatif.
C'est vrai **pour le MAD**, et c'est incomplet. Sur le **RMS**, le même test
apparié donne :

| composantes | RMS | gain | tirages donnant un gain |
|---|---|---|---|
| 1 | 2,107 | +0,079 | 91 % |
| 2 | 2,078 | +0,108 | 91 % |
| 3 | 2,086 | +0,100 | 90 % |
| 4 | 2,038 | +0,148 | 97 % |
| 6 | 2,072 | +0,114 | 95 % |
| 8 | 2,000 | +0,186 | **99 %** |

Le gain est **monotone** avec le nombre de composantes et significatif, là où
le MAD oscillait sans tendance.

**Pourquoi les deux se contredisent**, mesuré sur `n_earth = 8` en découpant
les poses par leur écart à la médiane :

| les poses | avant | après | gain |
|---|---|---|---|
| cœur (sous la médiane de \|v\|) | 0,629 | 0,874 | **−0,245** |
| 80e centile | 1,847 | 1,680 | +0,168 |
| les 20 % les plus écartés | 3,484 | 3,045 | **+0,440** |

**pca2d dégrade le cœur et redresse les extrêmes.** Le MAD ne regarde que le
cœur par construction, donc il voit la dégradation ; le RMS pèse les extrêmes,
donc il voit le gain. Les deux ont raison sur ce qu'ils mesurent.

La conséquence pratique dépend de l'usage : pour chercher une planète dans une
série déjà propre, la dégradation du cœur compte ; pour rendre utilisables des
poses prises dans de mauvaises conditions, le redressement des extrêmes compte.
