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
