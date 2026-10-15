# pca2d sur HARPS : ce que la correction fait, et ce qu'elle coûte

HD 10700 (HD 10700), 100 poses récupérées de DACE, 380-692 nm, 71 ordres,
7,3 ans, 71 nuits, SNR médian 88 par pose (247 à l'ordre 50). Mesuré le
2026-10-15.

## Avant toute chose : la base est saine maintenant

Toutes les mesures d'avant le 2026-10-15 étaient à jeter. La série **livrée**
sortait à 11,33 m/s de rms là où le DRS de DACE rend 1,28 sur les mêmes poses,
donc chaque comparaison de correction se jouait dans du bruit et ne disait
rien.

La cause était dans le masque de raies de LBL, et le correctif existait déjà
sur `origin/significant-mask-edges` :

| commit | ce qu'il change |
|---|---|
| `645f550` | les raies aux **extrêmums significatifs** du gabarit, et plus de coupure sur le poids (`MASK_SIGNIFICANT_EDGES`, `MASK_EDGE_NSIG = 3`) |
| `0f1c82d` | `COMPIL_MAX_PIXEL_WIDTH` 50 → **500 px**, soit 409 km/s au lieu de 41 avec les pixels HARPS |

À 41 km/s, Mg b, Na D, Ca I et Hα étaient retirées du masque. Et comme les
extrêmums du gabarit étaient dominés par le bruit (dérivée seconde positive sur
48,6 % d'entre eux, contre 51,4 % négative, un tirage au sort là où une naine
G8V doit donner une nette majorité de creux), la CCF sortait même **à l'envers**
(+16,20 km/s pour une étoile à −16,6).

Après : le masque passe de 25071 raies de `line_snr` médian 1,46 à **17124
raies de `line_snr` médian 8,83**, et la série livrée tombe à **1,57 m/s**.

| | avant | après | DACE |
|---|---|---|---|
| rms | 11,33 | **1,57** | 1,28 |
| erreur médiane | 11,07 | **0,39** | 0,34 |
| `vrad` | +282,85 km/s | **−16,30 km/s** | −16,6 |

Corrélation de **0,925** avec le DRS et 0,59 m/s d'écart pose à pose : les deux
chaînes mesurent la même étoile.

## Le réglage, mesuré pas à pas

64 poses communes, rms robuste de LBL.

| série | rms | erreur | chi2 réduit | NaN |
|---|---|---|---|---|
| **spectres livrés** | **1,56** | 0,39 | 27,7 | |
| 0-7, grille livrée (BERV nul) | 1,89 | 0,93 | 5,2 | 21,7 % |
| 0-7, repère observateur, unanimité | 5,19 | 1,84 | 10,6 | 43,0 % |
| 0-7, repère observateur, 75 % | 2,44 | 0,94 | 11,7 | 17,0 % |
| 0-7, repère observateur, 50 % | 2,42 | 0,92 | 9,9 | 14,2 % |
| 0-1, repère observateur, 50 % | 1,91 | 0,92 | 6,9 | |
| **0-2, repère observateur, 50 %** | **1,80** | 0,92 | **5,4** | 15,0 % |
| 0-3, repère observateur, 50 % | 1,79 | 0,91 | 5,6 | |

Trois réglages, trois raisons.

**`input.to_observer`.** La logique de pca2d veut que l'ajustement travaille
dans le repère de l'observateur **et** dans celui de l'étoile, comme sur SPIRou
et NIRPS. Sur HARPS seule la grille de départ diffère (le S2D de l'ESO porte
`WAVEDATA_VAC_BARY`, déjà barycentrique, et le `get_berv` de LBL rend zéro), et
cette option la ramène à leur convention. C'est tout l'intérêt de passer par
les fonctions de lecture de LBL : tout est uniformisé, et seule la grille
d'entrée change.

**`correct.mask_share: 0.50`.** Dans le repère du sol chaque pose couvre un
intervalle décalé de son BERV, et un masque commun à l'unanimité n'en garde que
l'intersection : 39,7 % des échantillons contre 61,5 % sur la grille livrée.
Le seuil demande qu'une fraction des poses, et non toutes, aient pondéré un
échantillon. **Le jeu de raies reste UNIQUE pour la campagne** ; il est
seulement décidé à la majorité. Ce n'est pas `mask: exposure`, qui donne à
chaque pose son propre jeu et qui vaut 693 m/s.

**`twoframe.n_earth: 2`.** Entre 380 et 700 nm le bloc observateur n'a à
décrire que les bandes de l'oxygène et quelques bandes d'eau, gouvernées par la
masse d'air et la colonne d'eau : deux paramètres. La première composante porte
**29,7 %** de la variance pondérée, la deuxième **4,4 %**, la troisième 2,0 %,
les suivantes moins. Dans le H et le K de SPIRou le méthane et l'eau couvrent
des ordres entiers avec des profils qui varient autrement, d'où les sept. Le
creux est large et net : 1,91 à une composante, 1,80 à deux, 1,79 à trois,
2,42 à sept.

## Ce que la correction fait, région par région

Sur une pose, l'écart-type de `ln(corrigé / original)` :

| région | correction | NaN |
|---|---|---|
| **O2 bande B** (686-691 nm) | **0,080** | 16,7 % |
| O2 gamma (626-630) | 0,013 | 8,2 % |
| eau (592-596) | 0,0043 | 5,2 % |
| eau (646-652) | 0,0046 | 6,9 % |
| Na D (588,5-590), raie **stellaire** | 0,0049 | 10,9 % |
| **Mg b** (516-519), **témoin sans tellurique** | **0,00091** | 14,4 % |

La correction est **88 fois plus forte dans la bande B de l'oxygène que sur le
témoin de magnésium**. Elle travaille là où il y a quelque chose à corriger et
laisse le reste tranquille, ce qui est exactement ce qu'on lui demande.

Les fenêtres du rapport visent d'ailleurs ces régions : choisies le 2026-10-15
en évaluant le modèle TAPAS sur une pose et en classant les tranches de 2 nm
par transmission au 2e centile, et non d'après des noms de bandes.

## Où cela laisse pca2d sur HARPS

**1,80 m/s corrigé contre 1,56 m/s livré**, soit 15 % de plus. Il n'y a donc
toujours pas de gain, mais l'écart est passé d'un facteur 33 (51,58 contre
1,56, le nominal d'avant sur la base bruitée) à 15 %.

Deux remarques pour la suite :

1. **Le chi2 réduit dit quelque chose.** Il vaut 27,7 sur les spectres livrés
   et 5,4 sur les corrigés : la série livrée est bien plus dispersée que ses
   propres barres d'erreur ne l'autorisent, les corrigées beaucoup moins. Une
   partie des 1,56 m/s de la référence n'est donc pas du bruit photonique, et
   l'écart réel entre les deux est plus petit qu'il n'y paraît.
2. **HD 10700 n'est pas le bon banc d'essai pour une correction tellurique.**
   C'est une naine G8V extrêmement calme observée à haute altitude, dont la
   série est déjà à 1,3 m/s. Ce que pca2d vise, ce sont les cibles où le résidu
   tellurique domine, et c'est sur une M faible qu'il faut le mesurer.

## Les figures

`outputs/_harps_obs_n2/HD10700_HARPS/0-2/HD10700_HARPS_0-2_67b975.pdf`,
**27 pages**, dont **18 de river plots** : trois panneaux par fenêtre, pour les
six fenêtres telluriques du tableau ci-dessus (685-690 nm pour la bande B de
l'oxygène, 626-630 pour la gamma, 592-596 et 646-652 pour l'eau, plus Na D et
Mg b).

La figure `oh_residual` manque, et c'est normal : elle lit l'extension
`OHLine`, que seul SPIRou écrit.

**Les river plots avaient disparu du premier rapport**, et c'était ma faute :
`--mask-share`, ajouté à `reconstruct.py`, n'existait pas dans
`figures/sequence.py`, qui s'arrêtait donc sur « unrecognized arguments ». Le
rapport se construisait quand même, avec une page « What did not build » que
j'avais lue de travers en l'attribuant à `oh_residual`. `sequence.py` connaît
maintenant le seuil et applique la même règle de majorité, pour que les
panneaux cachent exactement ce que les fichiers corrigés cachent.
