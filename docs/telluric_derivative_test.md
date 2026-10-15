# L'exposant H2O par la dérivée, testé sur HARPS

Demandé le 2026-10-15 : appliquer la méthode de convergence par la dérivée de
`telluric_fit` (`tellu_tools.optimize_exponents`, appelée par
`predict_abso.py`), avec l'exposant H2O libre par fichier, les espèces sèches
fixées à la masse d'air et la vitesse non libre.

Le document `docs/telluric_derivative_method.md` n'était pas dans
`papier_helium` ; la méthode a été reprise de
`papier_helium/referee/tellu_derivative.py` et de
`telluric_fit/tellu_tools.py:2413`, qui la portent tous deux.

## La méthode

On divise le spectre par la transmission modélisée et on retient l'exposant
qui rend le spectre corrigé **le plus lisse**, c'est-à-dire qui minimise la
dispersion pondérée du gradient pixel à pixel. Un exposant faux laisse à
chaque raie tellurique un résidu en forme de dérivée, qui gonfle cette
dispersion ; les raies stellaires y contribuent la même chose quel que soit
l'exposant, donc le critère ne dépend que des telluriques.

Suivant `optimize_exponents` : Nelder-Mead, bornes **larges** (0,001 à 20,
`tellu_tools_config.py:224`, contre 0,01 à 15 dans LBL), poids lisses de
`get_transmission_weights` (`depth_max` 0,6, `depth_saturated` 0,5,
`reject_saturated` 0,8, sigmoïde de 0,02) qui écartent les raies saturées et
leurs ailes, et sigmoïde d'écrêtage des aberrants à 6 sigma.

Ce qui est libre et ce qui ne l'est pas, comme demandé :

| | |
|---|---|
| H2O | **libre par fichier** : la colonne d'eau change trop vite pour un modèle par nuit |
| espèces sèches | **fixées à la masse d'air** : O2 et CO2 sont bien mélangés |
| vitesse | **non libre** : TAPAS évalué dans le repère du **sol**, pas sur la grille barycentrique |

L'ajustement est restreint aux trois bandes d'eau du domaine HARPS (586-600,
644-664, 686-696 nm), mesurées sur TAPAS : ailleurs l'eau ne dépasse pas
1 % de profondeur et n'informe rien.

## Ce que le test donne, sur 100 poses de HD 10700

| | LBL | dérivée |
|---|---|---|
| exposant H2O médian | 0,391 | 0,199 |
| poses sur la borne basse | 37 % | **35 %** |
| poses sur la borne haute | 0 | 0 |

**Les deux méthodes ne s'accordent pas** : corrélation de **0,155** entre
leurs exposants quand les deux sont libres, et chacune bloque sur un jeu de
poses différent (18 en commun, 19 débloquées par la dérivée, 17 bloquées par
elle seule).

**La borne basse n'est pas toujours un échec.** Pour la pose 35, la courbe de
coût monte de façon monotone depuis 0,001 : il n'y a vraiment pas d'eau
détectable, et 0,001 est la bonne réponse. Pour la pose 45, le minimum est net
à 0,8. Le critère distingue donc les deux cas, ce que l'extrapolation de LBL
ne faisait pas.

**Le juge, c'est le spectre.** Dispersion du gradient dans 593-595 nm, la
fenêtre où les résidus se voient, 40 poses :

| | dispersion | gain |
|---|---|---|
| brut | 364,7 | |
| LBL | 349,5 | −4,1 % |
| **dérivée** | **341,5** | **−6,4 %** |

## Deux améliorations, et la méthode marche

**1. Convoluer le modèle à la résolution de l'instrument.** TAPAS est bien
plus fin que HARPS (R = 115000, soit 2,6 km/s), et des raies trop étroites ne
peuvent pas annuler des raies élargies : il reste à chaque raie un résidu en
forme de dérivée, précisément ce que le critère cherche à annuler. La
puissance est prise **avant** la convolution, comme dans
`tellu_derivative.py`. Gain : **−4,6 %** de dispersion à elle seule.

**2. Deux groupes, pas quatre espèces** (consigne du 2026-10-15) : l'eau d'un
côté, tout le sec de l'autre. C'est d'ailleurs le découpage que TAPAS donne
dans LBL, dont `get_tapas_spl` rend `ABSO_WATER` et `ABSO_OTHERS`. Séparer
CO2, CH4 et O2 n'aurait pas de sens dans le visible, où seul O2 compte
vraiment et où les trois suivent la même masse d'air.

Avec les deux exposants libres, sur 40 poses :

| | |
|---|---|
| eau, médiane | 0,317 |
| sec, médiane | 1,094 |
| masse d'air, médiane | 1,20 |
| **sec / masse d'air** | **0,971** |
| eau sur la borne basse | **2 sur 40** |

Deux conclusions. D'abord **le blocage disparaît** : 2 poses sur 40 au lieu de
35 %. Ensuite **fixer le sec à la masse d'air est justifié** : ajusté
librement il vaut 0,971 fois la masse d'air, et le fixer ne coûte rien du tout
(+0,0 % de dispersion). C'est bien la physique, pas un pis-aller.

## Le résultat, sur les 100 poses

Dispersion du gradient dans 593-595 nm, la fenêtre où les résidus se voient :

| | dispersion | |
|---|---|---|
| brut | 344,8 | |
| **LBL** | 359,6 | **+4,3 %**, il DÉGRADE |
| **méthode par la dérivée** | **304,1** | **−11,8 %** |

## Ce qu'il a fallu pour que les bornes ne soient jamais atteintes

« Il y a toujours de l'eau », « les bornes ne doivent pas être atteintes »,
« basé sur des dizaines de milliers de spectres IR » (2026-10-15). C'est le
critère qui a guidé le réglage, et trois choses ont été nécessaires :

1. **La formulation du code qui marche** : la profondeur optique
   `tau = -ln(T)` sur une grille log uniforme, puis exponentiation, et la
   convolution à la résolution de l'instrument APRÈS la puissance. Une
   première version qui interpolait des transmissions déjà exponentiées sur
   la grille irrégulière du spectre tapait la borne sur 35 % des poses.
2. **L'ajustement en logarithme** de l'exposant : `alpha` reste strictement
   positif quoi que fasse le simplexe, et une variation relative compte
   autant en haut qu'en bas de la gamme. En linéaire, une pose tombait
   exactement à zéro.
3. **Un garde-fou sur les poses sans signal.** Il reste 9 poses sur 100 dont
   l'exposant s'effondre, et ce sont celles dont le SNR est faible : leur
   courbe de coût est **plate** (17,76 à 17,94 entre alpha 0,005 et 0,8, soit
   1 % de variation, là où une pose qui marche montre un minimum net). Le
   critère n'a pas de quoi trancher, donc les forcer inventerait une valeur :
   elles reçoivent la **médiane de la campagne**.

Après quoi l'exposant d'eau va de **0,091 à 1,097**, et **aucune pose n'est
sur une borne**.

Le réglage retenu :

| | |
|---|---|
| H2O | libre par fichier, ajusté en ln, sur les trois bandes d'eau |
| sec | **fixé à la masse d'air** : libre, il vaut 0,961 fois la masse d'air et ouvre un second mode d'échec (0,19 pour une masse d'air de 1,2) |
| vitesse | non libre, modèle dans le repère du sol |
| poses sans signal | médiane de la campagne, seuil à 0,05 |

`tools/tellu_deriv.py` porte l'implémentation : `fit_exponents` pour une
pose, `fit_campaign` pour une campagne avec le garde-fou, `model_grid` et
`transmission` pour le modèle, `objective` et `weighted_nanstd` pour le
critère de `telluric_fit` repris tel quel.
