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

## Le greffon dans LBL, et ce qu'il a fallu corriger

La méthode est maintenant **dans LBL** (`lbl/science/tellu_derivative.py`,
activée par `TELLUCLEAN_DERIVATIVE_WATER`) et seulement pour l'optique :
`Harps_ESO`, `HarpsN_ESO` et `Espresso`. SPIRou et NIRPS arrivent déjà
corrigés d'APERO et ne passent jamais par ce recipe.

**Le garde-fou devait descendre d'un étage.** Il vivait dans `fit_campaign`,
qui remplace par la médiane de campagne les poses dont l'exposant s'effondre.
Mais LBL appelle `fit_exponents` **une pose à la fois** et n'a aucune médiane
sous la main : sur les 34 premières poses de HD 10700, trois sont passées au
travers, dont deux à 0,0001, c'est-à-dire aucune correction d'eau.

Le test est donc maintenant **dans `fit_exponents`**, avant le simplexe : on
balaie neuf valeurs d'alpha de 0,02 à 2,5 et on regarde la **forme** de la
courbe de coût. Mesure sur deux poses de HD 10700 de la même nuit
(2003-09-22, 24 500 pixels dans les bandes d'eau) :

| pose | masse d'air | coût à alpha=0,005 | coût à alpha=2,0 | variation | minimum |
|---|---|---|---|---|---|
| 04:00 | 1,23 | 340,1 | 408,0 | **41 %** | intérieur, vers 0,5 |
| 06:26 | 1,03 | 17,76 | 18,23 | **2,6 %** | **monotone** |

La seconde n'a pas de minimum du tout : son coût monte dès la plus petite
valeur essayée, donc le simplexe la pousse à zéro quelle que soit la
paramétrisation. Les deux signatures sont séparées d'un ordre de grandeur,
donc une pose est refusée si sa courbe varie de moins de 5 % **ou** si son
minimum tombe sur un bord du balayage ; elle reçoit alors 0,35, la médiane
des poses qui convergent.

Le simplexe part désormais du meilleur point du balayage, et son résultat est
rejeté s'il sort de la plage balayée.

Vérification sur ces deux mêmes poses : celle qui marche garde **0,5936** (elle
donnait 0,6074 avec l'ancien point de départ), celle qui échoue reçoit
**0,3500** au lieu de 0,0001.

## Ce que les huit poses récalcitrantes étaient vraiment

Après le garde-fou, huit poses sur 100 recevaient encore le repli. J'ai
successivement supposé un manque de signal, un signe inversé dans le passage
au repère du sol, puis des raies stellaires qui noieraient le critère. **Les
trois hypothèses étaient fausses**, et chacune est écartée par une mesure :

- **le SNR est normal.** Six des sept poses au bord ont un SNR de 202 à 319,
  pour une médiane de campagne de 263, et une variation de coût saine de 31 %.
- **le signe est bon.** En balayant le sens du décalage, les poses qui
  marchent donnent le coût minimal en *divisant* la grille par le facteur
  Doppler : 505 contre 639 et 659 à BERV +23,7 km/s ; 394 contre 541 et 534 à
  +3,5 km/s. C'est bien le repère du sol.
- **les raies stellaires n'y sont pour rien.** Restreint aux 1045 pixels où
  TAPAS prédit l'eau la plus profonde, le critère trouve franchement 0,2 et
  0,6 sur les poses d'août et d'octobre, et reste collé au bord sur celles de
  décembre.

**Le balayage en vitesse a donné la réponse** : les poses qui marchent veulent
le modèle décalé de 0,0 km/s, celles qui échouent de **+16,0 km/s**, la même
valeur pour les trois testées. Et l'en-tête le confirme : `HIERARCH ESO INS
MODE` vaut **EGGS** sur ces poses, la fibre de 1,4" à haute efficacité, et
HARPS sur les autres. Le compte tombe juste : 8 poses EGGS sur 100, et ce sont
exactement les 8 qui échouent.

Ce n'est pas un élargissement des raies : baisser la résolution du modèle de
115 000 à 50 000 ne récupère rien, et leur coût résiduel reste de 1142 à 1593
contre 505 pour une pose HARPS.

`quality.drs_modes` les écarte donc en amont (`[HARPS]` pour cet instrument,
`null` par défaut). Elles sont encore nettoyées, puisque LBL ne connaît pas ce
réglage, mais `rms_table.py` compare sur les **poses communes** aux deux
séries, donc elles tombent de la comparaison de part et d'autre.

## Le bord du balayage cassait onze poses saines

Le garde-fou décrit plus haut mettait **19 poses sur 100** au repli, et le mode
EGGS n'en explique que 8. En étendant le balayage de 0,02 à **0,002**, les onze
poses HARPS restantes montrent toutes un minimum **franc et intérieur**, de
0,2 à 1,0, avec 25 à 45 % de variation de coût. Elles n'avaient aucun
problème : le bord du balayage tombait du mauvais côté, leur contraste de coût
(1,6 à 2,2) étant plus faible que celui des poses faciles (2,7). Un balayage en
vitesse le confirme : ces onze veulent le modèle à 0,0 km/s, comme les poses
qui marchent, donc leur repère n'est pas en cause.

Sur les 100 poses :

| | bord à 0,02 | bord à 0,002 |
|---|---|---|
| poses au repli | 19 | **2** |
| dont EGGS | 8 | 1 |
| dont HARPS | 11 | **1** |
| eau médiane | 0,332 | 0,301 |
| eau minimale | 0,0994 | 0,0039 |

La seule pose HARPS encore au repli est celle du 2003-09-22T06:26, SNR 6,6,
dont la courbe de coût ne varie que de 3 % : le vrai échec, celui pour lequel
le garde-fou a été écrit. Sur les poses HARPS seules, l'eau va de 0,0932 à sa
médiane de 0,304, sans toucher aucune borne.

## Correction: il n'y a aucun décalage en vitesse, nulle part

J'ai écrit plus haut que les poses EGGS voulaient le modèle TAPAS décalé de
+16 km/s. **C'était faux**, et `tools/tellu_shift.py` le montre.

L'erreur était d'estimateur : je lisais le minimum du coût par un `argmin`
**global** sur un balayage borné. Or la courbe de coût remonte moins haut d'un
côté que de l'autre, si bien que le dernier point du balayage peut passer sous
le creux. L'`argmin` lit alors le **bord**, pas le creux.

En prenant le minimum **local le plus proche de zéro**, qui est la quantité
qui a un sens :

| | argmin global | minimum local | coût à 0 / creux |
|---|---|---|---|
| 3 poses HARPS à bon SNR (183-337) | +16 à +18 km/s | **0,0** | 1,01 à 1,19 |
| 7 des 8 poses EGGS | +16 à +18 km/s | **0 à +3 km/s** | 1,00 à 1,03 |

La page 3 de la figure le confirme à l'œil : sur une pose à SNR 227, les raies
du modèle à 0,0 km/s tombent exactement sur celles du spectre, et celles à
+18 km/s sont visiblement à côté ; le spectre divisé par le modèle à zéro est
propre, divisé par celui à +18 il laisse des pics à 15.

**Le modèle est donc bien à sa place à vitesse zéro, et la vitesse ne doit pas
être un paramètre libre** (elle ne l'est pas : `fit_exponents` n'en a aucun).

### Ce qui reste vrai des poses EGGS, et ce qui reste à trouver

Leur ajustement d'eau échoue bel et bien : l'exposant y tombe entre **0,004 et
0,047**, cent fois moins que les 0,30 des poses HARPS, alors qu'il y a toujours
de l'eau au-dessus d'un télescope. Mais la cause n'est ni le repère, ni la
résolution (de R = 115 000 à 60 000, l'eau ne monte que de 0,0039 à 0,0089), ni
une pré-correction (les télluriques y sont aussi profondes que dans les poses
HARPS : creux de 0,135 contre 0,122 dans la bande B).

`quality.drs_modes: [HARPS]` les écarte donc parce que leur correction est
manifestement fausse, **pas** parce qu'on a compris pourquoi. C'est une
question ouverte.
