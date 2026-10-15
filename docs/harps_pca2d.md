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

## LE RIVER PLOT DIT QUE LA CORRECTION NE MARCHE PAS, ET IL A RAISON

Observation de l'utilisateur le 2026-10-15, en regardant les river plots : les
raies stellaires sont bien alignées et le panneau « minus the star block » les
enlève correctement, mais « la soustraction minus the observer block ne marche
pas », et il reste beaucoup de tellurique alors qu'on devrait être après la
correction tellurique de LBL.

**Les deux remarques sont justes, et la mesure va plus loin que l'oeil.**

D'abord, il n'y a eu AUCUNE correction tellurique de LBL : tous les runs ont
tourné avec `RUN_LBL_TELLUCLEAN=False`. Les telluriques sont donc intégralement
dans les spectres. Et pca2d ne les divise pas non plus : `weights.telluric_mode:
recon` ne fait que **pondérer**, il abaisse le poids là où la transmission est
basse.

Ensuite, sur la bande B de l'oxygène (686-691 nm), 25 poses appariées, en
écart-type de `ln(flux / continu)` colonne par colonne :

| | |
|---|---|
| profondeur tellurique moyenne | 0,469 |
| **variation d'une pose à l'autre, AVANT** | **0,0208** |
| **variation d'une pose à l'autre, APRÈS** | **0,0955** |
| correction appliquée | 0,0800 |

**La correction multiplie par 4,6 la variation tellurique au lieu de la
réduire.** Elle vaut 1,8 fois la variation qu'il fallait retirer : elle n'est
pas trop faible, elle est trop forte et mal placée.

Ce qui a été écarté en chemin, et qu'il faut noter pour ne pas y revenir :

- **le repère** : le cube est bien dans celui de l'observateur (85 BERV
  distincts de −27 à +25 km/s) et la tellurique y est fixe, corrélation 0,991
  entre chaque pose et le profil médian de la bande ;
- **la forme des composantes** : elles décrivent bien la bande (rms 0,0126
  dans la bande B contre 0,00025 sur le témoin Mg b, un facteur 50) ;
- **le rétrécissement** : il annule pourtant les composantes sur 63 à 72 % des
  colonnes, avec un facteur moyen de 0,17 à 0,24, mais le couper ne change
  rien aux fichiers écrits (0,08025 contre 0,08011) ;
- **la moyenne par parité** : elle est déjà soustraite des données du cube, et
  l'inclure ou non donne le même chiffre.

Un indice reste inexpliqué et mérite d'être creusé en premier : **plus il y a
de composantes, moins il se retire de variance tellurique pondérée** (6,6 % à
une composante, 4,4 % à deux, 3,4 % à trois, 2,9 % à sept). Un ajustement par
moindres carrés ne peut pas faire cela, donc quelque chose entre la solution et
son application ne tient pas.

## LA CHAÎNE COMPLÈTE : nettoyage tellurique, puis pca2d, puis LBL

Montée le 2026-10-15 à la demande de l'utilisateur, dans le bon ordre cette
fois. `lbl_telluclean` divise la transmission TAPAS ajustée pose par pose
(exposant d'eau et exposant sec, mesurés sur chaque spectre) et écrit dans
`<science>/<OBJECT>_tc/`. pca2d ajuste ensuite les **résidus** de ce nettoyage,
et LBL mesure les vitesses sur le tout.

| série | poses | rms | erreur | chi2 réduit |
|---|---|---|---|---|
| livrés, sans nettoyage | 67 | 1,57 | 0,39 | 27,2 |
| livrés + pca2d 0-2 | 65 | 1,79 | 0,92 | 5,3 |
| **nettoyés, sans pca2d** | 67 | **1,58** | 0,39 | 27,6 |
| **nettoyés + pca2d 0-2** | 66 | 2,09 | 0,93 | 6,6 |

Sur les 64 poses communes : 1,56 / 1,80 / 1,57 / 1,79 m/s.

**Le nettoyage tellurique ne change rien aux vitesses** (1,56 → 1,57 m/s), et
pca2d coûte la même chose avant comme après (1,80 contre 1,79).

### Mais pca2d fait bien son travail sur les résidus

Dans la bande B de l'oxygène, variation d'une pose à l'autre de
`ln(flux / continu)`, **sur les mêmes 3000 colonnes et les mêmes 94 poses** :

| | |
|---|---|
| nettoyés par LBL | 0,368 |
| **nettoyés + pca2d** | **0,094** |

**pca2d retire 74 % de la variation tellurique résiduelle.** Il fait donc
exactement ce qu'on lui demande dans la bande.

Et il faut corriger ici une mesure que j'avais annoncée plus haut : le
« ×4,6 » comparait des jeux de pixels différents, parce que la correction met
des échantillons à NaN. Sur des pixels identiques, le signe s'inverse.

### Ce que cela dit

Deux faits qui ne se contredisent pas : **pca2d nettoie très bien la bande
tellurique, et cela ne se voit pas dans les vitesses.** Sur HD 10700 le résidu
tellurique n'est tout simplement pas ce qui limite la mesure : une G8V très
calme observée à 1,3 m/s par le DRS, dont les raies profondes sont ailleurs
que dans les bandes de l'oxygène. Ce que la correction gagne dans la bande,
elle le perd en bruit ajouté sur le reste du domaine (l'erreur par pose passe
de 0,39 à 0,92 m/s).

C'est la même conclusion que `berv-span-decides-the-correction` sous un autre
angle : la correction se juge sur une cible où le résidu tellurique domine, et
HD 10700 n'en est pas une.

## POURQUOI L'INCERTITUDE TRIPLE : le masque, pas les pixels

« Ça n'a pas de sens d'avoir des incertitudes qui triplent en cours de route »
(2026-10-15). C'est juste, et la cause n'est pas celle qu'on croirait.

L'erreur par pose passe de 0,39 à 0,92 m/s, et elle se décompose en deux
facteurs qui se multiplient :

| | nettoyés | nettoyés + pca2d |
|---|---|---|
| raies utilisables | 10926 | **7239** (−33 %) |
| erreur par raie, médiane | 141 m/s | **244 m/s** (+73 %) |
| erreur combinée | 0,34 | 0,80 |

**Ce n'est pas le bruit des pixels.** Mesuré sur les mêmes pixels d'une même
pose : le bruit pixel à pixel ne bouge pas dans le bleu (−0 %), ni dans le
vert (−0 %), et **baisse de 4 %** dans la bande B. La correction n'abîme donc
pas les spectres.

**C'est le masque de raies.** Après pca2d, LBL en trouve PLUS mais de moins
bonne qualité :

| | nettoyés | nettoyés + pca2d |
|---|---|---|
| raies du masque | 17162 | **18799** |
| `line_snr` médian | 8,76 | **6,38** |
| raies de `line_snr` > 20 | 4766 | **3964** |

Et les raies en trop sont concentrées **dans le bleu**, là où il n'y a presque
pas de tellurique : +628 entre 380 et 430 nm, +344 entre 430 et 480, contre
+61 seulement entre 530 et 580.

**Le mécanisme, bout à bout.** Le `rms` du gabarit n'est pas une erreur
photonique mais la dispersion entre bacs de BERV (`lbl_template.py:556`,
`rms = (p84 − p16) / 2` sur 19 bacs). pca2d rend les poses plus semblables
entre elles, donc cette dispersion **baisse** : dans le bleu, 0,00296 contre
0,00315. Or `find_mask_lines` garde un extrêmum quand la dérivée y est
significative **devant ce rms**. Un seuil plus bas laisse donc passer des
extrêmums de bruit : 3954 extrêmums dans le bleu du gabarit corrigé contre
3111 dans celui des spectres nettoyés.

Autrement dit : **pca2d est pénalisé par le fait même de réduire la
dispersion**, parce que le seuil de sélection des raies est calculé sur cette
dispersion. Plus la correction marche, plus le masque se remplit de bruit.

**Ce qu'il faudrait.** Que le seuil de `find_mask_lines` ne soit pas lié à une
quantité que la correction modifie, ou que les deux séries soient mesurées sur
**le même masque** (celui des spectres non corrigés, ou un masque bâti sur le
modèle PHOENIX, voir `tools/mask_from_model.py`). C'est la première chose à
faire avant de rejuger une correction par ses vitesses, et cela vaut pour tous
les instruments, pas seulement HARPS.

## POURQUOI LE BLOC OBSERVATEUR NE MORD PAS : il ne pèse presque rien

« Après avoir enlevé le bloc stellaire, il reste une composante observateur
évidente » (2026-10-15, sur 686-690 nm). Vérifié, et l'observation est juste.

**Ce n'est pas un problème de décalage.** La structure résiduelle est bien
fixe dans le repère du cube : en corrélant chaque pose avec le profil médian
et en cherchant le décalage qui maximise l'accord, la pente contre le BERV
vaut **+0,0097 km/s par km/s** et le décalage médian est **nul**. Le cube est
dans le bon repère.

**C'est un problème d'ajustement.** Sur 686-690 nm, parité 0 (la seule qui
couvre cette bande), 94 lignes et 3087 colonnes :

| | |
|---|---|
| variance des données | 0,01231 |
| variance après le modèle pca2d | **0,01664** (+35 %) |
| corrélation données / modèle | **+0,367** |
| facteur d'échelle optimal du modèle | 0,762 |
| ce qu'une SVD locale à 2 composantes retirerait | **58 %** |

Le modèle **augmente** la variance au lieu de la réduire, et même ramené à son
échelle optimale il n'en retirerait que 13 %. Il a la bonne amplitude mais il
est mal placé : la structure résiduelle corrèle à 0,826 d'une pose à l'autre,
donc elle est bien là, et deux composantes correctement ajustées la
captureraient.

**La cause est dans les poids.** La bande ne pèse presque rien dans
l'ajustement global :

| région | poids médian | part du poids total |
|---|---|---|
| O2 bande B (686-690) | 1290 | **0,93 %** |
| O2 gamma (626-630) | 12900 | 1,84 % |
| vert (480-530), sans tellurique | 4160 | 8,81 % |

`weights.telluric_mode: recon` abaisse le poids là où la transmission est
basse, c'est-à-dire **exactement là où il faut corriger**. L'ajustement global
optimise donc le vert et ignore la bande tellurique, et les composantes qu'il
en tire n'y décrivent rien de bon.

**Une remarque au passage sur le rang.** Deux composantes n'expliquent que
58 % de la variance dans la bande B, mais aussi 50 % dans la gamma et **49 %
dans le vert sans tellurique**. Le rang est donc le même partout : ce qui
dépasse les deux premières composantes est du bruit réparti uniformément, et
non de la structure tellurique complexe. Cela confirme que deux composantes
sont le bon nombre, et que le problème est ailleurs.

**Ce qu'il faut essayer**, dans cet ordre :

1. **Ne pas écraser le poids des bandes telluriques** (`telluric_ramp_zero`
   plus bas, ou `telluric_mode` autre que `recon`) : tant que la bande pèse
   1 %, l'ajustement ne la verra pas.
2. **Ajuster par région** plutôt qu'une base unique sur 380-692 nm, pour que
   les composantes de la bande ne soient pas imposées par le vert.

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
13 pages. Les six fenêtres sont celles du tableau ci-dessus. La figure
`oh_residual` manque, et c'est normal : elle lit l'extension `OHLine`, que seul
SPIRou écrit.
