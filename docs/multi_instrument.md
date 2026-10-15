# Tout ce qui passe dans le LBL : HARPS, HARPS-N, ESPRESSO, SOPHIE

Journée du 2026-10-14. L'objectif était d'accepter d'autres instruments que
SPIRou et NIRPS, « à la limite en utilisant le module de lecture du LBL », et
d'aller jusqu'à une séquence LBL complète sur des données HARPS, avec et sans
correction. Les deux sont faits. Ce document dit comment, et surtout ce que la
mesure a répondu, qui n'est pas ce qu'on espérait.

## Résumé en trois lignes

1. **Le pont fonctionne.** pca2d lit un S2D de l'ESO par la classe de lecture
   de LBL, et LBL rend sur ces spectres **11,17 m/s de rms avec une erreur
   médiane de 11,04 et un chi2 réduit de 1,07** : un résultat de manuel pour
   tau Ceti sur HARPS.
2. **La correction, en l'état, dégrade :** 11,17 → 51,6 m/s. Ce n'est ni un
   réglage, ni un bogue.
3. **La raison est structurelle et était prévue.** Sur HARPS la grille est
   déjà barycentrique, donc le BERV vu par pca2d vaut zéro pour toutes les
   poses, et l'ajustement à deux repères n'a plus de levier. C'est le cas
   limite de ce que [`berv-span-decides-the-correction`] avait mesuré sur
   TOI-4552. Rendre le spectre au repère de l'observateur
   (`input.to_observer`) ramène le rms à **14,2 m/s** et confirme le
   diagnostic, mais la chaîne n'est pas finie : voir plus bas.

## Le pont : readers.py

Un seul module, `pca2d/readers.py`, qui rend pour n'importe quel instrument le
même dictionnaire que `tfits.read_tfits` : `flux`, `wave`, `blaze`, `recon`,
`sky`, `snr`, `meta`. Il n'implémente aucune lecture : il appelle la classe de
LBL, par `select.parse_args(..., parse=False)` puis `select.load_instrument`.
Vingt-deux combinaisons instrument/source s'instancient.

`read_any(path, config)` est le seul endroit où pca2d décide comment lire.
`cube.py`, `reconstruct.py` et `config.py` passent par là, donc ajouter un
instrument ne touche qu'à `reads_natively`.

| | APERO (t.fits) | ESO (S2D), SOPHIE (e2ds) |
|---|---|---|
| lecteur | `tfits.read_tfits` | la classe de LBL |
| transmission | extension `Recon` | modèle TAPAS de LBL, à la masse d'air |
| ciel par pose | extension `OHLine` | **rien** : `quality.max_sky_ratio` sans objet |
| SNR par ordre | `EXTSNxxx` de l'en-tête | mesuré sur le flux |
| blaze | extension | **vaut 1 partout**, et ce n'est pas une panne |

Le blaze plat mérite un mot : la classe HARPS de LBL ne déblaze pas, le S2D
est déjà plat de ce point de vue. Donc `input.blaze_min_frac: 0.0` sur ces
instruments, parce qu'une coupure relative au maximum du blaze ne couperait
pas un pixel.

## Les sept hypothèses d'APERO qu'il a fallu défaire

Chacune a coûté un arrêt de la chaîne, et chacune est maintenant un test.

**Le temps.** La classe HARPS de LBL fait pointer `KW_MJDATE` **et**
`KW_EXPTIME` sur le même mot-clé, `HIERARCH ESO QC BJD`, parce que LBL n'a
besoin que d'un axe temporel et se moque de la durée de pose. pca2d s'en sert
autrement : une pose de 60 s lue comme 2452859,8 fausserait toute pondération.
On lit donc `MJD-OBS` et `EXPTIME`, les vrais mots-clés, et la table de LBL ne
sert plus que de recours.

**Le contrat de `meta`.** `cube.py` lit `meta["bjd"]` sans savoir d'où vient le
spectre : le pont s'est arrêté dessus **après avoir lu 50 fichiers**. Les dix-
sept clés de `read_tfits` sont maintenant toutes là, NaN quand l'en-tête ne les
porte pas, et les sentinelles (-1 pour un seeing que l'ESO n'a pas, -9999,9
pour un mât météo en panne) deviennent NaN comme dans `read_tfits`.

**La transmission.** `order_quality` et `_order_sigma` indexent `recon` comme
un tableau : un `None` les ferait planter. Sans modèle, la transmission vaut 1
partout, c'est-à-dire « aucune absorption connue ici », et `meta` dit laquelle
des deux. Avec TAPAS, sur un S2D de tau Ceti à la masse d'air 1,135 :

| bande | transmission minimale | médiane |
|---|---|---|
| O₂ B (686-692 nm) | 0,000 | 0,974 |
| O₂ γ (627-632 nm) | 0,266 | 0,963 |
| H₂O (647-660 nm) | 0,491 | 0,978 |
| continu (500-520 nm) | 0,873 | 0,986 |

**`lbl.science` masqué.** Ce dépôt contient un dossier `lbl/` qui est une
arborescence de **données** de LBL (calib, masks, models, science...). Dès que
la racine du dépôt est sur `sys.path`, `import lbl` trouve ce dossier avant le
paquet installé et en fait un paquet-espace-de-noms : `lbl.__file__` vaut None
et `from lbl.science import tellu_clean` échoue avec « unknown location ». Et
c'est pire qu'un échec franc, parce que `lbl.instruments`, lui, n'est **pas**
masqué : l'instrument se chargeait, la transmission non, et rien ne le disait.
`_lbl_root()` retrouve le vrai paquet, y compris en installation éditable, où
il n'est sur aucune entrée de `sys.path` et où seul `importlib.metadata` sait
où il est.

**L'œuf et la poule du motif.** Le motif des fichiers (`*t.fits`,
`*_S2D_A.fits`, `*e2ds_A.fits`) est une propriété de l'instrument, qui se lit
dans un fichier, qu'il faut donc trouver avant de le connaître. `INSTRUME` est
dans tout FITS : on élargit à `*.fits` plutôt que de refuser.

**Deux noms pour une étoile.** `OBJECT = HD10700` en 2003, `tau-Cet` en 2004,
`HR509` en 2005 (SIMBAD les résout tous les trois vers tau Ceti). Cinq ans
d'archive, trois formulaires d'observation, et 36 poses sur 60 écartées en
silence tant qu'un seul nom était déclaré. `input.object_header` accepte
maintenant une liste.

**L'écriture.** Le fichier corrigé est une copie de l'original dont on remplace
le flux ; encore faut-il savoir quelle extension le porte. Retrouvée par la
forme des données, comme astropy l'a choisie, jamais devinée : `None`, et on
refuse d'écrire plutôt que d'écrire au mauvais endroit. Le suffixe que LBL
attend est préservé (`..._S2D_0-7_A.fits`, et non `..._S2D_At_0-7.fits`).

Vérifié sur un fichier corrigé : les huit extensions sont là, `WAVEDATA`,
`ERRDATA` et `QUALDATA` sont bit à bit identiques, seul `SCIDATA` a bougé.

## Le BERV : la mesure qui change tout

C'est l'erreur qu'il fallait trouver, et elle ne s'est vue que dans les
vitesses de LBL.

La classe `Harps_ESO` de LBL définit `get_berv` qui rend **zéro**, avec ce
commentaire : les fichiers de l'ESO sont déjà corrigés du barycentre. Le pont,
lui, lisait le BERV de l'en-tête (`HIERARCH ESO QC BERV`, jusqu'à 27 km/s) et
le donnait à pca2d, qui décalait donc les spectres une seconde fois.

Les sept poses au BERV le plus négatif sortaient alors à 1,3 et 4 km/s du
groupe, et le rms de LBL valait **1311 m/s**.

Corrigé : le BERV vient de `inst.get_berv()`, c'est-à-dire de LBL, qui a déjà
tranché classe par classe. SPIRou lit son en-tête, HARPS rend zéro. On pose la
question plutôt que d'y répondre, et le jour où un format change de convention
c'est LBL qui le saura.

**La vérification, par LBL lui-même :** la pente de `vrad` contre le BERV
d'en-tête vaut **+0,00002 km/s par km/s** (corrélation 0,018) sur 52 km/s
d'étendue, avec 11,17 m/s de résidu. La grille est bien déjà barycentrique.

## Ce que LBL mesure, avant et après

tau Ceti (HD 10700), 100 poses récupérées de DACE, 380-692 nm, 71 ordres,
7,3 ans, 71 nuits, SNR médian 88, pose médiane 30 s. Rms robuste, poses à plus
de 5 sigmas écartées de chaque série.

| série | poses | rms (m/s) | erreur (m/s) | chi2 réduit |
|---|---|---|---|---|
| **tels que livrés** | 67 | **11,17** | 11,04 | **1,07** |
| 0-7 nominal | 46 | 51,58 | 14,99 | 43,2 |
| 0-7 sans masque (`mask: none`) | 67 | 56,49 | 14,63 | 15,0 |
| 0-2 (deux composantes) | 67 | 314,41 | 15,43 | 435 |

Sur les 43 poses communes aux quatre : 11,54 → 52,93 (nominal), 50,14 (sans
masque), 265,04 (deux composantes).

**La ligne qui compte est la première.** Un chi2 réduit de 1,07 avec 11 m/s de
rms et 11 m/s d'erreur dit que le pont livre à LBL des spectres dont le modèle
de bruit est exact. Tout le travail de lecture est validé par là.

## Pourquoi la correction dégrade : ce qui a été éliminé, et ce qui reste

L'injection de la correction, pose par pose, vaut **55,8 m/s de rms** contre
une erreur photonique de 11. Elle ne corrèle ni avec le BERV (0,21) ni avec la
masse d'air (0,18) : ce n'est donc ni un résidu tellurique ni une signature
barycentrique, c'est du bruit.

Trois hypothèses testées et **écartées par la mesure** :

- **Le masque.** La correction met à NaN 35 % des pixels de raie profonde
  contre 4 % du continu, ce qui semblait désigner le coupable. Mais
  `mask: none`, qui ne blanchit rien, donne 50,1 m/s au lieu de 52,9 : le
  masque n'y est pour presque rien.
- **Le nombre de composantes.** Moins de composantes est bien **pire** (265
  m/s à deux contre 53 à sept) : le dégât suit la qualité de l'ajustement, pas
  la quantité retirée.
- **Le découpage en parités.** Les ordres de HARPS se recouvrent sur 28,5 % de
  leur largeur (1,56 nm sur 5,48), contre presque rien sur SPIRou, donc les
  deux parités décrivent souvent les mêmes photons. Mesuré sur 66 zones de
  recouvrement : les deux corrections ne diffèrent que de **0,042 % du flux**.
  Elles sont cohérentes.

**Ce qui reste, et que le chi2 de l'ajustement dit tout seul : le chi2 réduit
par spectre vaut 0,594.** L'ajustement explique les données *mieux* que le
modèle de bruit ne l'autorise, ce qui est la signature du surajustement. Avec
170 lignes et 1190 coefficients, et surtout **sans aucun contraste entre les
deux repères**, les sept composantes de l'observateur reproduisent la
réalisation de bruit de chaque pose, et la division la réécrit dans le
spectre.

Sur SPIRou le BERV s'étale sur des dizaines de km/s : une composante de
l'observateur est fixe en longueur d'onde pendant que l'étoile se déplace de
plusieurs pixels, et les deux blocs sont séparables. Ici les deux blocs
décrivent le même repère.

## La piste 1 essayée : remettre le spectre dans le repère de l'observateur

`input.to_observer: true` divise la grille par le Doppler **relativiste** du
BERV d'en-tête et rend ce BERV comme étant à appliquer. Le spectre redevient ce
que le spectrographe a vu, et l'ajustement retrouve son levier. Vérifié au
pixel : 23,7428 km/s appliqués, 23,7428 retrouvés.

L'ajustement change de nature : la première composante de l'observateur porte
**30,2 %** de la variance au lieu de 10,2 %, donc elle a trouvé une structure
réellement ancrée dans le repère du sol.

| série | poses | rms (m/s) | erreur (m/s) | chi2 réduit |
|---|---|---|---|---|
| tels que livrés | 67 | 11,17 | 11,04 | 1,07 |
| 0-7, BERV nul | 46 | 51,58 | 14,99 | 43,2 |
| **0-7, repère observateur** | 67 | **14,24** | 1,59 | 102,9 |
| 0-7, repère observateur, `mask: exposure` | 67 | 693,84 | 15,88 | 2259 |
| **0-7, repère observateur, sans shrink ni refit** | 67 | **12,84** | 1,50 | 88,2 |

Sur les poses communes : 11,54 → 13,73 m/s, contre 52,93 sans cette option.
**Le gros de la dégradation est donc bien venu de l'absence de contraste entre
les deux repères**, et l'hypothèse est confirmée.

**Mais ce n'est pas encore utilisable**, et il faut le dire : l'erreur tombe à
1,59 m/s pour un chi2 de 103, ce qui n'est pas croyable, et **44,9 % des pixels
finissent à NaN** contre 21,9 % avant.

### Pourquoi tant de NaN : l'intersection de 85 intervalles décalés

Ce n'est pas que le masque désigne les mauvaises raies, comme on pouvait le
croire : c'est qu'il en reste beaucoup moins. Dans le repère de l'observateur,
chaque pose couvre un intervalle de longueurs d'onde **décalé du sien propre**,
et 52 km/s séparent les extrêmes. `mask: common` ne garde qu'un échantillon que
**toutes** les poses ont pondéré, donc l'intersection de 85 intervalles
décalés :

| cube | échantillons vivants du masque commun |
|---|---|
| barycentrique (BERV nul) | 61,5 % |
| repère de l'observateur | **39,7 %** |

`mask: exposure`, qui ne demande rien aux autres poses, rend 63,3 % de vivants
**dans les deux repères**. Essayé, et c'est un désastre : **693,8 m/s de rms**.
Chaque pose garde alors un autre jeu de raies, et LBL moyenne des jeux
différents, exactement ce que [`corrected-files-equal-panel-3`] décrit et que
le nominal évite depuis le 2026-09-12. L'intersection est donc le prix à payer
pour un jeu de raies cohérent, pas un défaut à corriger.

### Les barres d'erreur de 1,6 m/s : ce que ce n'est PAS

Il faut raconter cette piste en entier, parce que la première réponse était
fausse et que la mesure qui l'a corrigée est simple.

Le bruit pixel à pixel du spectre corrigé vaut 0,85 de l'original, ce qui
semblait dire que l'ajustement absorbait la réalisation de bruit de chaque
pose. **C'est un artefact de comparaison.** La correction met 22 à 45 % des
pixels à NaN, et ce ne sont pas des pixels au hasard : ce sont les plus
bruyants. Mesurer le bruit avant sur tous les pixels et après sur les
survivants compare deux populations différentes.

**Sur les mêmes pixels**, le rapport vaut :

| série | pixels différents | **mêmes pixels** |
|---|---|---|
| 0-7, BERV nul | 0,851 | **0,999** |
| 0-7, repère observateur, sans shrink | 0,830 | 0,937 |

À 0,999, la correction ne lisse rien du tout. Le comptage des degrés de
liberté le disait d'ailleurs : sept composantes et une moyenne sur 85 lignes
par parité ne peuvent retirer que 5 % du bruit (sqrt(1 − 8/85) = 0,952), et
il aurait fallu 20 à 26 paramètres par colonne pour expliquer 0,85.

Ce qui reste vrai et inexpliqué : **LBL rend des incertitudes par raie 2,8 fois
plus petites** sur les spectres corrigés dans le repère de l'observateur (1249
contre 3517 m/s de médiane), alors qu'il y utilise *moins* de raies (8095
contre 11169). Les fichiers sont pourtant cohérents : `ERRDATA` est inchangée
et le rapport entre le bruit mesuré et l'erreur annoncée vaut 1,40 avant comme
après. La différence naît donc **dans LBL**, vraisemblablement dans le gabarit,
que chaque série reconstruit pour elle-même et qui est plus creux ici (72 % de
points finis contre 96,7 % pour la série livrée). C'est la piste à suivre, et
elle n'est pas dans pca2d.

Ni le rétrécissement ni le réajustement par pose n'y sont pour quelque chose :
les couper donne 12,84 m/s au lieu de 14,24, et laisse l'erreur à 1,50.

## Ce qu'il faut en retenir pour la suite

**Ne pas corriger un spectrographe dont la grille est déjà barycentrique, en
l'état.** HARPS, HARPS-N, ESPRESSO et SOPHIE sont tous dans ce cas (leur
`get_berv` rend zéro). La lecture, elle, est acquise et utile telle quelle.

Trois pistes, dans l'ordre où il faut les prendre :

1. **Comprendre les incertitudes par raie de LBL** sur un spectre à 28 % de
   NaN : c'est le seul point qui reste faux, et il est du côté de LBL, pas de
   pca2d. Le gabarit est le premier endroit à regarder.
2. **Réduire l'intersection du masque commun** sans passer à `mask: exposure`,
   qui est pire : n'exiger l'accord que d'une fraction des poses plutôt que de
   toutes rendrait une bonne part des 44,9 % de NaN.
3. **Ne garder que les composantes qui corrèlent avec une grandeur physique**
   (masse d'air, humidité, colonne d'eau), ce qui est déjà à moitié fait par
   `correlations.csv`.

## Ce qui est en place

- `pca2d/readers.py` : le pont, et `KNOWN`, la table des 22 combinaisons
- `config.yaml` : blocs `HARPS`, `HARPN`, `ESPRESSO`, `SOPHIE`, avec les
  résolutions de LBL (115 000, 115 000, 140 000, 75 000), vérifiées contre
  `APPROX_RESOLUTION`
- `output.model_directory` : où LBL garde `tapas_lbl.fits`
- `cube.ORDER_FORMATS` / `FLAT_FORMATS` : ce que `input.format` décide est la
  dimension, pas le nom du fichier
- `tests/test_readers.py` : 11 tests, qui tournent sans données et sans réseau
