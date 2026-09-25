# Masquer le résidu : ce qui a été essayé le 2026-09-25, et ce qui reste à faire

Une journée d'essais sur TOI-2120 à 0-7, tous sur **le même ajustement** (les
variantes ne changent que la correction, `reuse_fit: nominal`), tous mesurés
par LBL sur les mêmes 316 poses. Ce document existe pour que rien de tout cela
ne soit à refaire.

## Le résultat en une phrase

La seule méthode qui gagne franchement sur la dispersion (14%) est aussi celle
qui **déplace l'amplitude de la planète de 32%**, et toutes les méthodes fines,
mieux fondées, ne gagnent presque rien. Nous ne savons pas encore si le gain
est réel ou si c'est du signal mangé, et **c'est la question à trancher
avant d'adopter quoi que ce soit**.

## Les chiffres

316 poses communes, P = 5,7998164 j, ajustement circulaire pour K.

| essai | rms | robuste | par nuit | err | K ± σ |
| --- | --- | --- | --- | --- | --- |
| livré, sans correction | 47.70 | 35.14 | 46.86 | 7.29 | 3.90 ± 3.54 |
| **0-7 de départ** | **17.01** | **13.79** | **14.57** | 7.89 | **5.42 ± 1.23** |
| 4 sigma par échantillon | 16.67 | 13.14 | 14.25 | 7.94 | 5.37 ± 1.22 |
| régions observateur | 16.74 | 12.73 | 14.38 | 8.10 | 5.85 ± 1.22 |
| 4 sigma + régions (`0-7s40r`) | 16.70 | 12.91 | 14.34 | 8.14 | 5.82 ± 1.22 |
| 3 sigma | 15.64 | 12.97 | 13.00 | 8.22 | 5.29 ± 1.14 |
| **2,5 sigma** | **14.60** | **12.01** | **11.72** | 9.06 | **3.67 ± 1.08** |
| 2,5 sigma + régions | 14.64 | 12.31 | 11.72 | 9.30 | 4.34 ± 1.07 |
| 2 sigma | LBL refuse : moins de 50% du spectre valide dans le gabarit | | | | |
| excursions 6 sigma | 17.05 | 13.80 | — | 7.89 | — |

Gains en rms sur le 0-7 de départ, par bootstrap apparié (4000 rééchantillons
des poses, ce qui respecte la corrélation entre séries) :

    2,5 sigma            +2,41 m/s  [+1,52, +3,28]   P = 1,000
    3 sigma              +1,37 m/s  [+0,80, +1,94]   P = 1,000
    4 sigma              +0,34 m/s  [+0,09, +0,58]   P = 0,999
    régions observateur  +0,27 m/s  [-0,12, +0,66]   P = 0,908
    excursions 6 sigma   -0,04 m/s  [-0,13, +0,05]   P = 0,195

## Ce que chaque critère cherche, et ce qu'il vaut

### 1. La coupure par échantillon (`correct.nsig_cut`)

Tout échantillon dont le résidu dépasse N sigmas robustes glissants devient
NaN. Brutal, sans sélection, et **c'est ce qui gagne le plus**.

Le gain croît continûment quand le seuil baisse (4 → 3 → 2,5 sigma donne +0,34
→ +1,37 → +2,41 m/s) pendant que l'erreur rapportée par LBL monte (7,94 → 8,22
→ 9,06 m/s) et que le chi2 réduit descend vers 1 (4,02 → 3,28 → 2,31). C'est la
signature d'un **écrêtage d'une distribution à ailes lourdes**, pas d'une
détection de défauts.

**La borne dure est à 2 sigma** : LBL s'arrête, « Less than 50% of the data in
"flux" valid in the template », l'ordre 8 tombant à 50,5% de données valides.

**Le drapeau rouge : à 2,5 sigma, K passe de 5,42 à 3,67 m/s**, soit 1,6 σ plus
bas, et chaque vitesse individuelle bouge de 4,9 m/s en médiane (corrélation de
0,906 avec la série de départ). Les raies profondes, celles qui portent la
vitesse, sont aussi celles où le modèle est le moins bon : la coupure les mord.
À 3 sigma, K reste à 5,29, stable.

### 2. Les excursions corrélées (`correct.excursion_nsig`, `excursion_elements`)

Idée : un échantillon isolé au-delà de 3 sigma n'est pas un événement (0,27% du
bruit pur l'est, 1% de ces résidus l'étaient), alors qu'une suite d'échantillons
penchant du même côté sur la largeur d'une raie en est un. Chaque fenêtre, de 1
échantillon à deux éléments de résolution, est jugée sur la significativité
agrégée de sa somme, `Z(w) = |somme des z| / sqrt(V(w))`, et marquée entière
au-delà de 6 sigma.

**V(w) est mesuré, jamais supposé** (`outliers.window_variance`,
`noise_correlation`) :

    V(w) = w + 2 * somme sur k de (w - k) * rho_k

Sur TOI-2120 : rho_1 = 0,935, rho_2 = 0,766, rho_4 = 0,303. Donc V(4) = 13,8 et
non 4 : **quatre échantillons de la grille à 3 sigma valent 3,2 sigma, pas 6**,
parce qu'ils ne font qu'un seul pixel. Sur 17 échantillons (deux éléments de
résolution, quatre pixels) V = 83, soit 3,5 échantillons indépendants, et une
bosse plate de 3,2 sigma de profondeur atteint 6. C'est bien « quatre pixels
indépendants à 3 sigma », avec l'arithmétique faite.

Vérifications de calibration, toutes bonnes :
- longueur de corrélation intégrale 1 + 2Σrho = **4,42 échantillons = 2,21
  km/s**, contre un pixel SPIRou de 2,28 km/s. Le critère compte donc bien en
  pixels indépendants, ce qui était l'exigence.
- écart-type robuste de Z égal à 1,00 / 1,00 / 1,01 / 1,02 pour w = 1, 4, 9, 17.

**Et pourtant : gain nul** (−0,04 m/s). 236 469 échantillons marqués, soit
0,07%. Ce qui coûte des vitesses n'est donc **pas** la structure rare et forte.

### 3. Les régions ancrées dans le référentiel de l'observateur (`correct.excess_nsig`)

Ce que montrent les river plots du rapport LBL, et ce qui a été vu à l'œil à
1541,5 et 1249 nm. Une fuite ancrée dans le référentiel de l'observateur n'a pas
de signe privilégié, mais elle a toujours une **dispersion** que le bruit
n'explique pas : chi2(j), la moyenne des z² sur les poses à cette colonne,
au-dessus de 1. Moyenné sur deux éléments de résolution, ce chi2 est certain à
`sqrt(2 / (N * colonnes indépendantes))` = 0,04 avec 316 poses, donc l'excès
devient une significativité.

Une région part de **toutes** les poses quand elle dépasse `excess_nsig`
sigmas d'excès **et** un chi2 de `excess_chi2`. Pas de condition sur la
fraction de poses écrêtées : c'est cette condition-là qui stérilisait la
première version (`column_frac`), en exigeant qu'une colonne soit à la fois
mauvaise en moyenne et catastrophique dans 10% des poses.

Ce que ça trouve sur TOI-2120 (`scratchpad/observer_rms.py`, carte
`observer_rms.pdf`) :

| région | largeur | chi2 de colonne au pic |
| --- | --- | --- |
| 1259-1274 nm (bande O2 à 1,27 µm) | 17 à 39 échantillons | jusqu'à **222** |
| 1919,3-1919,5 nm (eau) | 62 échantillons | 63 |
| 1249 nm (vu à l'œil) | — | 1,77, soit 8 sigma d'excès |
| 1541,5 nm (vu à l'œil) | — | 1,46, soit 5 sigma d'excès |

5,7% des colonnes dépassent 20 sigmas d'excès, 1,2% dépassent 50, en 1368
suites de quatre échantillons ou plus. Au seuil retenu (20 sigmas et chi2 >
1,25), 4,75% des colonnes partent.

**Gain : +0,27 m/s sur la rms (non significatif), mais +1,06 m/s sur la rms
robuste** (13,79 → 12,73). Ça nettoie le cœur de la distribution, pas les poses
aberrantes. Et **K se déplace de +0,4 m/s systématiquement**, dans le sens
opposé à la coupure agressive ; c'est dans le bruit, mais c'est reproductible
(5,85 seul, 5,82 avec la coupure à 4 sigma).

### 4. Les combinaisons

**Elles ne s'additionnent pas.** 4 sigma seul donne 16,67, les régions seules
16,74, les deux ensemble 16,70. Idem à 2,5 sigma : 14,60 seul, 14,64 avec les
régions, et la rms robuste est même moins bonne (12,31 contre 12,01). Les deux
critères mordent les mêmes échantillons ; en ajouter un ne fait qu'enlever
davantage de données pour rien.

## Une structure qui échappe à tous ces critères

À **1267,845 nm**, le résidu médian vaut +0,57 sigma et il est positif dans
**72% des poses**, sur environ un élément de résolution, alors que **le chi2 de
la colonne reste à 1,0-1,2**. C'est un **biais**, pas un excès de dispersion :
ni la coupure ni les régions ne le voient, par construction.

Le ciel n'y contribue rien (rapport ciel/flux à 0,000) et la transmission
tellurique y vaut 0,956, donc 4% d'absorption. L'hypothèse la plus simple est
une résiduelle de la correction tellurique en bord de bande O2.

**Cela suggère un cinquième critère, non essayé : le biais médian par colonne**
plutôt que l'excès de dispersion. La carte existe déjà
(`scratchpad/observer_map.py`, `observer_map.pdf`) : sur 316 poses, une
structure à +0,57 sigma ressort à 10 sigma. Elle avait été écartée trop vite.

## Ce qu'il faut faire ensuite, par ordre d'importance

1. **Injecter une planète synthétique de K connu** et voir laquelle des
   variantes la restitue sans biais. C'est le seul instrument qui sépare
   « enlève du bruit » de « mange du signal », et sans lui le choix entre 3
   sigma et 2,5 sigma n'est pas décidable. Tout le reste est secondaire.
2. **Masquer pendant la boucle plutôt qu'après** (idée d'Étienne, 2026-09-25).
   Aujourd'hui l'ajustement tourne jusqu'au bout avec des échantillons qu'il ne
   sait pas décrire, et on nettoie à la fin : les composantes ont donc déjà été
   tirées par la bande O2 à 1267 nm, et cette distorsion reste dans la base
   même après masquage. Le résidu et son sigma local sont calculés à chaque
   balayage de toute façon.
   - **Un poids continu, pas un masque binaire** : `w / (1 + (z/z0)²)`, dans la
     logique de `correct.shrink`. Un seuil dans une boucle est un système qui
     bascule et qui peut osciller ; un poids continu est stable et ne crée pas
     de discontinuité là où le modèle est déjà fragile.
   - **Geler la pondération après le deuxième balayage**, sinon la boucle se
     mord la queue : un échantillon déprécié est moins bien modélisé, donc
     encore plus déprécié.
   - **Attention à la clé de cache** : un masque qui dépend de l'ajustement
     n'est plus décrit par la configuration, et deux passages nominalement
     identiques pourraient diverger.
   - **Le poids doit percoler jusqu'aux fichiers corrigés** (Étienne,
     2026-09-25). Il y a deux masques et il ne faut pas les confondre : le
     poids dans l'ajustement dit « cet échantillon ne doit pas tirer la base »,
     le NaN du fichier corrigé dit « LBL ne doit pas mesurer de vitesse
     là-dessus ». Baisser le poids sans propager laisserait un échantillon jugé
     indigne de contraindre le modèle être livré à LBL comme une mesure
     ordinaire, corrigé par un modèle qu'il n'a pas aidé à construire.
     - **Un seul seuil, à 50% du poids nominal** (Étienne, 2026-09-25) : les
       deux masques sont alors le même objet. Un échantillon écouté à moins de
       la moitié n'est pas une demi-mesure, c'est une mesure que le modèle n'a
       pas vraiment décrite, et elle devient NaN. Cela évite d'avoir un
       plancher de poids d'un côté et un seuil de sortie de l'autre : une
       valeur dit à la fois jusqu'où on écoute et à partir d'où on jette.
     - Le chemin existe déjà : `reconstruct.fit_weights_mask` met déjà à NaN ce
       que l'ajustement n'a pas pondéré (`PCA2WNAN`, 26 550 échantillons par
       fichier sur ce passage). Un poids sous le plancher rejoindrait ce même
       masque, donc le panneau 3 de la figure de séquence continuerait de
       montrer ce que les fichiers contiennent, ce qui est la règle du dépôt.
     - Le masque de sortie doit rester **commun à toutes les poses**
       (`correct.mask: common`, le nominal) : on masque là où le poids MOYEN
       sur les poses tombe sous les 50%, sinon chaque pose porte un jeu de
       raies différent et LBL mesure des séries incomparables.
     - Et le poids final de la boucle est sans doute **un meilleur critère de
       masquage** que les statistiques essayées ici : il agrège ce que
       l'ajustement a appris sur tous les balayages, au lieu de juger le
       résidu final une seule fois.
   - Gain attendu : **faible sur la rms** (ces régions ne valent que +0,27 m/s
     enlevées après coup), mais c'est la **stabilité de K** qu'on cherche, pas
     la dispersion.
3. **Le biais médian par colonne**, critère 5 ci-dessus.
4. **Vérifier sur un deuxième objet** avant toute adoption : Proxima et
   TOI-4552 se comportent à l'opposé de TOI-2120 sur d'autres réglages.

## Où sont les choses

| quoi | où |
| --- | --- |
| les variantes | `variants/0-7point.yaml`, `0-7s25`, `0-7s40`, `0-7s20`, `0-7chi2`, `0-7exc`, `0-7rms`, `0-7both`, `0-7s40r` |
| le code du résidu | `pca2d/outliers.py` : `excursions`, `window_variance`, `noise_correlation`, `excess_rms` |
| les réglages | `correct.nsig_cut`, `excursion_nsig`, `excursion_elements`, `excess_nsig`, `excess_chi2`, `excess_elements` dans `config.yaml` |
| les en-têtes | `PCA2RSIG`, `PCA2RNAN` (coupure), `PCA2XSIG`, `PCA2XWID`, `PCA2XRHO` (excursions), `PCA2ESIG`, `PCA2ECHI`, `PCA2EWID`, `PCA2ECOL` (régions) |
| comparaison des vitesses | `outputs/TOI2120/rv_comparison_0-7_clip.pdf` |
| avant/après par fenêtre | `outputs/TOI2120/residuals_0-7_vs_0-7s40r.pdf` (une parité par fenêtre, celle la plus proche du centre de l'ordre ; en jaune ce qui devient NaN) |
| rapports LBL complets | `outputs/lbl/lblreport/<objet>_<objet>/lbl_report_*.pdf`, 34 pages, recette `lbl_report` de la branche `developer-report` de LBL |
| cartes de diagnostic | `scratchpad/observer_rms.py` et `observer_map.py` (à refaire tourner : le scratchpad est effacé périodiquement) |

## Détail de méthode qui a coûté du temps

- **Une variante ne doit jamais écrire dans le dossier du passage dont elle
  réutilise l'ajustement.** Avec `--out-dir` explicite, le dossier `_<nom>` de
  la variante est ignoré et elle écrase les spectres corrigés de l'autre. Cela
  est arrivé (19 fichiers sur 316) et c'est maintenant refusé
  (`cli.refuse_to_overwrite_the_reused_fit`).
- **Un dossier nommé autrement que l'en-tête** a besoin de
  `input.object_header`, sauf si le suffixe est le nom de l'instrument, auquel
  cas c'est automatique depuis le 2026-09-18.
- **Les figures en river plot : une seule parité d'ordre par fenêtre**, celle
  dont on est le plus près du centre de l'ordre (sommet du blaze, distance
  relative à la demi-largeur). Empiler les deux parités fait des lignes
  blanches qui ressemblent à des poses perdues et n'en sont pas.
