# pca2d-preclean : un jeu de paramètres nominal, reproductible

*Ce document vivait dans `outputs/reports/`, qui est dans le `.gitignore` :
il enregistre le nominal et les décisions, donc il est suivi par git depuis le
12 septembre 2026. Son ancien chemin n'existe plus.*

Proposition du 2026-09-11, à confirmer ligne par ligne (☐ → ☑, ou corrigez la valeur).

## Principe

1. **On adopte un jeu de paramètres raisonnable,** écrit à un seul endroit, que tout passage reproduit exactement. Les autres options restent disponibles ; on ne retire rien tant qu'on n'a pas compris (section 4, pour plus tard).
2. **Toute variante est le nominal plus un écart nommé,** versionné dans le dépôt. Plus aucune configuration dans un dossier temporaire.
3. **Chaque passage enregistre de quoi le refaire :** la configuration résolue (déjà fait) et la version du code (absente aujourd'hui).
4. **Un choix scientifique n'est compris que mesuré aux deux extrêmes :**
   - TOI-2120 : SPIRou, étoile faible, SNR bas, 316 poses individuelles ;
   - Proxima : NIRPS, SNR bien meilleur, 782 poses empilées en 259 nuits.

---

## 1. Le jeu nominal : ce que vaut chaque choix

### Ajustement

| Choix | Valeur adoptée | Autres valeurs disponibles | Mesuré sur | OK |
|---|---|---|---|---|
| composantes étoile / observateur | 1 / 3 | tout nombre | TOI-2120 (2 étoile dégénère : 67 m/s) | ☐ |
| partie statique (`mean`) | `star` : un spectre stellaire par parité | offset, full, iterate | TOI-2120 : 20,36 contre 31,74 (offset). Proxima : offset injecte ~47 m/s ; star en cours | ☐ |
| représentation de l'étoile | spline cubique exacte | grille Lanczos | TOI-2120 : égalité en LBL, plus rapide. Proxima en cours | ☐ |
| Savitzky-Golay sur le gabarit stellaire (`star_smooth`) | **aucun pour l'instant** | fraction d'élément de résolution ; propagé au gabarit LBL | TOI-2120 : ¼ et ½ élément en cours. Proxima : à faire | ☐ |
| terme de vitesse | désactivé | activé | par choix (option gardée pour plus tard) | ☐ |
| règle d'arrêt | **à adopter** (section 3) | — | les deux cibles se comportent à l'opposé | ☐ |

Valeurs techniques adoptées telles quelles (toujours modifiables) :
- **Rejets :** 10 σ sur les coefficients, 2 tours ; poses sous 0,5 × le SNR médian ; écrêtage doux à 3 σ.
- **Gabarit :** groupement par BERV de 3000 m/s, 3 poses au minimum par groupe.
- **Calcul :** étoile mise à jour d'abord ; amplitudes partagées par les deux parités ; noyau Lanczos de 8 échantillons ; garde des trous activée ; float32.

☐

### Correction des spectres

| Choix | Valeur adoptée | Autres valeurs disponibles | Mesuré sur | OK |
|---|---|---|---|---|
| rétrécissement des composantes non significatives | activé | désactivé | TOI-2120 : robuste 16,16 contre 17,28. Proxima en cours (6ᵉ panneau) | ☐ |
| significativité lissée | désactivée | activée | TOI-2120 : 21,17 / 16,39, moins bon | ☐ |
| lissage de composantes observateur | aucun | liste | TOI-2120 : 21,00 / 16,61, moins bon ; contraire à la règle (défauts de pixels) | ☐ |
| composantes divisées | toutes les composantes observateur, aucune étoile | M-N quelconque | par principe | ☐ |
| NaN sur les résidus au-delà de N σ | désactivé | N | jamais mesuré récemment | ☐ |

### Cube

| Choix | Valeur adoptée | Autres valeurs disponibles | OK |
|---|---|---|---|
| pas de la grille | 0,5 km/s | « smart », à partir des pixels (section 3) | ☐ |
| passe-haut | Savitzky-Golay de 100 km/s, ordre 2, sur ln f | autre largeur | ☐ |
| empilement par nuit | auto ; NIRPS toujours | oui, non | ☐ |
| rejet du ciel | là où le ciel soustrait dépasse 4 × l'étoile | autre seuil, aucun | ☐ |

Valeurs techniques adoptées telles quelles :
- **Poids :** photon ; bruit mesuré sur 201 échantillons (thermique en K) ; transmission lue dans Recon, rampe de 0,5 à 1 ; raies sous −0,5 en ln f sans poids.
- **Rejets :** échantillons isolés (3) ; 1 pixel au bord des trous ; 90 % de NaN au maximum ; colonnes avec moins de 20 % de spectres ; blaze sous 5 % ; ordres sous 5 % de pixels finis.

☐

### LBL

- Les deux objets sont mesurés : les spectres livrés et les spectres corrigés.
- L'objet corrigé est mesuré contre le gabarit de l'ajustement, lissé comme l'étoile.
- Teff lue dans les en-têtes ; chaque objet a son propre gabarit ; liens symboliques.

☐

---

## 2. Reproductibilité : ce qui manque aujourd'hui, à corriger maintenant

1. **Deux sources de vérité qui divergent.** Les valeurs par défaut du code et `config.yaml` diffèrent sur 12 clés (tableau ci-dessous). Proposition : les défauts du code deviennent le nominal, et `config.yaml` ne garde que le site, l'instrument et l'objet. ☐
2. **Des variantes lancées depuis un dossier temporaire.** C'est le cas des configurations ss025, ss05 et star*, et du script des variantes de correction, qui disparaîtront avec la session. Proposition : un dossier `variants/` dans le dépôt, un fichier par variante qui ne contient que ses écarts au nominal, et un objet LBL qui porte son nom (`TOI2120_PCA2D_1-3_<variante>`). Une variante qui ne touche que la correction réutilise l'ajustement. ☑ **Fait (aacc6b1).** `--variant NOM`, `reuse_fit`, et les 12 passages de la série TOI-2120 dans `variants/` ; chacun redonne sa configuration et le cube a92a1b985cdd.
3. **La version du code n'est enregistrée nulle part.** Proposition : le commit git (et s'il y avait des modifications non commises) dans la configuration résolue, dans `fit.npz`, dans l'en-tête des fichiers corrigés et dans celui du gabarit LBL. ☑ **Fait (aacc6b1)** : `provenance` dans la configuration résolue, `pca2d_code` dans `fit.npz`, `PCA2GIT` dans les fichiers corrigés et le gabarit LBL.
4. **Les anciennes orthographes** (`star_resolution`, `highpass.window`, `--max-mad`) restent lues pour refaire les anciens passages, mais ne sont plus écrites ni documentées. ☐
5. **Un gabarit qui change reçoit un nouveau nom d'objet LBL.** On n'écrase ni masques ni vitesses. ☑ **Fait**, par le nom de variante.
6. **La référence TOI-2120 est à refaire avec exactement le nominal adopté.** L'actuelle (spl_shr) a été faite avec 6 itérations en gardant la dernière, depuis une configuration temporaire. ☐
7. **Pas de modification du code pendant qu'un passage tourne.** Les étapes importent le code au fil de l'eau : un passage peut donc mélanger deux versions alors que sa provenance n'en nomme qu'une. Proposition : le pipeline revérifie le commit à la fin, et le signale s'il a changé. ☐

| Clé | défaut du code | `config.yaml` (nominal) |
|---|---|---|
| composantes | 2, 7 | 1, 3 |
| `mean` | offset | star |
| stockage | float64 | float32 |
| `highpass.mode` | divide | log_sub |
| `registration.frame` | barycentric | observer |
| `max_nan_fraction` | 0,5 | 0,9 |
| `max_sky_ratio` | aucun | 4 |
| `isolated_window` | aucun | 3 |
| `empirical_noise_box` | aucun | 201 |
| `max_memory_gb` | 12 | 16 |
| `lbl.run` | false | true |
| arrêt | 16, 2, meilleure | 16, (2, meilleure) ; mais les mesures TOI-2120 : 6, 99, dernière |

---

## 2 bis. Résultat du 11 septembre au soir : une garde sur l'étendue en BERV

Sur TOI‑4552, deux tranches de 14 nuits de la même campagne, au même SNR médian, ne diffèrent que par leur étendue en BERV, et la correction change de signe :

| tranche | étendue BERV | livré (rms / robuste) | corrigé (rms / robuste) |
|---|---:|---|---|
| bervwide | 42,8 km/s | 7,80 / 8,48 | 6,18 / 4,52 |
| bervnarrow | 0,7 km/s | 8,96 / 9,14 | 21,76 / 21,76 |

Les deux saisons confirment : s1 (52,6 km/s) passe de 12,28 à 8,63 en robuste, s2 (9,8 km/s) se dégrade de 11,19 à 14,93. Écartés comme causes : l'empilement par nuit, le nombre de composantes, le rétrécissement, le gabarit et le masque.

**À ajouter au nominal :** une garde sur l'étendue en BERV de l'ensemble ajusté, qui refuse ou signale une correction sous un seuil à choisir. ☐

## 2 ter. Résultat du 12 septembre : trois choses mesurées, une à décider

**Retirer la composante stellaire aide partout où c'est mesuré** (`n_star: 0`, le
passage étiqueté `0-3_k0`). Une amplitude libre devant un terme en logarithme est
un exposant sur le flux, et un exposant sur le spectre moyen d'une étoile ne
décrit aucune étoile :

| cible | livré (rms / robuste) | 1-3 | **0-3 (k0)** |
|---|---|---|---|
| TOI-2120 | 47,71 / 34,99 | 18,34 / 14,37 | **17,40 / 14,07** |
| GJ 1 | 2,70 / 2,74 | 2,85 / 2,78 | **2,61 / 2,38** |
| TOI-4552 | 12,64 / 11,37 | 15,86 / 15,86 | **14,24 / 11,20** |

GJ 1 passe ainsi *sous* la série livrée sur les deux mesures, et TOI-4552 sous
elle en robuste. Le k0 de PROXIMA est dans la chaîne en cours.

**Le gabarit du LBL, pas le nôtre** (`lbl.star_template: false`). Sur Proxima, le
nôtre faisait ajuster au LBL des raies 20 % plus larges et doublait l'erreur par
pose, de 0,97 à 1,67 m/s. ☑ adopté

**Le masque commun** (`correct.mask: common`). Un jeu de raies par campagne et non
par époque : sur TOI-4552 le masque par pose coûtait à lui seul 1,8 m/s. ☑ adopté

**Le fit conjoint exige des nuits partagées, pas des SNR semblables.** Mesuré sur
PROXIMA+GJ1+GJ3090 : 2 nuits communes aux trois, sur des campagnes de 258, 147 et
99 nuits. PROXIMA y gagne (2,91 / 2,42 contre 3,04 / 2,52 en solo) parce qu'il
pèse 773 des 1453 lignes et que la base partagée est la sienne, simplement
empêchée de le suivre ; **GJ 1 y perd ce que son solo avait gagné** (3,04 / 2,85
contre 2,61 / 2,38). L'expérience qui tranche est écrite en deux variantes,
`nightsshared` (les 47 nuits que GJ 1 et GJ 3090 partagent) et `sameseason` (la
même période, pas les mêmes nuits), non lancées. ☐

### Options ajoutées le 12 septembre

- `quality.nights` : une **liste** de nuits (rjd entiers) là où `min_rjd`/`max_rjd`
  est une plage. Les nuits que deux campagnes partagent ne sont pas contiguës.
  Hachée dans la clé du cube quand elle est définie.
- **`pca2d-gui`** : la fenêtre qui lance le tout, avec l'index par dossier de
  données (SNR médian, pose médiane, magnitude et sa bande), la coche par cible,
  le refus de deux instruments, le compte des nuits partagées et la fenêtre LBL.
  `docs/gui.md`.
- Une **variante dans un passage conjoint** écrit désormais sous son propre nom
  (`outputs/_<variante>/joint/...`). Avant, deux variantes conjointes partageaient
  un dossier science du LBL, qui aurait mesuré le mélange sans un mot.

## 3. À décider par vous

0. **Le nominal de `config.yaml` : `twoframe.n_star` à 0 ?** Le fichier dit encore
   `1`, et tous les passages k0 passent par le drapeau ou la variante. Le mettre à
   `0` ferait du nominal ce qui est mesuré comme meilleur partout, et ne périme
   aucun cube (`twoframe` n'est pas dans la clé). Je ne l'ai pas fait : c'est votre
   nominal. ☐


1. **Règle d'arrêt.** Sur TOI-2120 (spline), le χ² baisse à chaque itération. Sur l'ancien passage Proxima, R² plafonne à l'itération 2 (0,9861) puis redescend (0,983 à l'itération 11). Je propose **16 itérations au plus, arrêt après 2 itérations moins bonnes, garder la meilleure**. Cela vaut « la dernière » quand le χ² baisse toujours, et protège quand il remonte ; c'est aussi ce que le passage Proxima en cours utilise. ☐
2. **Lissage du gabarit stellaire.** Aucun pour l'instant. On choisit une largeur après ss05, ss025s et un essai sur Proxima : à SNR élevé, un lissage peut coûter plus. ☐
3. **Pas de la grille :** 0,5 km/s ou « smart ». NIRPS n'a ni la résolution ni les pixels de SPIRou. ☐
4. **NaN sur les résidus :** garder désactivé. ☐

---

## 4. Plus tard, quand on aura compris : candidats au retrait

Rien n'est retiré maintenant. Ces voies ne servent pas au nominal :
- **Partie statique :** `offset`, `full`, `iterate`.
- **Représentation et décalage :** base grille de l'étoile, décalage par FFT, ancien gabarit médian unique, `star_resolution`.
- **Correction :** significativité lissée, lissage de composantes observateur.
- **Entrée et poids :** l'entrée s1d (et `grid_source: native`, `highpass.frame`), les modes de transmission autres que Recon, les poids uniformes, `ln_clip_high`, `sigma_floor_frac`, `snr_pixel_scale`, `input.object`, `input.orders`.
- **Divers :** `--instrument`, et `strpca` (utile seulement avec au moins 2 composantes étoile).
