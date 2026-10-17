# Plan de nuit : pourquoi Proxima et GJ 1 font moins bien

Préparé le 2026-10-16 au soir, à discuter avant le lancement.

## 1. Ce qu'on sait déjà (ne pas le refaire)

Tout ce qui suit date de septembre, sur l'ancien code (avant les correctifs
des 15 et 16 octobre).

| | Livré | Corrigé 1-3 | Corrigé 0-3 (k0) |
|---|---|---|---|
| Proxima, rms / robuste (765 poses) | 2.86 / 2.22 | 3.34 / 2.64 | 2.93 / 2.37 |
| Proxima, erreur médiane | 0.92 | **1.67** | ? |
| GJ 1, rms / robuste (solo) | 2.70 / 2.74 | ~2.56 (par nuit) | **2.61 / 2.38** |

- **L'envergure BERV n'est pas en cause** : 42.6 km/s pour Proxima, 48.1
  pour GJ 1. Le mode d'échec connu (TOI-4552, 0.7 km/s) ne s'applique pas.
- **L'injection en forme de BERV a déjà disparu** sur Proxima (Spearman -0.22
  entre la différence injectée et le BERV) : ce qui reste est un coût, pas un
  biais.
- **Les barres d'erreur presque doublent** sur Proxima (0.92 vers 1.67). En
  quadrature, cela fait 1.39 m/s ajoutés, contre 1.73 m/s ajoutés au rms : la
  hausse des erreurs explique **à peu près les deux tiers de la variance
  ajoutée**. Ordre de grandeur à refaire proprement, mais le gros de la perte
  ressemble à du bruit injecté, pas à un biais.
- **Sans composante stellaire (0-3), les deux cibles font mieux** qu'avec
  (1-3) : Proxima à peu près à égalité avec le livré, GJ 1 gagne en robuste.
- **Le template donné à LBL coûte cher sur Proxima** : sur les spectres
  *livrés*, mesurer avec le template de pca2d au lieu de celui de LBL faisait
  passer le rms de 2.98 à 3.97 m/s (attribution du 2026-09-11, ancien fit).
  Sur TOI 2120, le template de pca2d gonflait les erreurs par raie de 20 %,
  « cause toujours ouverte ». Sur TOI-4552 (faible SNR), c'était l'inverse :
  le template de pca2d faisait mieux (15.86 contre 16.82).

## 2. Ce que HARPS change à la lecture

- **Trois correctifs touchent les deux cibles et n'y ont jamais été
  appliqués** : les fonds de raies profondes ne sont plus coupés (7b88be6,
  sur HD 10700 l'ancien seuil tuait 13.7 % des pixels et un tiers des raies
  LBL), le passe-haut contre template (masques qui se déplacent), et la
  classe HE de LBL pour NIRPS HE + CADC (lbl-rapide 95bf4f4, hier soir). Tous
  les produits NIRPS de septembre ont tourné avec le profil HA.
- **Le SNR seul n'ordonne pas le gain.** Sur HARPS, GJ 581 à SNR 13 n'a
  gagné que 6 %, là où la loi de puissance de la présentation lui prédisait
  plus d'un facteur 2. La présentation (diapos 13 et 14) devra être revue
  après la nuit.
- **Le gonflement des barres d'erreur croît avec le SNR**, sur le code
  actuel aussi : GJ 581 (SNR 13) +8 %, TOI 2120 (SNR 33) +8 %, HD 69830
  (SNR 72) +39 %, HD 10700 (SNR 90) +29 %, et Proxima (SNR 287, ancien code)
  +82 %. C'est ce que ferait un bruit injecté d'amplitude fractionnaire fixe :
  négligeable devant le bruit de photons d'une cible faible, dominant devant
  celui d'une cible brillante.

## 3. Hypothèses, par ordre de vraisemblance

**H1. Le template.** Le template de pca2d (P1 à amplitude moyenne, sur la
grille de 0.5 km/s) est meilleur que celui de LBL quand les spectres sont
bruités et moins bon quand ils ne le sont pas. Comme LBL tire son bruit par
pixel du résidu au template (`USE_NOISE_MODEL` à False), un template imparfait
gonfle directement les erreurs. *Prédit* un croisement en SNR, à lui seul.
*Test* : la matrice 2 x 2 {livré, corrigé} x {template LBL, template pca2d},
qui ne demande que des passes LBL.

**H2. La version.** Les trois correctifs ci-dessus. *Test* : tout refaire
sur une version figée, ce qui est de toute façon obligatoire.

**H3. La composante stellaire.** 0-3 fait mieux que 1-3 sur les deux cibles.
Mais avec 0-3 le template LBL change aussi, donc H3 peut se réduire à H1.
La matrice le dira.

**H4. Un plancher de bruit injecté par b.Q.** Q estimé sur l'étoile elle-même
porte son propre bruit et un peu de l'étoile ; divisé dans chaque pose, il
ajoute un motif fixe dont l'effet relatif croît avec le SNR. *Tests* :
mesurer directement le bruit par pixel corrigé contre livré sur les colonnes
sans telluriques, en fonction du SNR de la pose, et en tirer le plancher
fractionnaire.

**H5. La métrique.** Le rms de Proxima (2.86 m/s pour 0.92 d'erreur) est
dominé par l'étoile et ses planètes. Retirer 1 m/s de systématique d'un rms de
2.9 ne le ramène qu'à 2.7 : le rapport des rms est presque aveugle chez une
cible brillante, alors que le coût en bruit s'y voit en entier.

**H6. Les données.** Mélange de modes HA et HE, jamais vérifié sur ces
deux cibles (la leçon HAM et EGGS de HARPS). Et les éruptions de Proxima,
qui peuvent tirer une composante à elles.

## 4. Hors des sentiers battus

- **Correction nulle.** Passer par l'étape `correct` avec b = 0 : les fichiers
  doivent ressortir identiques aux livrés (masque NaN mis à part), et LBL doit
  redonner le rms livré. Si ce n'est pas le cas, c'est la chaîne qui injecte,
  pas la correction. C'est exactement la classe de bogue trouvée aujourd'hui
  dans `reconstruct.py`.
- **b permutés.** Appliquer à la nuit i la correction de la nuit j. Même Q,
  même distribution d'amplitudes, mais sans rapport avec la nuit. Si le rms
  sort pareil qu'avec la vraie correction, celle-ci n'apporte rien sur cette
  cible.
- **Base validée croisée.** Ajuster Q sur les nuits paires, corriger les
  impaires, et l'inverse. L'écart avec la correction « dans l'échantillon »
  mesure le surajustement.
- **Correction restreinte au support tellurique.** Q mis à zéro là où la
  transmission (TAPAS, sur toute la plage de BERV) dépasse 0.995. Si le bruit
  injecté vient des colonnes propres, cela le coupe net. Attention : sur TOI
  2120 la composante 1 était significative même hors telluriques (un motif
  instrumental), donc ce test peut aussi perdre quelque chose.
- **Injection de planète.** Décaler (en relativiste) la solution en longueur
  d'onde des fichiers livrés d'une keplérienne connue, une à 35 jours et une à
  365 jours, puis refaire toute la chaîne. Les deux équipes demanderont si le
  signal est préservé, surtout à un an. Plutôt sur GJ 1, plus petit et calme.
- **NIRPS traité comme HARPS.** Partir des spectres *non* corrigés des
  telluriques, faire le nettoyage LBL par la dérivée avec template, puis
  pca2d. Compare la correction d'APERO à la nôtre sur la même étoile.
  Seulement si les fichiers non corrigés sont sur irrisor.
- **Une base d'étoiles chaudes.** Ajuster Q sur les standards telluriques de
  NIRPS (pas de raies stellaires, donc rien de l'étoile ne peut fuir dans Q),
  et seulement b sur Proxima. La machinerie `--objects` le permet presque.
  Seulement si les standards sont disponibles.

## 5. Mesures qui tranchent (pas seulement le rms)

- rms et robuste, sur les poses communes, bootstrap apparié (comme d'habitude)
- **dispersion intra-nuit** : différences entre poses d'une même nuit, sur
  laquelle ni les planètes ni l'activité n'ont le temps de bouger (Proxima a
  environ 3 poses par nuit, GJ 1 environ 2)
- **rms après retrait de Proxima b et d**, éphémérides tirées de la NASA
  Exoplanet Archive pendant la nuit, pas de mémoire
- **le biais en V_tot** (la figure de ce soir), avant et après
- **les barres d'erreur**, et le plancher de bruit par pixel (H4)
- **par raie** (fichiers `lblrv`) : variation du rms de chaque raie contre la
  profondeur tellurique maximale qu'elle croise, sa propre profondeur, la bande
  Y/J/H. Distingue un effet global (template, passe-haut) d'un effet local
  (telluriques).

## 6. Déroulé et budget

Machine : 16 Go, 10 coeurs. Proxima et GJ 1 en séquence, pas en parallèle.
Durées estimées depuis TOI 2120 (321 poses, 44 min) ; je mesure GJ 1 d'abord
et je replanifie.

| Étape | Contenu | Durée |
|---|---|---|
| 0 | irrisor monté, `caffeinate`, version figée (hash pca2d et lbl-rapide, refus si non commité), `VERSIONS.txt`, inventaire HA/HE par cible | 15 min |
| 1 | GJ 1 : cube, fits 1-3 et 0-3, puis la matrice 2 x 2 de templates pour chacun | 1 h 30 |
| 2 | Proxima : idem | 2 h 30 à 3 h |
| 3 | Diagnostics sur ces produits (section 5), sans nouveau fit | 1 h |
| 4 | Variantes choisies selon l'étape 3 (arbre ci-dessous) | 3 h |
| 5 | Rapport LaTeX, ouvert à la fin | 30 min |

**Arbre de décision après l'étape 3** :

- si la matrice dit que le template porte la perte : variantes sur le
  template (LBL sur les corrigés avec son propre template, template pca2d
  sans passer par la grille de 0.5 km/s) ;
- si le livré avec b = 0 ne redonne pas le livré : chercher le bogue dans la
  chaîne, tout le reste attend ;
- si le plancher de bruit vient des colonnes sans telluriques : support
  tellurique, rétrécissement plus fort, moins de composantes (n_earth 1 et 2) ;
- s'il vient de l'estimation de Q : base validée croisée, base d'étoiles
  chaudes ;
- dans tous les cas, s'il reste du temps : b permutés et injection de planète
  sur GJ 1.

## 7. Garde-fous

- une seule version pour toute la nuit ; **toute option de diagnostic doit
  être codée, testée et commitée avant le gel**, sinon il faut tout refaire ;
- les produits NIRPS de septembre mis de côté, jamais effacés, et chaque
  variante sous un nouveau nom d'objet (pas de `rm`, refusé deux fois ce soir) ;
- un mode NIRPS par cible, le majoritaire ;
- toute vitesse par LBL, jamais par projection ;
- ne jamais comparer une variante à un nominal dépassé ;
- journal horodaté avec le code de couleurs, figures en PDF ;
- rapport du matin : la compilation des séries RV contre l'original, variantes
  en couleur, les river plots des deux cibles, et un verdict par hypothèse.

## 8. À décider ensemble

1. **Le gel** : HEAD de pca2d et lbl-rapide avec la classe HE corrigée.
   D'accord pour invalider tous les produits NIRPS de septembre, livrés
   compris ?
2. **Le verdict** : je propose dispersion intra-nuit, rms sans planètes et
   biais en V_tot comme mesures principales, rms et robuste gardés pour la
   continuité.
3. **Le code à écrire avant le gel** (1 h 30 à 2 h ce soir, avec tests) :
   correction nulle (peut-être déjà possible par `correct.n_earth: 0`), b
   permutés, base validée croisée, support tellurique, liste d'éruptions à
   exclure du fit, outil d'injection de planète. Lesquels ?
4. **« Les updates de HARPS » pour Proxima** : je comprends les trois
   correctifs plus le filtre de mode. Vouliez-vous aussi dire la voie HARPS
   complète, spectres non corrigés et nettoyage LBL par la dérivée ? Les
   fichiers non corrigés sont-ils sur irrisor ?
5. **Les standards telluriques NIRPS** : disponibles, et où ?
6. **TOI-1452 et GJ 328** : les refaire aussi, pour que la figure du SNR soit
   sur une seule version (1 h 30 de plus) ? Au passage, TOI-1452 a gagné 1.46
   fois avec une envergure BERV de seulement 3.1 km/s, ce qui contredit la
   règle de l'envergure : à regarder.
