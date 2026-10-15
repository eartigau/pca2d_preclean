# Le modèle TAPAS est décalé de tout le BERV sur les S2D de l'ESO

Mesuré le 2026-10-15 sur HD 10700, bande B de l'oxygène (686-690 nm).

## La mesure

En corrélant le modèle TAPAS avec le spectre sur une grille log commune, et en
cherchant le décalage qui maximise l'accord :

| grille lue | accord | décalage à appliquer |
|---|---|---|
| `WAVEDATA_VAC_BARY` (vide) | **0,974** | −23,79 km/s |
| `WAVEDATA_AIR_BARY` (air) | 0,431 | +17,49 km/s |

**La grille du vide est la bonne** : l'accord y est de 0,974 contre 0,431.
C'est cohérent avec LBL, dont la classe `Harps_ESO` lit l'extension 4, qui est
précisément `WAVEDATA_VAC_BARY`, et ne fait qu'une conversion d'unités. Les
classes `Harps` et `Harps_ORIG`, elles, lisent une solution **en air** et la
convertissent (`harps.py:600` et `1085`, « HARPS wave solution is in air -
convert to vacuum », avec `mp.air_index`). Chaque classe fait donc ce qu'il
faut pour SA source ; il n'y a pas de confusion air/vide dans le chemin ESO.

Mais il reste un décalage, et sur 10 poses il est sans ambiguïté :

| | |
|---|---|
| pente du décalage contre le BERV | **−1,0003** |
| corrélation | **−1,0000** |
| ordonnée à l'origine | −0,03 km/s |

## La cause

`WAVEDATA_VAC_BARY` est dans le vide **et barycentrique**. TAPAS décrit
l'atmosphère terrestre, donc le repère du sol. Il faut l'évaluer à
`lambda / D(BERV)`, et personne ne le fait :

- `lbl_telluclean.py:179` passe `e2ds_dict['BERV'] = inst.get_berv(sci_hdr)` ;
- le `get_berv` de `Harps_ESO` rend **zéro**, parce que la grille est déjà
  barycentrique et que c'est la bonne réponse pour les vitesses ;
- `tellu_clean.py:398` écrit `wave_trans = np.array(wave_vector)` sous le
  commentaire « get a wave solution for the transmission shifted by the berv
  if required », mais ne décale rien.

Le modèle est donc évalué aux longueurs d'onde barycentriques, décalé de tout
le BERV, jusqu'à 27 km/s sur des raies de 2,6 km/s de large.

Vérifié : en évaluant TAPAS à `lambda / D(BERV)`, le décalage résiduel tombe à
**0,00 km/s**, la pente contre le BERV à **0,0000**, et l'accord monte à 0,955.

## Ce que cela explique

- **« on dirait que les telluriques n'ont pas été corrigées »** : la division
  appliquée plafonne à ×2,7 là où le modèle demande ×7,5, parce qu'elle tombe
  à côté des raies ;
- **« à 593-595 c'est tout beau »** : la transmission y reste vers 0,77, et un
  décalage y coûte beaucoup moins qu'au fond d'une bande saturée ;
- et le résidu qui ne descend pas dans la bande B quoi qu'on fasse.

## Qui est touché

Tout instrument dont la grille livrée est barycentrique et dont le `get_berv`
rend zéro : **HARPS/ESO, HARPS-N/ESO, ESPRESSO**. SPIRou et NIRPS ne le sont
pas, leur grille étant dans le repère de l'observateur.

C'est en amont de pca2d : la correction tellurique de LBL est touchée pour
elle-même, que pca2d tourne ensuite ou non.

## Vérification par le graphique

`/Users/eartigau/tmp/diagnostic_harps/abso_check.pdf` trace, pour quatre
fenêtres, le spectre brut, le spectre nettoyé, et le modèle TAPAS des deux
façons : évalué dans le repère du sol (`lambda / D(BERV)`) et sur la grille
barycentrique, c'est-à-dire ce que LBL fait aujourd'hui.

Accord entre le modèle et le spectre, par fenêtre :

| fenêtre | repère du sol | grille barycentrique | part tellurique du signal |
|---|---|---|---|
| **O2 bande B** (686-690) | **0,973** | **−0,035** | 1,12 |
| O2 gamma (626-630) | 0,282 | −0,001 | 0,55 |
| eau (593-595) | 0,183 | 0,035 | 0,26 |
| eau (646-652) | 0,087 | −0,008 | 0,15 |

Dans la bande B, l'accord passe de **−0,035**, c'est-à-dire aucune corrélation,
à **0,973**. Le modèle était simplement ailleurs.

Les accords faibles des trois autres fenêtres ne cachent pas un second
problème : ils suivent la **part tellurique du signal**, donnée dans la
dernière colonne comme le rapport entre la profondeur du modèle et celle du
spectre. Là où les telluriques ne font que 15 % de la structure (eau à
649 nm), les raies stellaires dominent la corrélation et la noient. L'ordre
des quatre valeurs est le même dans les deux colonnes, ce qui est le
comportement attendu.
