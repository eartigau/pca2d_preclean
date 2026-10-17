# La nuit du 2026-10-17: pourquoi Proxima et GJ 1 faisaient moins bien

Le plan: `plan_nuit_proxima_gj1.md`. Le rapport (12 pages, 12 Mo, non versionne):
`outputs/reports/nuit17/rapport_nuit17.pdf`, reconstruit par `report.py`.

- `runner.sh`, `queue.txt`: le lanceur et sa file, telle qu'elle a fini la nuit
  (relue avant chaque tache, modifiee en cours de route; chaque modification
  est commentee dans le fichier). `done.txt`: chaque tache, son code de sortie,
  ses heures, sa version. `VERSIONS.txt`: le gel.
- `analyse.py`: toutes les series LBL contre la livree, sur les memes poses;
  `metrics_*.json` en sortie.
- `amplitudes.py`: les en-tetes des fichiers corriges (coupures par pose, et la
  part de la vitesse ajoutee intra-nuit qui suit les amplitudes b).
- `report.py` et `findings.tex`: le rapport.

Deux defauts du lanceur, a corriger avant de le reprendre: la derniere ligne de
la file est ignoree si elle n'a pas de saut de ligne (G09L lance a la main), et
`caffeinate` s'arrete avec le lanceur, donc rien ne garde la machine eveillee
apres.
