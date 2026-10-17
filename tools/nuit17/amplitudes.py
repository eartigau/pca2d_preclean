#!/usr/bin/env python
"""Ce que les en-tetes des fichiers corriges disent du bruit ajoute.

    python outputs/reports/nuit17/amplitudes.py --run outputs/GJ1_N17/1-3 \
        --target GJ1_N17 --corrected-object GJ1_N17_PCA2D_1-3

Deux questions, sans nouveau run:

1. LES COUPURES ATTEIGNENT-ELLES TOUTES LES POSES D'UNE NUIT? Sur un cube
   empile, la ligne d'une nuit porte le nom de sa premiere pose, et l'etape
   correct cherche les coupures du residu par nom de fichier. PCA2RNAN (NaN
   des coupures) et PCA2CCOL (colonnes retirees) par pose, premiere pose de la
   nuit contre les autres.

2. LA DISPERSION INTRA-NUIT AJOUTEE SUIT-ELLE LES AMPLITUDES? Chaque pose
   est corrigee avec ses propres amplitudes b, reajustees contre la base. Si
   ce reajustement met du bruit dans b, la vitesse ajoutee (corrigee moins
   livree) doit, A L'INTERIEUR D'UNE NUIT, suivre les ecarts de b a la moyenne
   de la nuit. Regression des ecarts intra-nuit de la vitesse ajoutee sur ceux
   des b: la part de variance expliquee.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits

REPO = "/Users/eartigau/GitHubProjects/pca2d_preclean"
sys.path.insert(0, REPO)
from pca2d.lblscan import rdb_rows          # noqa: E402
from pca2d.logger import log                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(REPO, "outputs", "lbl", "lblrdb")


def stem(name):
    """Le nom de la pose livree, quel que soit le suffixe: tout avant 't'."""
    base = os.path.basename(str(name))
    return base.split("t_")[0].split("t.fits")[0]


def series(obj):
    path = os.path.join(TREE, "lbl_%s_%s.rdb" % (obj, obj))
    t, v, e, table = rdb_rows(path)
    names = {c.lower(): c for c in table.colnames}
    files = [stem(f) for f in table[names["local_file_name"]]]
    return {s: (ti, vi, ei) for s, ti, vi, ei in zip(files, t, v, e)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--corrected-object", required=True)
    args = p.parse_args(argv)

    heads = {}
    for path in sorted(glob.glob(os.path.join(REPO, args.run, "corrected", "*.fits"))):
        h = fits.getheader(path, 0)
        n = int(h.get("PCAOBS_N", 0))
        heads[stem(path)] = {
            "b": [float(h["PCAOBS%02d" % (j + 1)]) for j in range(n)],
            "wnan": int(h.get("PCA2WNAN", 0)), "rnan": int(h.get("PCA2RNAN", 0)),
            "ccol": int(h.get("PCA2CCOL", 0)), "npix": int(h.get("PCA2NPIX", 0))}
    deliv = series(args.target)
    corr = series(args.corrected_object)
    keys = sorted(set(heads) & set(deliv) & set(corr), key=lambda k: deliv[k][0])
    log("%s: %d exposures with a header, a delivered and a corrected velocity"
        % (args.corrected_object, len(keys)), "value")
    t = np.array([deliv[k][0] for k in keys])
    added = np.array([corr[k][1] - deliv[k][1] for k in keys])
    added -= np.median(added)
    b = np.array([heads[k]["b"] for k in keys])
    rnan = np.array([heads[k]["rnan"] for k in keys])
    wnan = np.array([heads[k]["wnan"] for k in keys])
    night = np.floor(t - 0.5).astype(int)

    # 1. qui recoit les coupures
    first = np.zeros(len(keys), dtype=bool)
    for n in np.unique(night):
        idx = np.flatnonzero(night == n)
        first[idx[0]] = True
    multi = np.array([np.sum(night == n) > 1 for n in night])
    out = {"n": len(keys), "nights": int(np.unique(night).size)}
    for name, sel in (("first", first & multi), ("others", ~first & multi)):
        out["rnan_" + name] = float(np.median(rnan[sel])) if sel.any() else np.nan
        out["rnan_zero_" + name] = float(np.mean(rnan[sel] == 0)) if sel.any() else np.nan
        out["wnan_" + name] = float(np.median(wnan[sel])) if sel.any() else np.nan
    log("residual-cut NaN per exposure, nights of 2+: first exposure median %s"
        " (%.0f%% at zero), the others %s (%.0f%% at zero)"
        % (out["rnan_first"], 100 * out["rnan_zero_first"], out["rnan_others"],
           100 * out["rnan_zero_others"]), "value")

    # 2. la vitesse ajoutee suit-elle les b, dans la nuit
    dv, db = [], []
    for n in np.unique(night):
        s = night == n
        if s.sum() < 2:
            continue
        dv.append(added[s] - added[s].mean())
        db.append(b[s] - b[s].mean(axis=0))
    dv = np.concatenate(dv)
    db = np.vstack(db)
    coef, *_ = np.linalg.lstsq(db, dv, rcond=None)
    explained = 1 - np.var(dv - db @ coef) / np.var(dv)
    out.update(intra_added_rms=float(np.std(dv)), intra_explained=float(explained),
               intra_coef=[float(c) for c in coef],
               intra_b_rms=[float(x) for x in db.std(axis=0)])
    log("intra-night added velocity: rms %.2f m/s, of which %.0f%% of the variance"
        " follows the exposures' own b (coefficients %s m/s per unit b)"
        % (out["intra_added_rms"], 100 * explained,
           ", ".join("%.2f" % c for c in coef)), "value")
    path = os.path.join(HERE, "amplitudes_%s.json" % args.corrected_object)
    with open(path, "w") as handle:
        json.dump(out, handle, indent=1)
    log("wrote %s" % path, "info")


if __name__ == "__main__":
    main()
