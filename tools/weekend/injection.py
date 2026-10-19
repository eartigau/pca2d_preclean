#!/usr/bin/env python
"""Ce que la correction laisse d'une keplerienne injectee.

    python tools/weekend/injection.py TOI210 [--orig TOI210_W] [--joint]

Pose par pose: (serie injectee - serie originale) contre la vitesse injectee
(injected.csv, a cote des fichiers injectes). Sur les spectres livres la pente
doit valoir 1 (le controle de l'outil d'injection); sur les corriges, elle est
la fraction du signal que la correction conserve. Avec --joint, les series
corrigees sont celles de la base commune.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "outputs", "reports", "nuit17"))
import analyse as A                                    # noqa: E402

A.TREE = "/Volumes/irrisor/pca2d_preclean/lbl/lblrdb"


def injected(target):
    path = "/Volumes/irrisor/pca2d_preclean/science/%s_INJ/injected.csv" % target
    out = {}
    for row in csv.DictReader(open(path)):
        out[round(float(row["bjd"]) - 2400000.0, 5)] = float(row["v_inj_ms"])
    return out


def series(name, template=None):
    run = A.load(name, template or name)
    if run is None:
        raise SystemExit("pas de serie LBL pour %s" % name)
    return {round(float(t), 5): (float(v), float(e)) for t, v, e in zip(run["t"], run["v"], run["e"])}


def slope(a, b, vinj):
    """pente de (b - a) contre la vitesse injectee, et son erreur."""
    keys = sorted(set(a) & set(b) & set(vinj))
    d = np.array([b[k][0] - a[k][0] for k in keys])
    x = np.array([vinj[k] for k in keys])
    d = d - np.median(d)
    A_ = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A_, d, rcond=None)
    resid = d - A_ @ coef
    err = np.std(resid) / np.sqrt(np.sum((x - x.mean()) ** 2))
    return coef[0], err, np.std(resid), len(keys)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("target", help="TOI210 ou TOI782")
    p.add_argument("--orig", help="l'objet original (TOI210_W, TOI782_F)")
    p.add_argument("--joint", action="store_true")
    args = p.parse_args()
    orig = args.orig or {"TOI210": "TOI210_W", "TOI782": "TOI782_F"}[args.target]
    inj = "%s_INJ" % args.target
    vinj = injected(args.target)
    s, e, r, n = slope(series(orig), series(inj), vinj)
    print("%s livre:  pente %.3f +/- %.3f (residu %.2f m/s, %d poses)  <- doit valoir 1"
          % (args.target, s, e, r, n))
    suffix = "_PCA2D_1-3_joint" if args.joint else "_PCA2D_1-3"
    s, e, r, n = slope(series(orig + suffix), series(inj + suffix), vinj)
    print("%s corrige%s: pente %.3f +/- %.3f (residu %.2f m/s, %d poses)  <- fraction conservee"
          % (args.target, " (base commune)" if args.joint else "", s, e, r, n))


if __name__ == "__main__":
    main()
