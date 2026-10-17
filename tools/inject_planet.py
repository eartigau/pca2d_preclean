#!/usr/bin/env python
"""Injecter une ou plusieurs keplériennes circulaires dans une campagne t.fits.

    python tools/inject_planet.py --source data/GJ1 \
        --dest /Volumes/irrisor/pca2d_preclean/science/GJ1_INJ17 \
        --planet 37.3:3.0 --planet 365.25:3.0 --seed 17

Chaque fichier est recopie avec sa solution en longueur d'onde (WaveA)
multipliee par le facteur Doppler RELATIVISTE de la vitesse injectee a son BJD,
sqrt((1 + v/c) / (1 - v/c)), v positive en eloignement. Une raie qui etait au
pixel p sort donc a la longueur d'onde lambda_p * D(v): l'etoile a pris la
vitesse v. Rien d'autre n'est touche, ni le flux, ni le blaze, ni la
transmission.

LA LIMITE, dite une fois: decaler la solution decale TOUT ce que le fichier
porte, les residus telluriques et le ciel compris, et pas seulement l'etoile.
A quelques m/s contre des raies telluriques larges de quelques km/s, l'effet
sur la composante observateur est de l'ordre de 1e-3 de son amplitude; c'est
la facon habituelle d'injecter, et c'est la seule sans separer l'etoile de
l'atmosphere, ce qui est precisement la question posee.

Le BJD est lu par le lecteur du pipeline (tfits.read_tfits), donc avec la
convention exacte du cube. Les vitesses injectees sont ecrites dans
injected.csv a cote des fichiers, et sur chaque en-tete primaire (INJ_V,
INJ_DESC). Un fichier deja present a la destination n'est jamais remplace.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from pca2d import tfits as sptf        # noqa: E402
from pca2d.logger import log           # noqa: E402

C_MS = 299792458.0


def doppler(v_ms):
    """Le facteur relativiste d'une vitesse en m/s, positive en eloignement."""
    beta = np.asarray(v_ms, dtype=float) / C_MS
    return np.sqrt((1.0 + beta) / (1.0 - beta))


def planets_velocity(bjd, planets):
    """Somme de sinusoides: planets = [(P jours, K m/s, phase rad), ...]."""
    bjd = np.asarray(bjd, dtype=float)
    v = np.zeros_like(bjd)
    for period, amplitude, phase in planets:
        v += amplitude * np.sin(2 * np.pi * bjd / period + phase)
    return v


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", required=True, help="le dossier des t.fits livres")
    p.add_argument("--dest", required=True, help="le dossier a creer")
    p.add_argument("--planet", action="append", required=True,
                   help="P:K, periode en jours et demi-amplitude en m/s;"
                        " repetable")
    p.add_argument("--seed", type=int, default=17,
                   help="la graine des phases, tirees uniformement")
    p.add_argument("--wave-ext", default="WaveA",
                   help="l'extension de la solution en longueur d'onde")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rng = np.random.default_rng(args.seed)
    planets = []
    for spec in args.planet:
        period, amplitude = (float(x) for x in spec.split(":"))
        planets.append((period, amplitude, float(rng.uniform(0, 2 * np.pi))))
    files = sorted(glob.glob(os.path.join(args.source, "*t.fits")))
    if not files:
        raise SystemExit("aucun t.fits dans %s" % args.source)
    os.makedirs(args.dest, exist_ok=True)
    desc = " + ".join("P=%.4gd K=%.3gm/s phi=%.3f" % p for p in planets)
    log("injecting %s into %d files, %s -> %s"
        % (desc, len(files), args.source, args.dest), "info")
    table = os.path.join(args.dest, "injected.csv")
    rows, skipped = [], 0
    for i, path in enumerate(files):
        out = os.path.join(args.dest, os.path.basename(path))
        if os.path.exists(out):
            skipped += 1
            continue
        bjd = float(sptf.read_tfits(path)["meta"]["bjd"])
        v = float(planets_velocity([bjd], planets)[0])
        with fits.open(path) as hdul:
            hdul[args.wave_ext].data = (hdul[args.wave_ext].data.astype(np.float64)
                                        * doppler(v))
            hdul[0].header["INJ_V"] = (v, "injected stellar velocity, m/s")
            hdul[0].header["INJ_DESC"] = (desc[:68], "inject_planet.py")
            hdul.writeto(out)
        rows.append((os.path.basename(path), bjd, v))
        if (i + 1) % 50 == 0:
            log("  %d of %d written" % (i + 1, len(files)), "info")
    new = not os.path.exists(table)
    with open(table, "a", newline="") as handle:
        writer = csv.writer(handle)
        if new:
            writer.writerow(["filename", "bjd", "v_inj_ms"])
        writer.writerows(rows)
    log("wrote %d files, %d already there and left alone; velocities in %s"
        % (len(rows), skipped, table), "value")
    if rows:
        v = np.array([r[2] for r in rows])
        log("injected velocity: rms %.3f m/s, from %.3f to %.3f"
            % (np.std(v), v.min(), v.max()), "value")
    return None


if __name__ == "__main__":
    main()
