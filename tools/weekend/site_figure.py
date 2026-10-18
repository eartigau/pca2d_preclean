#!/usr/bin/env python
"""La figure du site: dix naines M NIRPS tirees au hasard, le gain contre
l'exces de dispersion du livre et contre le SNR.

    python tools/weekend/site_figure.py [--format svg|pdf]

Lit outputs/reports/weekend/weekend.json (analyse_weekend.py) et ecrit
docs/figures/weekend_gain.<format>. L'exces est la dispersion de la serie
livree au-dela de ses propres barres d'erreur, sqrt(rms^2 - sigma^2) / sigma,
sur les poses que la serie corrigee partage avec elle.
"""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
NAMES = {"GJ1": "GJ 1", "GL699": "Barnard", "YZ_CMI": "YZ CMi", "GL406": "Wolf 359",
         "GJ643": "GJ 643", "GJ3737": "GJ 3737", "GJ3618": "GJ 3618",
         "GJ3622": "GJ 3622", "TOI5524": "TOI-5524", "TOI210": "TOI-210"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--format", default="svg")
    args = p.parse_args()
    t = json.load(open(os.path.join(REPO, "outputs/reports/weekend/weekend.json")))
    rows = []
    for r in t:
        n, c = r["nominal"], r["cutsoff"]
        ref = n["ref"]
        excess = np.sqrt(max(ref["rms"] ** 2 - ref["median_error"] ** 2, 0)) / ref["median_error"]
        rows.append(dict(name=NAMES.get(r["object"], r["object"]), snr=r["snr"],
                         excess=excess, g=n["gain"], glo=n["gain_16_84"][0],
                         ghi=n["gain_16_84"][1], c=c["gain"]))
    x_ex = np.array([r["excess"] for r in rows])
    x_snr = np.array([r["snr"] for r in rows])
    g = np.array([r["g"] for r in rows])
    rho_ex, p_ex = spearmanr(x_ex, g)
    rho_snr, p_snr = spearmanr(x_snr, g)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True)
    for ax, x, label, rho, pv in ((axes[0], x_ex,
                                   "excess scatter of the delivered series,\n"
                                   r"$\sqrt{\mathrm{rms}^2-\sigma^2}\,/\,\sigma$", rho_ex, p_ex),
                                  (axes[1], x_snr, "median SNR per exposure", rho_snr, p_snr)):
        ax.errorbar(x, g, yerr=[g - [r["glo"] for r in rows], [r["ghi"] for r in rows] - g],
                    fmt="o", ms=7, color="#26a1ff", capsize=2.5, lw=1.1, zorder=3)
        # l'amas autour de 1.5 en exces: les noms alternent au-dessus et
        # au-dessous des points, dans l'ordre de x, pour ne pas se chevaucher
        order = np.argsort(x)
        for k, i in enumerate(order):
            up = (k % 2 == 0)
            ax.annotate(rows[i]["name"], (x[i], g[i]), fontsize=7.5,
                        xytext=(4, 7 if up else -11), textcoords="offset points",
                        color="0.25")
        ax.axhline(1.0, color="0.35", lw=0.9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(label)
        ax.set_title(r"Spearman $\rho = %+.2f$ ($p = %.3f$)" % (rho, pv), fontsize=10)
        ax.grid(alpha=0.25, which="both")
    axes[0].set_ylabel("gain in RV rms, delivered / precleaned")
    nul = matplotlib.ticker.NullFormatter()
    plain = matplotlib.ticker.FormatStrFormatter("%g")
    for ax, ticks in ((axes[0], [1, 1.5, 2, 3, 5, 8]), (axes[1], [30, 50, 100, 200])):
        ax.set_ylim(0.8, 3.0)
        ax.set_yticks([0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])
        ax.set_xticks(ticks)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(plain)
            axis.set_minor_formatter(nul)
    fig.suptitle("%d NIRPS M dwarfs drawn at random: the gain follows what the delivered"
                 " series carries beyond its photon noise, not the SNR" % len(rows), fontsize=10.5)
    fig.tight_layout()
    out = os.path.join(REPO, "docs/figures/weekend_gain.%s" % args.format)
    fig.savefig(out, bbox_inches="tight")
    print("wrote %s; excess rho %+.2f p %.4f; SNR rho %+.2f p %.3f" % (out, rho_ex, p_ex, rho_snr, p_snr))


if __name__ == "__main__":
    main()
