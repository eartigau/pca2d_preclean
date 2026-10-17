#!/usr/bin/env python
"""L'analyse de la nuit du 2026-10-17: pourquoi Proxima et GJ 1 font moins bien.

    python outputs/reports/nuit17/analyse.py            # tout ce qui existe
    python outputs/reports/nuit17/analyse.py --target GJ1_N17

Lit les rdb de LBL (arbre outputs/lbl), compare chaque serie a la serie
livree sur les poses communes, et ecrit:

- metrics_<cible>.json, tous les nombres
- compilation_<cible>.pdf, les series RV contre le temps, la livree en noir
  et chaque variante dans sa couleur
- injected_<cible>.pdf, ce que chaque correction a AJOUTE aux vitesses (la
  serie moins la livree) contre V_tot

Les mesures, par serie:
  rms, sigma robuste, erreur mediane, rms des moyennes par nuit
  dispersion INTRA-NUIT (poolee): ni planete ni activite n'y bougent
  rms par nuit APRES retrait des planetes connues (Proxima b et d, periodes de
    config.yaml, amplitude et phase ajustees sur chaque serie)
  biais en V_tot (pca2d.bervbias, meme graine que le rapport)
  le gain en rms contre la livree, bootstrap apparie 2000 tirages
  ce que la correction a ajoute: rms, et Spearman contre le BERV et V_tot

Rien n'est ecrit hors de ce dossier.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import zlib

import numpy as np

REPO = "/Users/eartigau/GitHubProjects/pca2d_preclean"
sys.path.insert(0, REPO)

import matplotlib                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from scipy.stats import spearmanr                   # noqa: E402

from pca2d import bervbias                          # noqa: E402
from pca2d.lblscan import rdb_rows, velocity_stats  # noqa: E402
from pca2d.logger import log                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(REPO, "outputs", "lbl", "lblrdb")

PLANETS = {"PROXIMA_N17": [11.18465, 5.12338]}

# (cle, libelle, science, template) ; {T} est la cible
SERIES = [
    ("delivered", "delivered, LBL template", "{T}", "{T}"),
    # le nominal est mesure contre le template que LBL tire des fichiers
    # corriges (lbl.star_template: false depuis le 2026-09-12); _DFT mesure les
    # spectres LIVRES contre ce meme template
    ("delivered_ft", "delivered, corrected files' template", "{T}_DFT", "{T}_PCA2D_1-3"),
    ("nominal", "pca2d 1-3", "{T}_PCA2D_1-3", "{T}_PCA2D_1-3"),
    ("cutsoff", "1-3, residual cuts off", "{T}_PCA2D_1-3_n17_cutsoff", None),
    ("nonan", "1-3, no NaN added", "{T}_PCA2D_1-3_n17_nonan", None),
    ("k0", "0-3, no star component", "{T}_PCA2D_0-3_n17_k0", None),
    ("k0cut", "0-3, residual cuts off", "{T}_PCA2D_0-3_n17_k0cut", None),
    ("k0c", "0-3, residual cuts off", "{T}_PCA2D_0-3_n17_k0c", None),
    ("permb", "1-3, b of another night", "{T}_PCA2D_1-3_n17_permb", None),
    ("qsupport", "1-3, telluric support only", "{T}_PCA2D_1-3_n17_qsupport", None),
    ("cvodd", "1-3 cuts off, odd nights held out", "{T}_PCA2D_1-3_n17_cvodd", None),
]

COLOURS = {"delivered": "black", "delivered_ft": "0.55", "nominal": "#e0564a",
           "nominal_lt": "#f4a259", "cutsoff": "#26a1ff", "nonan": "#6ce2ad",
           "k0": "#9b59b6", "k0cut": "#e84393", "k0c": "#e84393", "permb": "#c0392b", "qsupport": "#16a085",
           "cvodd": "#2c3e50"}


def load(science, template):
    path = os.path.join(TREE, "lbl_%s_%s.rdb" % (science, template or science))
    if not os.path.exists(path):
        return None
    t, v, e, table = rdb_rows(path)
    names = {c.lower(): c for c in table.colnames}
    berv = (np.asarray(table[names["berv"]], float) if "berv" in names else None)
    if berv is not None and berv.size != t.size:
        berv = None
    return {"path": path, "t": t, "v": v, "e": e, "berv": berv}


def on_common(ref, run):
    key = lambda x: np.round(x, 6)                  # noqa: E731
    shared = np.intersect1d(key(ref["t"]), key(run["t"]))
    out = []
    for r in (ref, run):
        keep = np.isin(key(r["t"]), shared)
        order = np.argsort(r["t"][keep])
        cut = {k: (np.asarray(r[k])[keep][order] if r.get(k) is not None else None)
               for k in ("t", "v", "e", "berv")}
        out.append(cut)
    return out


def nights_of(t):
    """La nuit d'un rjd, au sens du cube: floor(BJD - 0.5), BJD = rjd + 2400000."""
    return np.floor(np.asarray(t, dtype=float) - 0.5).astype(int)


def held_nights(target):
    """Les nuits que le fit de validation croisee a tenues a l'ecart, lues dans
    ses propres coefficients (lignes marquees `rejected`), et non recalculees:
    le rang d'une nuit depend de l'ensemble des nuits, qui n'est pas le meme
    dans le cube et dans la serie LBL."""
    import csv
    path = os.path.join(REPO, "outputs", "_n17_cvodd", target, "1-3",
                        "coefficients.csv")
    if not os.path.exists(path):
        log("%s: no %s, the held-out nights are unknown" % (target, path), "warn")
        return None
    bjd = []
    for row in csv.DictReader(open(path)):
        if str(row.get("rejected", "")).strip() in ("True", "1", "true"):
            value = float(row["bjd"])
            bjd.append(value - 2400000.0 if value > 2400000.0 else value)
    return np.unique(nights_of(bjd))


def intra_night(t, v, e):
    """Dispersion poolee autour de la moyenne ponderee de chaque nuit."""
    nights = nights_of(t)
    dev, dof = [], 0
    for n in np.unique(nights):
        s = nights == n
        if s.sum() < 2:
            continue
        w = 1 / e[s] ** 2
        dev.append(v[s] - np.sum(w * v[s]) / w.sum())
        dof += s.sum() - 1
    if not dev:
        return np.nan, 0
    dev = np.concatenate(dev)
    return float(np.sqrt(np.sum(dev ** 2) / dof)), int(dof)


def nightly_means(t, v, e):
    nights = nights_of(t)
    rows = []
    for n in np.unique(nights):
        s = nights == n
        w = 1 / e[s] ** 2
        rows.append((t[s].mean(), np.sum(w * v[s]) / w.sum(), 1 / np.sqrt(w.sum())))
    return np.array(rows).T


def minus_planets(t, v, e, periods):
    """Residu apres ajustement pondere d'une constante et d'une sinusoide par
    periode (amplitude et phase libres)."""
    cols = [np.ones_like(t)]
    for p in periods:
        cols += [np.sin(2 * np.pi * t / p), np.cos(2 * np.pi * t / p)]
    A = np.array(cols).T
    w = 1 / e
    coef, *_ = np.linalg.lstsq(A * w[:, None], v * w, rcond=None)
    amps = [float(np.hypot(coef[1 + 2 * i], coef[2 + 2 * i]))
            for i in range(len(periods))]
    return v - A @ coef, amps


def paired_gain(ref_v, run_v, draws=2000, seed=11):
    rng = np.random.default_rng(seed)
    n = ref_v.size
    g = [np.std(ref_v[s]) / np.std(run_v[s])
         for s in (rng.integers(0, n, n) for _ in range(draws))]
    return float(np.std(ref_v) / np.std(run_v)), [float(x) for x in
                                                   np.percentile(g, [16, 84])]


def describe(target, key, label, run, ref, seed):
    ref_c, run_c = on_common(ref, run)
    s = velocity_stats(run_c["t"], run_c["v"], run_c["e"])
    intra, dof = intra_night(run_c["t"], run_c["v"], run_c["e"])
    nt = nightly_means(run_c["t"], run_c["v"], run_c["e"])
    out = {"label": label, "path": run["path"], "n": int(run_c["t"].size),
           "rms": s["rms"], "robust": s["robust"],
           "median_error": s["median_error"], "nightly_rms": s["nightly_rms"],
           "intra_night": intra, "intra_dof": dof}
    if target in PLANETS:
        res, amps = minus_planets(nt[0], nt[1], nt[2], PLANETS[target])
        out["nightly_rms_minus_planets"] = float(np.std(res))
        out["planet_amplitudes"] = amps
    if key != "delivered":
        # LA LIVREE SUR LES MEMES POSES, pour que chaque ligne se lise contre
        # sa propre reference et non contre la livree entiere
        r = velocity_stats(ref_c["t"], ref_c["v"], ref_c["e"])
        r_intra, _ = intra_night(ref_c["t"], ref_c["v"], ref_c["e"])
        out["ref"] = {"rms": r["rms"], "robust": r["robust"],
                      "median_error": r["median_error"],
                      "nightly_rms": r["nightly_rms"], "intra_night": r_intra}
        if target in PLANETS:
            rnt = nightly_means(ref_c["t"], ref_c["v"], ref_c["e"])
            out["ref"]["nightly_rms_minus_planets"] = float(np.std(
                minus_planets(rnt[0], rnt[1], rnt[2], PLANETS[target])[0]))
        out["gain"], out["gain_16_84"] = paired_gain(ref_c["v"], run_c["v"])
        added = (run_c["v"] - np.median(run_c["v"])) - (ref_c["v"] - np.median(ref_c["v"]))
        out["added_rms"] = float(np.std(added))
        if run_c["berv"] is not None:
            out["added_vs_berv"] = float(spearmanr(added, run_c["berv"]).statistic)
            vtot = bervbias.total_velocity(run_c["v"], run_c["berv"])
            out["added_vs_vtot"] = float(spearmanr(added, vtot).statistic)
    if run_c["berv"] is not None:
        vtot = bervbias.total_velocity(run_c["v"], run_c["berv"])
        fitted = bervbias.fit(vtot, run_c["v"], run_c["e"], t=run_c["t"], seed=seed)
        out["vtot_bias"] = bervbias.summary(fitted)
        if fitted is not None:
            out["vtot_jitter"] = float(fitted["jitter"][0])
            out["vtot_significance"] = float(fitted["significance"])
    return out, ref_c, run_c


def analyse(target):
    ref = load(target, target)
    if ref is None:
        log("%s: no delivered series yet" % target, "warn")
        return None
    seed = zlib.crc32(target.encode()) & 0xffffffff
    found, curves = {}, {}
    for key, label, sci, tmpl in SERIES:
        sci = sci.format(T=target)
        tmpl = tmpl.format(T=target) if tmpl else None
        run = ref if key == "delivered" else load(sci, tmpl)
        if run is None:
            continue
        if key == "cvodd":
            # seules les nuits TENUES A L'ECART disent quelque chose
            held = held_nights(target)
            if held is None:
                continue
            keep = np.isin(nights_of(run["t"]), held)
            run = {k: (np.asarray(v)[keep] if isinstance(v, np.ndarray) else v)
                   for k, v in run.items()}
        out, ref_c, run_c = describe(target, key, label, run, ref, seed)
        found[key] = out
        curves[key] = (run_c, ref_c)
        line = ("%-34s n %4d | rms %6.2f robust %6.2f err %5.2f | nightly %6.2f"
                " | intra %5.2f" % (label, out["n"], out["rms"], out["robust"],
                                    out["median_error"], out["nightly_rms"],
                                    out["intra_night"]))
        if "nightly_rms_minus_planets" in out:
            line += " | -planets %5.2f" % out["nightly_rms_minus_planets"]
        if "gain" in out:
            line += " | gain %.3f [%.3f, %.3f]" % (out["gain"], *out["gain_16_84"])
        log(line, "value")
        if "ref" in out:
            r = out["ref"]
            line = ("%-34s   delivered on these %d: rms %6.2f robust %6.2f err"
                    " %5.2f | nightly %6.2f | intra %5.2f" % ("", out["n"],
                    r["rms"], r["robust"], r["median_error"], r["nightly_rms"],
                    r["intra_night"]))
            if "nightly_rms_minus_planets" in r:
                line += " | -planets %5.2f" % r["nightly_rms_minus_planets"]
            log(line, "value")
        if "vtot_bias" in out:
            log("%-34s V_tot: %s; jitter %.2f" % ("", out["vtot_bias"],
                                                   out.get("vtot_jitter", np.nan)),
                "value")
    # pour la validation croisee, la reference juste: cutsoff sur les memes nuits
    if "cvodd" in found and "cutsoff" in curves:
        held_n = held_nights(target)
        cut = load("%s_PCA2D_1-3_n17_cutsoff" % target, None)
        cv = load("%s_PCA2D_1-3_n17_cvodd" % target, None)
        a = {k: (np.asarray(v)[np.isin(nights_of(cut["t"]), held_n)]
                 if isinstance(v, np.ndarray) else v) for k, v in cut.items()}
        b = {k: (np.asarray(v)[np.isin(nights_of(cv["t"]), held_n)]
                 if isinstance(v, np.ndarray) else v) for k, v in cv.items()}
        held_nights_count = np.unique(nights_of(b["t"])).size
        a_c, b_c = on_common(a, b)
        found["cv_comparison"] = {
            "held_nights": int(held_nights_count), "n": int(a_c["t"].size),
            "in_sample_rms": float(np.std(a_c["v"])),
            "held_out_rms": float(np.std(b_c["v"])),
            "in_sample_nightly": velocity_stats(a_c["t"], a_c["v"], a_c["e"])["nightly_rms"],
            "held_out_nightly": velocity_stats(b_c["t"], b_c["v"], b_c["e"])["nightly_rms"],
            "in_sample_error": float(np.median(a_c["e"])),
            "held_out_error": float(np.median(b_c["e"]))}
        c = found["cv_comparison"]
        log("cross-validation, %d held-out nights, %d exposures: rms in-sample"
            " %.2f vs held-out %.2f; nightly %.2f vs %.2f; error %.2f vs %.2f"
            % (c["held_nights"], c["n"], c["in_sample_rms"], c["held_out_rms"],
               c["in_sample_nightly"], c["held_out_nightly"],
               c["in_sample_error"], c["held_out_error"]), "value")
    with open(os.path.join(HERE, "metrics_%s.json" % target), "w") as handle:
        json.dump(found, handle, indent=1, default=float)
    plot(target, curves, found)
    return found


def plot(target, curves, found):
    keys = [k for k, *_ in SERIES if k in curves]
    fig, ax = plt.subplots(figsize=(11, 1.1 + 0.9 * len(keys)))
    step = 6 * max(found["delivered"]["robust"], 1.0)
    for i, key in enumerate(keys):
        run = curves[key][0]
        nt = nightly_means(run["t"], run["v"], run["e"])
        off = -i * step
        ax.errorbar(nt[0], nt[1] - np.median(nt[1]) + off, yerr=nt[2], fmt="o",
                    ms=2.5, lw=0.6, color=COLOURS[key], alpha=0.85)
        ax.text(1.005, off, "%s\nnightly rms %.2f" % (found[key]["label"],
                                                      found[key]["nightly_rms"]),
                transform=ax.get_yaxis_transform(), fontsize=7,
                color=COLOURS[key], va="center")
    ax.set_xlabel("rjd (BJD - 2400000)")
    ax.set_ylabel("nightly velocity - median, offset per series (m/s)")
    ax.set_title("%s: every series of the night, nightly means" % target)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = os.path.join(HERE, "compilation_%s.pdf" % target)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    log("wrote %s" % path, "info")

    others = [k for k in keys if k != "delivered"]
    if not others:
        return
    fig, axes = plt.subplots(len(others), 1, figsize=(7.5, 1.9 * len(others)),
                             sharex=True, squeeze=False)
    for ax, key in zip(axes[:, 0], others):
        run, ref = curves[key]
        if run["berv"] is None:
            continue
        added = (run["v"] - np.median(run["v"])) - (ref["v"] - np.median(ref["v"]))
        vtot = bervbias.total_velocity(run["v"], run["berv"])
        ax.plot(vtot, added, ".", ms=2.5, color=COLOURS[key], alpha=0.6)
        ax.axhline(0, color="0.4", lw=0.7)
        ax.set_ylabel("added (m/s)", fontsize=8)
        ax.set_title("%s: added rms %.2f m/s, Spearman vs V_tot %+.2f"
                     % (found[key]["label"], found[key].get("added_rms", np.nan),
                        found[key].get("added_vs_vtot", np.nan)), fontsize=8)
        ax.grid(alpha=0.25)
    axes[-1, 0].set_xlabel("V_tot = vrad - BERV (km/s)")
    fig.tight_layout()
    path = os.path.join(HERE, "injected_%s.pdf" % target)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    log("wrote %s" % path, "info")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target", action="append")
    args = p.parse_args(argv)
    for target in args.target or ["GJ1_N17", "PROXIMA_N17"]:
        log("=== %s" % target, "info")
        analyse(target)


if __name__ == "__main__":
    main()
