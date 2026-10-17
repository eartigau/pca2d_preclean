"""Le gain en RMS contre le SNR par pose, pour la presentation NIRPS/SPIRou.

Les chiffres viennent des rdb du depot, apparies par epoque; le script les
remesure a chaque execution plutot que de les recopier, de sorte qu'une
figure ne puisse pas vieillir sans qu'on le voie.
"""
import os, csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, "/Users/eartigau/GitHubProjects/pca2d_preclean")
os.chdir("/Users/eartigau/GitHubProjects/pca2d_preclean")
from pca2d.lblscan import read_rdb

R = "lbl/lblrdb"
CASES = [
    ("TOI-2120",  "SPIRou", "outputs/TOI2120/1-3", "outputs/lbl/lblrdb/lbl2_TOI2120_TOI2120.rdb",
     "outputs/lbl/lblrdb/lbl2_TOI2120_PCA2D_1-3_TOI2120_PCA2D_1-3.rdb"),
    ("TOI-1452",  "NIRPS", "outputs/TOI1452/1-3", R + "/lbl2_TOI1452_TOI1452.rdb",
     R + "/lbl2_TOI1452_PCA2D_1-3_TOI1452_PCA2D_1-3.rdb"),
    ("GJ 328",    "NIRPS", "outputs/GJ328/1-3",   R + "/lbl2_GJ328_GJ328.rdb",
     R + "/lbl2_GJ328_PCA2D_1-3_GJ328_PCA2D_1-3.rdb"),
    ("GJ 1",      "NIRPS", "outputs/GJ1/1-3",     R + "/lbl2_GJ1_GJ1.rdb",
     R + "/lbl2_GJ1_PCA2D_1-3_GJ1_PCA2D_1-3.rdb"),
    ("Proxima",   "NIRPS", "outputs/PROXIMA/1-3", R + "/lbl2_PROXIMA_PROXIMA.rdb",
     R + "/lbl2_PROXIMA_PCA2D_1-3_PROXIMA_PCA2D_1-3.rdb"),
]

rng = np.random.default_rng(11)
rows = []
for name, inst, run, before, after in CASES:
    meta = [r for r in csv.DictReader(open(os.path.join(run, "coefficients.csv")))
            if r["rejected"] in ("", "False", "0")]
    snr = float(np.median([float(r["snr_band"]) for r in meta]))
    A, B = read_rdb(before), read_rdb(after)
    ka = {round(float(x), 3): i for i, x in enumerate(A[0])}
    kb = {round(float(x), 3): i for i, x in enumerate(B[0])}
    common = sorted(set(ka) & set(kb))
    va = A[1][[ka[k] for k in common]]
    vb = B[1][[kb[k] for k in common]]
    gain = np.std(va) / np.std(vb)
    # BARRE D'ERREUR: bootstrap apparie sur les epoques, 2000 tirages
    draws = [np.std(va[s]) / np.std(vb[s])
             for s in (rng.integers(0, len(va), len(va)) for _ in range(2000))]
    lo, hi = np.percentile(draws, [16, 84])
    rows.append(dict(name=name, inst=inst, snr=snr, gain=gain, lo=lo, hi=hi,
                     n=len(common), rms_in=np.std(va), rms_out=np.std(vb)))
    print("%-10s %-7s SNR %5.0f  n %4d  %7.2f -> %7.2f  gain %.2f [%.2f, %.2f]"
          % (name, inst, snr, len(common), np.std(va), np.std(vb), gain, lo, hi))

snr = np.array([r["snr"] for r in rows])
gain = np.array([r["gain"] for r in rows])
rho, p = spearmanr(snr, gain)
slope, icept = np.polyfit(np.log10(snr), np.log10(gain), 1)
print("Spearman rho = %+.3f (p = %.3f) | gain ~ SNR^%.2f" % (rho, p, slope))

fig, ax = plt.subplots(figsize=(7.6, 4.8))
xs = np.logspace(np.log10(25), np.log10(360), 100)
ax.plot(xs, 10 ** (icept + slope * np.log10(xs)), lw=1.2, color="0.55", ls="--",
        zorder=1, label=r"power law, gain $\propto$ SNR$^{%.2f}$" % slope)
ax.axhline(1.0, lw=1.0, color="0.3", zorder=1)
ax.fill_between([25, 360], 1.0, 0.72, color="#e0564a", alpha=0.07, zorder=0)
ax.text(28, 0.755, "precleaning makes it worse", fontsize=8, color="#b03a2e",
        ha="left", va="center")

for r in rows:
    colour = "#26a1ff" if r["inst"] == "NIRPS" else "#e0564a"
    marker = "o" if r["inst"] == "NIRPS" else "s"
    ax.errorbar(r["snr"], r["gain"],
                yerr=[[r["gain"] - r["lo"]], [r["hi"] - r["gain"]]],
                fmt=marker, ms=8, color=colour, capsize=3, lw=1.2, zorder=3)
    dy = 1.055 if r["name"] != "Proxima" else 1.055
    ax.annotate("%s\n%d epochs" % (r["name"], r["n"]),
                (r["snr"], r["gain"] * dy), fontsize=8, ha="center",
                va="bottom" if dy > 1 else "top", color="0.25")

ax.plot([], [], "s", color="#e0564a", label="SPIRou")
ax.plot([], [], "o", color="#26a1ff", label="NIRPS")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(25, 360); ax.set_ylim(0.72, 4.0)
ax.set_xticks([30, 50, 100, 200, 300])
ax.set_yticks([0.8, 1.0, 1.5, 2.0, 3.0])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel("median SNR per exposure (band used by the fit)")
ax.set_ylabel(r"gain in RV rms,  $\sigma_{\rm delivered} / \sigma_{\rm precleaned}$")
ax.set_title("Where two-frame precleaning pays: faint targets\n"
             r"Spearman $\rho = %+.2f$ over 5 campaigns ($p = %.3f$)" % (rho, p),
             fontsize=10)
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.25, which="both")
fig.tight_layout()
out = "/private/tmp/claude-502/-Users-eartigau-GitHubProjects-pca2d-preclean/2a66dacd-d683-4457-838a-971de3c44460/scratchpad/talk/figs/gain_vs_snr"
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".pdf")
