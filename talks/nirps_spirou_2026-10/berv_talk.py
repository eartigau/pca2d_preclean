"""La figure BERV de TOI 2120, version diapo: titres courts, ratio 16:9.

Le rapport a la sienne (pca2d.texreport.figure_berv) avec des titres longs
qui portent les chiffres; ici les chiffres sont sur la diapo suivante, donc
on garde la figure et on raccourcit les titres. On passe par les memes
fonctions pour que la figure de la presentation et celle du rapport ne
puissent pas diverger.
"""
import os, sys, zlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

sys.path.insert(0, "/Users/eartigau/GitHubProjects/pca2d_preclean")
os.chdir("/Users/eartigau/GitHubProjects/pca2d_preclean")

from pca2d import texreport as tr
from pca2d.lbl import object_names

RUN = "outputs/TOI2120/1-3"
config = yaml.safe_load(open(os.path.join(RUN, "resolved_config.yaml")))
tag = os.path.basename(os.path.normpath(RUN))
tree = tr.find_tree(config, RUN, None)
star = tr.stars_of(RUN, config)[0]
bn, an = object_names(config, star, tag)
before = tr.load(tree, bn, "delivered")
after = tr.load(tree, an, "corrected")
before, after = tr.common(before, after)
numbers = tr.star_numbers(before, after, seed=zlib.crc32(star.encode()) & 0xffffffff)
fits = {"before": numbers["before"].get("bias"), "after": numbers["after"].get("bias")}

vb, va = tr.vtot_of(before), tr.vtot_of(after)
grid = np.linspace(min(vb.min(), va.min()), max(vb.max(), va.max()), 300)

fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
resid = []
for ax, run, x, key, colour, title in zip(
        axes, (before, after), (vb, va), ("before", "after"),
        (tr.BEFORE, tr.AFTER),
        ("as delivered", "after pca2d precleaning")):
    f = fits[key]
    y = tr.without_dc(run, f)
    resid.append(y)
    tr._close_band(ax)
    tr._points(ax, x, y, run["e"], colour, "exposures", alpha=0.35)
    c, m, e = tr.binned(x, y)
    if c.size:
        ax.errorbar(c, m, yerr=e, fmt="s", ls="none", ms=5, mfc="white",
                    mec=tr.INK, mew=0.9, ecolor=tr.INK, elinewidth=0.8,
                    capsize=0, zorder=5, label="2 km/s bins, median")
    tr._bias_band(ax, f, grid, colour, label=r"fitted bias, 1$\sigma$")
    ax.axhline(0, color=tr.MUTED, lw=0.7)
    ax.set_title(title, fontsize=11, color=tr.INK)
    ax.set_xlabel(tr.VTOT_LABEL, fontsize=10, color=tr.INK)
    tr._style(ax)

axes[0].set_ylim(*tr._limits(*resid))
axes[0].set_ylabel("velocity - fitted DC level (m/s)", fontsize=10, color=tr.INK)
axes[0].legend(fontsize=8, frameon=False, loc="upper right")

# L'ANNOTATION QUI PORTE LE MESSAGE, a la place du titre long
axes[0].annotate(r"$144.8^{+9.2}_{-8.7}$ m/s" "\n" r"$15.8\,\sigma$",
                 xy=(0.5, 120), xytext=(11, 118), fontsize=11,
                 color="#b03a2e", fontweight="bold", ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="#b03a2e", lw=1.3))
axes[1].annotate("none detected\n" r"$1.4\,\sigma$", xy=(14, 95),
                 fontsize=11, color="#147846", fontweight="bold",
                 ha="left", va="center")

fig.tight_layout()
out = "/private/tmp/claude-502/-Users-eartigau-GitHubProjects-pca2d-preclean/2a66dacd-d683-4457-838a-971de3c44460/scratchpad/talk/figs/berv_talk.pdf"
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
