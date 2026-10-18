#!/usr/bin/env python
"""Le rapport du week-end du 2026-10-17, en PDF.

    python tools/weekend/report.py

Tableaux: outputs/reports/weekend/{sample,faint_sample,weekend,cuts_paired}.json.
Figures: tools/weekend/report_figures.py (a lancer avant). Texte des
conclusions: outputs/reports/weekend/findings.tex, ecrit a la main. Liste des
runs: done.txt et VERSIONS.txt.
"""

from __future__ import annotations

import json
import os
import subprocess

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
W = os.path.join(REPO, "outputs", "reports", "weekend")
FIG = os.path.join(W, "figures")
NAMES = {"GJ1": "GJ 1", "GL699": "Barnard (GJ 699)", "YZ_CMI": "YZ CMi",
         "GL406": "Wolf 359 (GJ 406)", "GJ643": "GJ 643", "GJ3737": "GJ 3737",
         "GJ3618": "GJ 3618", "GJ3622": "GJ 3622", "TOI5524": "TOI-5524",
         "TOI210": "TOI-210", "GJ3693": "GJ 3693", "LP_938M71": "LP 938-71",
         "TOI2407": "TOI-2407", "TOI756": "TOI-756",
         "SIPS_J1259M4336": "SIPS J1259-4336", "TOI206": "TOI-206",
         "TOI782": "TOI-782", "PROXIMA": "Proxima", "SMETHELLS_20": "Smethells 20",
         "G_99M49": "G 99-49", "GJ682": "GJ 682"}


def tex(text):
    out = []
    for ch in str(text):
        out.append({"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$",
                    "&": r"\&", "#": r"\#", "^": r"\^{}", "_": r"\_", "%": r"\%",
                    "~": r"\~{}"}.get(ch, ch))
    return "".join(out).replace("--", "-{}-").replace("/", "/\\allowbreak{}")


def f(x, d=2):
    return "--" if x is None else "%.*f" % (d, x)


def name(o):
    return NAMES.get(o, o)


def table(caption, cols, head, rows, small=True):
    return ("\\begin{table}[h]\\centering%s\n\\caption{%s}\n\\begin{tabular}{%s}\\toprule\n"
            "%s \\\\\n\\midrule\n%s\n\\bottomrule\\end{tabular}\\end{table}\n"
            % ("\\footnotesize" if small else "", caption, cols, head,
               "\n".join(r + " \\\\" for r in rows)))


def figure(stem, caption, width="\\textwidth"):
    path = os.path.join(FIG, stem + ".pdf")
    if not os.path.exists(path):
        return ""
    return ("\\begin{figure}[h]\\centering\\includegraphics[width=%s]{%s}"
            "\\caption{%s}\\end{figure}\n" % (width, path, caption))


def sample_tables():
    s = json.load(open(os.path.join(W, "sample.json")))
    chosen = {r["object"] for r in s["chosen"]}
    rows = []
    for c in sorted(s["candidates"], key=lambda c: c["H"]):
        rows.append("%s%s & %s & %s & %s & %d & %d" % (
            "\\textbf{" if c["object"] in chosen else "", tex(name(c["object"])) +
            ("}" if c["object"] in chosen else ""), tex(c["sp_type"] or "--"),
            f(c["H"]), f(c.get("teff"), 0), int(c["n_nights"]), int(c["n_files"])))
    out = table("The fifteen NIRPS targets with more than 100 nights in their majority mode"
                " (all HE), counted on rali. The ten drawn are in bold: five equal-width bins"
                " of $H$ between %.2f and %.2f, two drawn per bin, seed %d."
                % (s["edges"][0], s["edges"][-1], s["seed"]),
                "l l r r r r", "Target & SIMBAD type & $H$ & $T_{\\rm eff}$ & nights & exposures",
                rows)
    fs = os.path.join(W, "faint_sample.json")
    if os.path.exists(fs):
        rows = []
        for c in json.load(open(fs)):
            rows.append("%s & %s & %s & %s & %d & %d" % (
                tex(name(c["object"])), tex(c["sp_type"] or "--"), f(c["H"]),
                f(c["extsn060"], 1), int(c["n_nights"]), int(c["n_files"])))
        out += table("The faint set: every NIRPS target with more than 50 nights and a"
                     " pipeline SNR (EXTSN060, median of its LBL rdb) below 40, TOI-210"
                     " excepted, already in the drawn sample.",
                     "l l r r r r", "Target & SIMBAD type & $H$ & EXTSN060 & nights & exposures",
                     rows)
    return out


def results_table():
    t = json.load(open(os.path.join(W, "weekend.json")))
    rows = []
    for r in sorted(t, key=lambda r: (r.get("group", "random") != "random", r["H"])):
        n, c = r.get("nominal"), r.get("cutsoff")
        if not n:
            continue
        ref = n["ref"]
        rows.append("%s & %s & %s & %s & %s & %s & %s & %s & %.2f\\,[%.2f, %.2f] & %s & %s" % (
            tex(name(r["object"])), f(r["H"]), f(r.get("snr"), 0), f(r.get("berv_span"), 0),
            f(ref["rms"]), f(ref["robust"]), f(ref["median_error"]),
            f(n["rms"]), n["gain"], *n["gain_16_84"], "%+.2f" % n["removed_rms"],
            "%.2f" % c["gain"] if c else "--"))
    return table("Every target on one code version: the delivered series and the"
                 " precleaned one (1-3), on the exposures they share, in m/s. ``gain'' is"
                 " the rms ratio delivered / precleaned, paired bootstrap 16--84\\%;"
                 " ``removed'' is $\\pm\\sqrt{|\\mathrm{rms}_{\\rm del}^2-\\mathrm{rms}_{\\rm"
                 " pre}^2|}$, what the correction took out (positive) or added (negative),"
                 " insensitive to the star's own signal; the last column is the gain without"
                 " the residual cuts.",
                 "l r r r r r r r l r r",
                 "Target & $H$ & SNR & BERV & rms & robust & error & rms 1-3 & gain 1-3 &"
                 " removed & no cuts", rows)


def runs():
    lines = []
    for line in open(os.path.join(W, "done.txt")):
        p = line.rstrip("\n").split("|")
        if len(p) < 6:
            continue
        cmd = p[4].split("pca2d.cli")[-1].replace("--fits-dir /Volumes/irrisor/pca2d_preclean", "")
        cmd = cmd.replace("--lbl-dir /Volumes/irrisor/pca2d_preclean/lbl", "")
        if p[0].startswith("W"):
            cmd = "(wait for 150 GB free on irrisor)"
        if p[0].startswith("D_"):
            cmd = "rsync from rali, HE"
        lines.append("%s & %s & %s & %s & \\texttt{%s} \\\\" % (
            tex(p[0]), p[1].replace("code=", ""), p[2][7:15], p[3][7:15], tex(cmd.strip()[:130])))
    return ("\\begin{longtable}{l l l l p{11.5cm}}\\toprule ID & exit & start & end & command"
            " \\\\ \\midrule\n%s\n\\bottomrule\\end{longtable}\n" % "\n".join(lines))


def main():
    versions = open(os.path.join(W, "VERSIONS.txt")).read().split()
    body = [r"""\documentclass[11pt]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{graphicx,booktabs,longtable,xcolor,hyperref,amsmath,float}
\usepackage[T1]{fontenc}\usepackage{lmodern}
\title{Two-frame precleaning on seventeen NIRPS M dwarfs:\\ what the correction does, and when}
\author{pca2d-preclean}
\date{Weekend of 2026-10-17}
\begin{document}\maketitle
"""]
    body.append(r"\noindent\textbf{Code:} pca2d \texttt{%s}, lbl-rapide \texttt{%s}, frozen at"
                r" the start, %s %s. One NIRPS mode (HE), one version, every target under a"
                r" fresh name so that no earlier LBL product is read back." %
                (versions[1][:10], versions[3][:10], versions[5], versions[6]))
    findings = os.path.join(W, "findings.tex")
    body.append("\n\\input{%s}\n" % findings if os.path.exists(findings)
                else "\n\\emph{findings.tex not written yet.}\n")
    body.append("\\clearpage\\section{The sample}\n" + sample_tables())
    body.append("\\clearpage\\section{Results}\n" + results_table())
    body.append(figure("f1_predictor", "The gain against the excess scatter of the delivered"
                       " series (left: rms; middle: robust sigma) and against the SNR (right)."))
    body.append(figure("f2_removed", "What the correction removes in quadrature, against"
                       " the delivered scatter. The dashed line is everything removed."))
    body.append(figure("f3_errors", "The cost in LBL error bars, with and without the"
                       " residual cuts.", "0.7\\textwidth"))
    body.append(figure("f4_cuts", "The residual cuts measured directly: rms with / rms"
                       " without, same exposures, paired bootstrap.", "0.75\\textwidth"))
    body.append(figure("f5_core_tails", "Removed from the core (robust sigma) against removed"
                       " from the rms: the correction trims the tails and widens the core.",
                       "0.65\\textwidth"))
    body.append("\\clearpage" + figure("f6_compilation", "Every target, nightly means:"
                                       " delivered, precleaned 1-3, and without the residual"
                                       " cuts, on the exposures the three share."))
    for target, label in (("TOI210_W", "TOI-210, the largest gain"),
                          ("GJ643_W", "GJ 643, a bright star that gains")):
        seq = os.path.join(REPO, "outputs", target, "1-3", "report", "figures",
                           "03-the-sequence-step-by-step.pdf")
        if os.path.exists(seq):
            body.append("\\begin{figure}[H]\\centering\\includegraphics[width=0.8\\textwidth]"
                        "{%s}\\caption{%s: the sequence figure of its own run report.}"
                        "\\end{figure}\n" % (seq, label))
    body.append("\\clearpage\\section{What ran}\n" + runs())
    body.append("\\end{document}\n")
    out = os.path.join(W, "rapport_weekend.tex")
    open(out, "w").write("\n".join(body))
    for _ in range(2):
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", W, out],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pdf = out.replace(".tex", ".pdf")
    print("wrote", pdf if os.path.exists(pdf) else "(no PDF, see the .log)")


if __name__ == "__main__":
    main()
