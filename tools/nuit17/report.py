#!/usr/bin/env python
"""The morning report of the night of 2026-10-17, as a LaTeX PDF.

    python outputs/reports/nuit17/report.py

Tables come from metrics_<target>.json (analyse.py), the run list from
done.txt and VERSIONS.txt, the prose from findings.tex, written by hand once
the numbers are in. Figures: the compilation and the added-velocity plots of
analyse.py, and the sequence figure of each nominal run's own report.
"""

from __future__ import annotations

import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/Users/eartigau/GitHubProjects/pca2d_preclean"
TARGETS = [("GJ1_N17", "GJ 1"), ("PROXIMA_N17", "Proxima")]
ORDER = ["delivered", "delivered_ft", "nominal", "cutsoff", "nonan", "k0",
         "k0cut", "k0c", "permb", "qsupport", "cvodd"]


def fmt(x, d=2):
    return "--" if x is None else ("%.*f" % (d, x))


def table(target, name):
    path = os.path.join(HERE, "metrics_%s.json" % target)
    if not os.path.exists(path):
        return "\\emph{No measurements for %s.}\n" % name
    m = json.load(open(path))
    planets = any("nightly_rms_minus_planets" in v for v in m.values()
                  if isinstance(v, dict))
    cols = "l r r r r r r" + (" r" if planets else "") + " l"
    head = ("Series & $N$ & rms & robust & error & nightly & intra"
            + (" & $-$planets" if planets else "") + " & gain \\\\")
    rows = []
    for key in ORDER:
        if key not in m:
            continue
        v = m[key]
        gain = ("%.3f [%.3f, %.3f]" % (v["gain"], *v["gain_16_84"])
                if "gain" in v else "")
        row = ("%s & %d & %s & %s & %s & %s & %s" % (
            v["label"].replace("_", "\\_"), v["n"], fmt(v["rms"]),
            fmt(v["robust"]), fmt(v["median_error"]), fmt(v["nightly_rms"]),
            fmt(v["intra_night"])))
        if planets:
            row += " & %s" % fmt(v.get("nightly_rms_minus_planets"))
        rows.append(row + " & %s \\\\" % gain)
        if "ref" in v:
            r = v["ref"]
            ref = ("\\quad\\footnotesize delivered, same exposures & & %s & %s"
                   " & %s & %s & %s" % (fmt(r["rms"]), fmt(r["robust"]),
                                        fmt(r["median_error"]),
                                        fmt(r["nightly_rms"]),
                                        fmt(r["intra_night"])))
            if planets:
                ref += " & %s" % fmt(r.get("nightly_rms_minus_planets"))
            rows.append("\\color{gray}" + ref + " & \\\\")
    cv = m.get("cv_comparison")
    extra = ""
    if cv:
        extra = ("\n\\noindent Cross-validation on the %d held-out nights (%d"
                 " exposures): rms %.2f in-sample against %.2f held-out; nightly"
                 " %.2f against %.2f; median error %.2f against %.2f.\n"
                 % (cv["held_nights"], cv["n"], cv["in_sample_rms"],
                    cv["held_out_rms"], cv["in_sample_nightly"],
                    cv["held_out_nightly"], cv["in_sample_error"],
                    cv["held_out_error"]))
    return ("\\begin{table}[h]\\centering\\footnotesize\n"
            "\\caption{%s: every series of the night, velocities in m/s."
            " Each corrected series is followed, in grey, by the delivered"
            " series on the same exposures. ``intra'' is the pooled scatter"
            " about each night's weighted mean; the gain is the rms ratio to"
            " the delivered series, paired bootstrap 16--84\\%%.}\n"
            "\\begin{tabular}{%s}\\toprule\n%s\n\\midrule\n%s\n\\bottomrule"
            "\\end{tabular}\\end{table}\n%s" % (name, cols, head,
                                              "\n".join(rows), extra))


def tex(text):
    """Du texte brut en LaTeX: caracteres speciaux echappes, et '--' en deux
    tirets que la police ne fusionne pas."""
    out = []
    for ch in text:
        out.append({"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}",
                    "$": r"\$", "&": r"\&", "#": r"\#", "^": r"\^{}",
                    "_": r"\_", "%": r"\%", "~": r"\~{}"}.get(ch, ch))
    return "".join(out).replace("--", "-{}-").replace("/", "/\\allowbreak{}")


def figure(path, caption, width="\\textwidth"):
    if not os.path.exists(path):
        return ""
    return ("\\begin{figure}[h]\\centering\\includegraphics[width=%s]{%s}"
            "\\caption{%s}\\end{figure}\n" % (width, path, caption))


def runs():
    path = os.path.join(HERE, "done.txt")
    lines = []
    for line in open(path):
        parts = line.rstrip("\n").split("|")
        if len(parts) < 6:
            continue
        ident, status, start, end, cmd, commit = parts[:6]
        short = cmd.split("pca2d.cli")[-1].split("lbl_extra.py")[-1]
        short = short.replace("--fits-dir /Volumes/irrisor/pca2d_preclean", "")
        if ident.startswith("W"):
            short = "(wait until the internal disk has 18 GB free)"
        lines.append("%s & %s & %s & %s & \\texttt{%s} \\\\"
                     % (ident, status.replace("code=", ""), start[7:15],
                        end[7:15], tex(short.strip()[:160])))
    return ("\\begin{longtable}{l l l l p{11cm}}\\toprule ID & exit & start"
            " & end & command \\\\ \\midrule\n%s\n\\bottomrule\\end{longtable}\n"
            % "\n".join(lines))


def main():
    versions = open(os.path.join(HERE, "VERSIONS.txt")).read().strip()
    findings = os.path.join(HERE, "findings.tex")
    body = [r"""\documentclass[11pt]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{graphicx,booktabs,longtable,xcolor,hyperref,amsmath}
\usepackage[T1]{fontenc}\usepackage{lmodern}
\title{Proxima and GJ~1 on the current code:\\ the night of 2026-10-17}
\author{pca2d-preclean}
\date{\today}
\begin{document}\maketitle
"""]
    parts = versions.split()
    body.append(r"\noindent\textbf{Code:} pca2d \texttt{%s}, lbl-rapide \texttt{%s},"
                r" frozen at the start of the night, %s %s."
                % (parts[1][:10], parts[3][:10], parts[5], parts[6]))
    body.append("\n\n\\input{%s}\n" % findings if os.path.exists(findings)
                else "\n\\emph{findings.tex not written yet.}\n")
    for target, name in TARGETS:
        body.append("\\clearpage\\section{%s}\n" % name)
        body.append(table(target, name))
        body.append(figure(os.path.join(HERE, "compilation_%s.pdf" % target),
                           "%s: nightly means of every series, the delivered"
                           " one first, offset for clarity." % name))
        body.append(figure(os.path.join(HERE, "injected_%s.pdf" % target),
                           "%s: what each correction added to the velocities"
                           " (series minus delivered), against $V_{\\rm tot}$."
                           % name, "0.8\\textwidth"))
        seq = os.path.join(REPO, "outputs", target, "1-3", "report", "figures",
                           "03-the-sequence-step-by-step.pdf")
        body.append(figure(seq, "%s, nominal run: the sequence figure of its"
                           " own report." % name, "0.75\\textwidth"))
    body.append("\\clearpage\\section{What ran}\n" + runs())
    body.append("\\end{document}\n")
    tex = os.path.join(HERE, "rapport_nuit17.tex")
    with open(tex, "w") as handle:
        handle.write("\n".join(body))
    for _ in range(2):
        subprocess.run(["pdflatex", "-interaction=nonstopmode",
                        "-output-directory", HERE, tex],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pdf = tex.replace(".tex", ".pdf")
    print("wrote", pdf if os.path.exists(pdf) else "(no PDF, see the .log)")


if __name__ == "__main__":
    main()
