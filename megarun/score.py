#!/usr/bin/env python
"""What every run of a sweep was worth, in one table.

    python megarun/score.py <output root> [<output root> ...] \
        --csv megarun/results.csv --status megarun/status.md

The metric of the optimisation is the one the plan names: the robust sigma
(1.4826 MAD) of the LBL velocities after the correction, against the same
star's delivered velocities on the exposures both have. This reads the runs
themselves (pca2d.runs), finds each one's LBL velocities through its own
resolved configuration, and writes:

    results.csv   one row per star per run: the settings that run was given
                  and every number, so a later question can be asked of the
                  file rather than of the runs again
    status.md     the same, sorted by what each run did to the robust sigma,
                  for reading on another machine after a `git pull`

A run still going, or one whose LBL has not finished, is listed with what it
has. Nothing here reads a cube or a spectrum.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pca2d import runs as runs_module                          # noqa: E402
from pca2d import texreport as tr                              # noqa: E402
from pca2d.config import setting_value                         # noqa: E402
from pca2d.lbl import object_names                             # noqa: E402
from pca2d.lblscan import velocity_stats                       # noqa: E402

#: the settings a sweep moves, as they are named in a resolved configuration
SWEPT = (("weight", "correct.weight"), ("n_star", "twoframe.n_star"),
         ("n_earth", "twoframe.n_earth"),
         ("velocity_term", "twoframe.velocity_term"),
         ("high_pass_kms", "highpass.width_kms"), ("dv", "domain.dv"),
         ("shrink", "correct.shrink"), ("mean", "twoframe.mean"),
         ("star_basis", "twoframe.star_basis"),
         ("nightly_stack", "input.nightly_stack"))
FIELDS = (["star", "run", "started", "targets", "joint", "tag"]
          + [name for name, _path in SWEPT]
          + ["n", "nights", "rms_before", "rms_after", "robust_before",
             "robust_after", "nightly_before", "nightly_after",
             "median_error_before", "median_error_after", "gain_robust",
             "folder", "report"])


def stars_of(run):
    """The objects a run fitted: its members when it is joint, else its own."""
    members = sorted(glob.glob(os.path.join(run["folder"],
                                            "cube_config_*.yaml")))
    if members:
        return [os.path.basename(m)[len("cube_config_"):-len(".yaml")]
                for m in members]
    return [run["objects"]] if run["objects"] else []


def scored(run, lbl_dir=None):
    """One row per star of one run, or none when LBL has measured nothing."""
    config = run["config"]
    tree = tr.find_tree(config, run["folder"], lbl_dir)
    rows = []
    for star in stars_of(run):
        before_name, after_name = object_names(config, star, run["tag"])
        before = tr.load(tree, before_name, "delivered") if tree else None
        after = tr.load(tree, after_name, "corrected") if tree else None
        if before is None or after is None:
            continue
        before, after = tr.common(before, after)
        if before["t"].size < 4:
            continue
        b = velocity_stats(before["t"], before["v"], before["e"])
        a = velocity_stats(after["t"], after["v"], after["e"])
        row = {"star": star, "run": run["name"], "started": run["started"],
               "targets": run["objects"], "joint": len(stars_of(run)) > 1,
               "tag": run["tag"], "folder": run["folder"],
               "report": run["report"] or "", "n": b["n"],
               "nights": a["nights"]}
        for name, path in SWEPT:
            row[name] = setting_value(config, path)
        for key, short in (("rms", "rms"), ("robust", "robust"),
                           ("nightly_rms", "nightly"),
                           ("median_error", "median_error")):
            row[short + "_before"] = round(b[key], 4)
            row[short + "_after"] = round(a[key], 4)
        row["gain_robust"] = round(b["robust"] / a["robust"], 4) \
            if a["robust"] > 0 else ""
        rows.append(row)
    return rows


def table(roots, lbl_dir=None):
    """Every run under every root, scored, the best gain first."""
    rows = []
    for root in roots:
        for run in runs_module.listed(root):
            rows.extend(scored(run, lbl_dir))
    rows.sort(key=lambda r: -(r["gain_robust"] or 0))
    return rows


def write_csv(rows, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_status(rows, path, roots, waiting=()):
    """The same table as markdown, for reading after a `git pull`."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    lines = ["# The optimisation megarun, as it stands", "",
             "Written by `megarun/score.py` at %s, from %s."
             % (now, ", ".join("`%s`" % r for r in roots)), "",
             "The metric is the robust sigma (1.4826 MAD) of the LBL"
             " velocities: `gain` is the delivered one divided by the"
             " corrected one, so above 1 is better.", "",
             "| star | run | n* | nE | weight | high pass | vel. term |"
             " robust before | robust after | gain | report |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
             " --- |"]
    for row in rows:
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %.2f | %.2f |"
                     " **%.2f** | %s |"
                     % (row["star"], row["run"], row["n_star"],
                        row["n_earth"], row["weight"], row["high_pass_kms"],
                        row["velocity_term"], row["robust_before"],
                        row["robust_after"], row["gain_robust"] or 0,
                        "yes" if row["report"] else ""))
    if waiting:
        lines += ["", "## Runs with no velocities yet", ""]
        lines += ["- %s (%s), started %s" % (w["objects"], w["name"],
                                             (w["started"] or "")[:19])
                  for w in waiting]
    lines.append("")
    with open(path, "w") as handle:
        handle.write("\n".join(lines))
    return path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("roots", nargs="+", help="output roots to read")
    p.add_argument("--lbl-dir", default=None,
                   help="LBL's tree, when it is not the one the runs name")
    p.add_argument("--csv", default="megarun/results.csv")
    p.add_argument("--status", default="megarun/status.md")
    args = p.parse_args(argv)
    rows = table(args.roots, args.lbl_dir)
    scored_folders = {row["folder"] for row in rows}
    waiting = [run for root in args.roots for run in runs_module.listed(root)
               if run["folder"] not in scored_folders]
    write_csv(rows, args.csv)
    write_status(rows, args.status, args.roots, waiting)
    print("%d row(s) from %d run(s) scored, %d waiting -> %s, %s"
          % (len(rows), len(scored_folders), len(waiting), args.csv,
             args.status))
    return 0


if __name__ == "__main__":
    sys.exit(main())
