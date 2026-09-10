"""The one command.

    pca2d-preclean --object TOI-2120

That is the whole interface. The object names a folder under the input root,
the instrument is read from the INSTRUME keyword of the files there, and
`config.yaml` supplies everything else in three layers: what is general, what
belongs to that spectrograph, and what belongs to that target.

Two roots, both in the config and both overridable for one run:

    input.directory    read from <root>/<object>/, never written to
    output.directory   written to <root>/<object>/<M>-<N>/, never read from

Four stages, in order, each skippable and each announcing how long it took:

    cube      the spectra registered onto one log-uniform grid
    fit       the two-frame decomposition, star frame and observer frame
    figures   one multipage PDF, parameters first, then every plot
    correct   the observer block divided out, written as t.fits
    lbl       both sets of spectra handed to LBL, delivered and corrected, as
              two objects in one tree. Prepared by default and run when asked.

Nothing is silent. A run takes tens of minutes and the difference between a
working one and a wedged one has to be visible from across the room.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
import time

from .config import cache_key, load_config, spectra_dir
from .logger import log
from .progress import human, stage

STAGES = ("cube", "fit", "figures", "correct", "lbl")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="pca2d-preclean", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--object", required=True, metavar="NAME",
                   help="the target, and the name of its directory under"
                        " --data-dir. Case matters: it is also the key looked"
                        " up under `objects:` in the config")
    p.add_argument("--config", default="config.yaml",
                   help="default: config.yaml beside the repository root")
    p.add_argument("--data-dir", default=None, metavar="DIR",
                   help="override input.directory from the config. It is the"
                        " input ROOT: spectra are read from <DIR>/<object>/")
    p.add_argument("--out-dir", default=None, metavar="DIR",
                   help="override output.directory from the config. It is the"
                        " output ROOT: everything a run writes lands in"
                        " <DIR>/<object>/<M>-<N>/")
    p.add_argument("--stages", default=",".join(STAGES),
                   help="comma-separated subset of %s, in this order"
                        % ",".join(STAGES))
    p.add_argument("--instrument", default=None,
                   help="override the INSTRUME keyword. Almost never right:"
                        " reading the wrong extensions raises no error, it"
                        " returns different photons")
    p.add_argument("--n-star", type=int, default=None)
    p.add_argument("--n-earth", type=int, default=None)
    p.add_argument("--windows", nargs="+", default=None,
                   help="centre:width in nm; default is the list in the config")
    p.add_argument("--rebuild-cube", action="store_true",
                   help="ignore any cached cube for this configuration")
    p.add_argument("--clean-cache", action="store_true",
                   help="empty the cache first: cubes, figure snippets and"
                        " telluric maps, everything this package rebuilds on"
                        " its own. With --dry-run it only lists what would go")
    p.add_argument("--run-lbl", action="store_true",
                   help="have the lbl stage run LBL, not only prepare it."
                        " Hours. Same as lbl.run: true in the config")
    p.add_argument("--dry-run", action="store_true",
                   help="resolve everything, print the plan, touch nothing")
    return p.parse_args(argv)


def resolve(args):
    """The config for this object, and where its pieces will land."""
    config = load_config(args.config, object_name=args.object,
                         data_dir=args.data_dir, out_dir=args.out_dir,
                         instrument=args.instrument)
    directory = spectra_dir(config)
    if not os.path.isdir(directory):
        log("no directory %s. The object names a folder under the input root,"
            " which is input.directory in %s (%s) unless --data-dir overrides"
            " it." % (directory, args.config, config["input"]["directory"]),
            "error")
        raise SystemExit(2)
    if args.n_star is not None:
        config["twoframe"]["n_star"] = args.n_star
    if args.n_earth is not None:
        config["twoframe"]["n_earth"] = args.n_earth
    if args.windows:
        config["output"]["windows"] = list(args.windows)
    if args.rebuild_cube:
        config["output"]["use_cache"] = False
    if args.run_lbl:
        config.setdefault("lbl", {})["run"] = True

    # `v` when the velocity term is on, so a run with it and a run without it
    # do not write into the same folder, and their LBL objects do not either.
    # Comparing the two is the reason the knob exists at all.
    tag = "%d-%d%s" % (config["twoframe"]["n_star"],
                       config["twoframe"]["n_earth"],
                       "v" if config["twoframe"].get("velocity_term", True) else "")
    # Everything this run writes hangs off one root, the corrected spectra
    # included: they used to land in a subfolder of the input directory, which
    # made the input tree both read and written and meant a shared or
    # read-only archive of spectra could not be used as it stands.
    outdir = os.path.join(config["output"]["directory"], args.object, tag)
    plan = {
        "config": config,
        "directory": directory,
        "files": sorted(glob.glob(os.path.join(
            directory, config["input"].get("pattern", "*t.fits")))),
        "key": cache_key(config),
        "tag": tag,
        "outdir": outdir,
        "corrdir": os.path.join(outdir, "corrected"),
    }
    plan["cube"] = os.path.join(config["output"]["cache_directory"],
                                "cube_%s_%s" % (config["input"]["format"],
                                                plan["key"]))
    return plan


def announce(args, plan):
    """Say what is about to happen, in the terms the config used."""
    cfg = plan["config"]
    log("object      %s" % args.object, "value")
    log("instrument  %s   (read from INSTRUME, not chosen)"
        % cfg["input"].get("instrument", "?"), "value")
    log("spectra     %d files in %s" % (len(plan["files"]), plan["directory"]),
        "value")
    source = ("measured, %.0f%% of the %.2f km/s finest pixel; domain.dv in the"
              " config is not used" % (100 * cfg["domain"]["dv"]
                                       / cfg["domain"]["pixel_dv"],
                                       cfg["domain"]["pixel_dv"])
              if cfg["domain"].get("pixel_dv") else "domain.dv")
    log("domain      %.1f to %.1f nm at %.2f km/s per sample   (%s)"
        % (cfg["domain"]["wave_min"], cfg["domain"]["wave_max"],
           cfg["domain"]["dv"], source), "value")
    log("components  %d in the star's frame, %d in the observer's"
        % (cfg["twoframe"]["n_star"], cfg["twoframe"]["n_earth"]), "value")
    log("cube        %s" % plan["cube"], "value")
    log("outputs     %s" % plan["outdir"], "value")
    log("corrected   %s" % plan["corrdir"], "value")


def run_cube(plan):
    from .build import main as build_main

    if os.path.isdir(plan["cube"]) and plan["config"]["output"]["use_cache"]:
        log("a cube for this exact configuration is already on disk, reusing"
            " it: %s" % plan["cube"], "warn")
        log("--rebuild-cube forces it to be built again", "warn")
        return
    cfg = plan["config"]
    log("building the cube: %d spectra onto %.0f-%.0f nm at %.2f km/s"
        % (len(plan["files"]), cfg["domain"]["wave_min"],
           cfg["domain"]["wave_max"], cfg["domain"]["dv"]), "info")
    log("this reads every file once and takes a few minutes", "info")
    argv = [plan["written_config"]]
    if not cfg["output"]["use_cache"]:
        argv.append("--no-cache")
    build_main(argv)


def run_fit(plan):
    from .twoframe import main as fit_main

    cfg = plan["config"]
    log("fitting %d star and %d observer components by block coordinate"
        " descent, at most %d sweeps"
        % (cfg["twoframe"]["n_star"], cfg["twoframe"]["n_earth"],
           cfg["twoframe"]["iters"]), "info")
    log("this is the long one: it prints an R2 after every sweep, and stops"
        " when chi2 has clearly turned over", "info")
    fit_main(["--config", plan["written_config"], "--cube", plan["cube"],
              "--outdir", plan["outdir"]])


def run_figures(plan):
    from .config import check_windows
    from .figures.bundle import main as bundle_main

    check_windows(plan["config"]["output"]["windows"], plan["config"]["domain"])
    log("drawing %d windows and binding everything into one PDF"
        % len(plan["config"]["output"]["windows"]), "info")
    log("each figure script reads the cube once; expect a few minutes", "info")
    bundle_main(["--config", plan["written_config"], "--cube", plan["cube"],
                 "--outdir", plan["outdir"],
                 "--windows", *[str(w) for w in
                                plan["config"]["output"]["windows"]]])


def run_correct(plan):
    from .reconstruct import main as apply_main

    cfg = plan["config"]
    n_star = cfg["correct"]["n_star"]
    n_earth = cfg["correct"]["n_earth"]
    if n_earth is None:
        n_earth = cfg["twoframe"]["n_earth"]
    log("writing %d corrected spectra to %s"
        % (len(plan["files"]), plan["corrdir"]), "info")
    log("dividing out the observer block, its per-parity mean and %d"
        " components%s, and setting to NaN every sample the fit gave no weight:"
        " exactly panel 3 of the sequence figure"
        % (n_earth, " with %d star components" % n_star if n_star else ""), "info")
    apply_main([
        "--correct", "--all",
        "--fits", os.path.join(plan["outdir"], "twoframe_components.fits"),
        "--n-star", str(n_star), "--n-earth", str(n_earth),
        "--max-sky-ratio", str(cfg["quality"]["max_sky_ratio"] or 0),
        "--source-dir", plan["directory"], "--cube", plan["cube"],
        "--corrected-dir", plan["corrdir"], "--overwrite"])


def run_lbl(plan):
    from . import lbl as splbl

    cfg = plan["config"]
    block = cfg.get("lbl") or {}
    if not block.get("prepare", True) and not block.get("run", False):
        log("lbl.prepare and lbl.run are both off, so nothing to do here",
            "warn")
        return
    ok, detail = splbl.available()
    log("LBL %s" % ("is installed: %s" % detail if ok
                    else "cannot be imported here (%s). The files below are"
                         " still written; `conda env update -f"
                         " environment.yml` puts LBL in this environment."
                         % detail), "value" if ok else "warn")
    # With output.fits_directory set, every folder LBL writes in is a link to
    # the external disk before it starts; lbl/science stays here (storage.py)
    from . import storage as _storage
    for sub in _storage.LBL_FOLDERS:
        _storage.link_dir(plan["config"],
                          os.path.join(block.get("directory") or "lbl", sub))
    prepared = splbl.prepare(plan)
    if block.get("run", False):
        if not prepared.get("readable", True):
            log("not running LBL: the profile above cannot read these spectra,"
                " and it would take a few hundred megabytes of downloads and a"
                " template to find that out again", "error")
            raise SystemExit(2)
        splbl.run(prepared["script"])
        return
    log("LBL is not run by this stage unless asked. Both objects are staged"
        " and everything it needs is written; to run it:", "info")
    log("    python %s" % prepared["script"], "value")
    log("or set lbl.run: true in the config, or pass --run-lbl", "info")


def main(argv=None):
    args = parse_args(argv)
    wanted = [s.strip() for s in args.stages.split(",") if s.strip()]
    unknown = [s for s in wanted if s not in STAGES]
    if unknown:
        log("unknown stage %s; known ones are %s"
            % (", ".join(unknown), ", ".join(STAGES)), "error")
        raise SystemExit(2)

    started = time.time()
    log("pca2d-preclean", "info")
    plan = resolve(args)
    announce(args, plan)
    from . import storage as _storage
    where = _storage.check(plan["config"], dry_run=args.dry_run)
    if where:
        log("kept on     %s   (output.fits_directory; %s links there)"
            % (_storage.external(plan["config"], plan["outdir"]),
               plan["outdir"]), "value")
    if args.clean_cache:
        from . import cache as _cache
        removed, _ = _cache.clean(plan["config"], dry_run=args.dry_run)
        if removed and not args.dry_run and "cube" not in wanted:
            log("the cube stage is not in --stages, so the stages after it have"
                " no cube to read until it runs again", "warn")
    if not plan["files"]:
        log("no spectra matched %s in %s"
            % (plan["config"]["input"].get("pattern"), plan["directory"]),
            "error")
        raise SystemExit(2)
    if args.dry_run:
        log("dry run: stopping here, nothing written", "warn")
        return None

    # The resolved config is written out and every stage is handed THAT file,
    # not the one on the command line: the three layers have been merged and
    # the instrument resolved, so this is what actually ran, and it sits next
    # to the outputs it produced.
    # With output.fits_directory set, the run folder is a link to the external
    # disk, made now, before the first file goes into it (storage.py).
    if where:
        _storage.link_dir(plan["config"], plan["outdir"])
    os.makedirs(plan["outdir"], exist_ok=True)
    plan["written_config"] = os.path.join(plan["outdir"], "resolved_config.yaml")
    import yaml
    with open(plan["written_config"], "w") as fh:
        yaml.safe_dump(plan["config"], fh, sort_keys=False,
                       default_flow_style=False)
    log("resolved configuration written to %s" % plan["written_config"], "info")

    runners = {"cube": run_cube, "fit": run_fit, "figures": run_figures,
               "correct": run_correct, "lbl": run_lbl}
    log("%d stages to run: %s" % (len([s for s in STAGES if s in wanted]),
                                  ", ".join(s for s in STAGES if s in wanted)),
        "info")
    log("on a few hundred exposures the whole thing takes tens of minutes,"
        " most of it in the fit", "info")
    for name in STAGES:
        if name not in wanted:
            log("stage %s: skipped" % name, "warn")
            continue
        with stage("stage %s" % name):
            runners[name](plan)

    log("done in %s" % human(time.time() - started), "info")
    bundle = os.path.join(plan["outdir"], "%s_%s.pdf" % (args.object, plan["tag"]))
    if os.path.exists(bundle):
        log("everything this run produced: %s" % bundle, "value")
    return None


if __name__ == "__main__":
    sys.exit(main())
