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
import copy
import glob
import os
import sys
import time

from .config import VARIANT_META, cache_key, load_config, read_yaml, spectra_dir
from .logger import log
from .progress import human, stage

STAGES = ("cube", "fit", "figures", "correct", "lbl")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="pca2d-preclean", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--object", default=None, metavar="NAME",
                   help="the target, and the name of its directory under"
                        " --data-dir. Case matters: it is also the key looked"
                        " up under `objects:` in the config. Required unless"
                        " --objects names several")
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
    p.add_argument("--objects", default=None, metavar="A,B,C",
                   help="fit these objects TOGETHER against one observer basis,"
                        " each keeping its own star spectrum per order parity"
                        " (pca2d.joint). The atmosphere and the instrument are"
                        " the same for all of them; the stars, their BERV"
                        " coverage and their systemic velocities are not, which"
                        " is what makes the shared basis cleaner. Everything"
                        " lands under <output root>/joint/<A+B+C>/, and each"
                        " object is measured by LBL on its own")
    p.add_argument("--variant", default=None, metavar="NAME",
                   help="variants/NAME.yaml, beside the config, on top of it:"
                        " the nominal plus what the variant changes. Its"
                        " products go to <output root>/_NAME and its LBL object"
                        " is <object>_PCA2D_<M-N>_NAME, unless the file says"
                        " otherwise; with reuse_fit it corrects an existing fit")
    p.add_argument("--dry-run", action="store_true",
                   help="resolve everything, print the plan, touch nothing")
    return p.parse_args(argv)


def load_variant(config_path, name):
    """variants/<name>.yaml, beside the config file: the nominal plus what the
    variant changes. None when no variant is named."""
    if not name:
        return None
    path = os.path.join(os.path.dirname(os.path.abspath(config_path)),
                        "variants", "%s.yaml" % name)
    if not os.path.exists(path):
        log("no variant %s: there is no %s" % (name, path), "error")
        raise SystemExit(2)
    with open(path, "r") as handle:
        return read_yaml(handle) or {}


def name_variant(config, name, variant, out_dir=None):
    """Where a variant's products go and what its LBL object is called, unless
    its file says: beside the nominal's, under the variant's name, the way
    outputs/_star_spl and TOI2120_PCA2D_1-3_star_spl were named. Returns the
    output root the nominal's runs are under, where a reused fit is found."""
    root = config["output"]["directory"]
    if variant is None:
        return root
    if not out_dir and not (variant.get("output") or {}).get("directory"):
        config["output"]["directory"] = os.path.join(root, "_" + name)
    if not (variant.get("lbl") or {}).get("suffix"):
        config["lbl"]["suffix"] = "%s_%s" % (config["lbl"].get("suffix")
                                             or "_PCA2D_{tag}", name)
    config["variant"] = dict({"name": name},
                             **{k: variant[k] for k in VARIANT_META if k in variant})
    return root


def reused_fit(root, base, object_name, tag):
    """The run folder a correction-only variant takes its fit from: the
    nominal's (reuse_fit: nominal) or another variant's."""
    parent = root if base == "nominal" else os.path.join(root, "_" + str(base))
    return os.path.join(parent, object_name, tag)


def check_reused_fit(plan):
    """Refuse a reused fit that is not there, or that was made from another
    cube than the one this variant's settings give."""
    fitdir = plan["fitdir"]
    for name in ("fit.npz", "twoframe_components.fits", "resolved_config.yaml"):
        if not os.path.exists(os.path.join(fitdir, name)):
            log("this variant reuses the fit in %s, which has no %s: run the"
                " variant it names first" % (fitdir, name), "error")
            raise SystemExit(2)
    base = cache_key(load_config(os.path.join(fitdir, "resolved_config.yaml")))
    if base != plan["key"]:
        log("the fit in %s was made from cube %s and this variant's settings"
            " give cube %s: a variant repeats its base's cube settings"
            % (fitdir, base, plan["key"]), "error")
        raise SystemExit(2)


def run_tag(config):
    """<M>-<N>, and `v` when the velocity term is on: the name of a run."""
    return "%d-%d%s" % (config["twoframe"]["n_star"], config["twoframe"]["n_earth"],
                        "v" if config["twoframe"].get("velocity_term", False) else "")


def without_object(config):
    """The config with the object taken out, so two objects built the same way
    give the same cube key."""
    other = copy.deepcopy(config)
    other["input"]["object"] = ""
    return other


def joint_members(args, variant):
    """One entry per object of a joint run: its config, its spectra, its cube.

    Every object is resolved exactly as a solo run resolves it, so each keeps
    its own folder and its own cube; the joint fit adds the one basis they
    share. They must agree on the domain, the step and the high pass, or their
    rows could not sit on one grid.
    """
    members = []
    for name in args.objects:
        cfg = load_config(args.config, object_name=name, data_dir=args.data_dir,
                          out_dir=args.out_dir, instrument=args.instrument,
                          variant=variant)
        directory = spectra_dir(cfg)
        if not os.path.isdir(directory):
            log("no directory %s: %s has no spectra under the input root"
                % (directory, name), "error")
            raise SystemExit(2)
        key = cache_key(cfg)
        members.append({
            "object": name, "config": cfg, "directory": directory,
            "files": sorted(glob.glob(os.path.join(
                directory, cfg["input"].get("pattern", "*t.fits")))),
            "key": key,
            "shared_key": cache_key(without_object(cfg)),
            "cube": os.path.join(cfg["output"]["cache_directory"],
                                 "cube_%s_%s" % (cfg["input"]["format"], key)),
        })
    return members


def joint_plan(args, variant):
    """The plan of a run that fits several objects against one observer basis."""
    from . import joint as _joint

    members = joint_members(args, variant)
    if len({m["shared_key"] for m in members}) > 1:
        log("these objects are not built the same way, so their rows cannot"
            " share a grid: %s" % ", ".join("%s %s" % (m["object"], m["shared_key"])
                                            for m in members), "error")
        raise SystemExit(2)
    config = copy.deepcopy(members[0]["config"])
    # the same command-line overrides a solo run takes
    if args.n_star is not None:
        config["twoframe"]["n_star"] = args.n_star
    if args.n_earth is not None:
        config["twoframe"]["n_earth"] = args.n_earth
    if args.windows:
        config["output"]["windows"] = list(args.windows)
    if args.rebuild_cube:
        config["output"]["use_cache"] = False
        for member in members:
            member["config"]["output"]["use_cache"] = False
    if args.run_lbl:
        config.setdefault("lbl", {})["run"] = True
    name = _joint.joint_name(args.objects)
    config["input"]["object"] = name
    # its own LBL objects, so a joint measurement is never taken for a solo one
    config["lbl"]["suffix"] = "%s_joint" % (config["lbl"].get("suffix")
                                            or "_PCA2D_{tag}")
    tag = run_tag(config)
    outdir = os.path.join(config["output"]["directory"], "joint", name, tag)
    plan = {
        "config": config, "objects": list(args.objects), "members": members,
        "directory": ", ".join(m["directory"] for m in members),
        "files": [f for m in members for f in m["files"]],
        "key": members[0]["shared_key"], "tag": tag, "outdir": outdir,
        "corrdir": os.path.join(outdir, "corrected"),
        "cube": _joint.cube_path(config["output"]["cache_directory"],
                                 config["input"]["format"],
                                 members[0]["shared_key"], args.objects),
    }
    for k, member in enumerate(members):
        member["star_group"] = 2 * k
        member["corrdir"] = os.path.join(plan["corrdir"], member["object"])
    return plan


def resolve(args):
    """The config for this object, and where its pieces will land."""
    name = getattr(args, "variant", None)
    variant = load_variant(args.config, name)
    config = load_config(args.config, object_name=args.object,
                         data_dir=args.data_dir, out_dir=args.out_dir,
                         instrument=args.instrument, variant=variant)
    root = name_variant(config, name, variant, args.out_dir)
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
                       "v" if config["twoframe"].get("velocity_term", False) else "")
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
    if variant and variant.get("reuse_fit"):
        plan["fitdir"] = reused_fit(root, variant["reuse_fit"], args.object, tag)
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
    if plan.get("objects"):
        log("objects     %s, fitted together against ONE observer basis, each"
            " with its own star spectrum per order parity"
            % ", ".join(plan["objects"]), "value")
    if cfg.get("variant"):
        log("variant     %s: variants/%s.yaml on top of the config%s"
            % (cfg["variant"]["name"], cfg["variant"]["name"],
               "; the fit is %s's" % plan["fitdir"] if plan.get("fitdir") else ""),
            "value")


def run_joint_cube(plan, build_main):
    """Each object's own cube, then the one cube their rows share."""
    import yaml

    from . import joint as _joint

    if os.path.isdir(plan["cube"]) and plan["config"]["output"]["use_cache"]:
        log("a joint cube for these objects and this configuration is already"
            " on disk, reusing it: %s" % plan["cube"], "warn")
        return
    for member in plan["members"]:
        if os.path.isdir(member["cube"]) and plan["config"]["output"]["use_cache"]:
            log("  %-10s its own cube is already built: %s"
                % (member["object"], member["cube"]))
            continue
        path = os.path.join(plan["outdir"],
                            "cube_config_%s.yaml" % member["object"])
        with open(path, "w") as handle:
            yaml.safe_dump(member["config"], handle, sort_keys=False,
                           default_flow_style=False)
        log("  %-10s %d spectra of its own" % (member["object"],
                                               len(member["files"])), "info")
        build_main([path] + ([] if plan["config"]["output"]["use_cache"]
                             else ["--no-cache"]))
    _joint.build([m["cube"] for m in plan["members"]],
                 [m["object"] for m in plan["members"]],
                 plan["cube"], plan["written_config"])


def run_cube(plan):
    from .build import main as build_main

    if plan.get("members"):
        return run_joint_cube(plan, build_main)
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
    log("dividing out the observer block, its per-parity mean and %d"
        " components%s, and setting to NaN every sample the fit gave no weight:"
        " exactly panel 3 of the sequence figure"
        % (n_earth, " with %d star components" % n_star if n_star else ""), "info")
    mode = correct_mode(plan)
    if "--refit" in mode:
        log("the fit's rows are nights: every exposure of every fitted night is"
            " corrected with its own coefficients, solved for that file against"
            " the fixed basis", "info")
    # one target, or each object of a joint fit with its own star spectra
    targets = plan.get("members") or [{"object": cfg["input"].get("object"),
                                       "directory": plan["directory"],
                                       "corrdir": plan["corrdir"],
                                       "files": plan["files"], "star_group": 0}]
    for target in targets:
        log("writing %d corrected spectra to %s%s"
            % (len(target["files"]), target["corrdir"],
               " (star spectra %d and %d of the joint fit)"
               % (target["star_group"], target["star_group"] + 1)
               if len(targets) > 1 else ""), "info")
        apply_main([
            "--correct", *mode,
            "--fits", os.path.join(plan.get("fitdir") or plan["outdir"],
                                   "twoframe_components.fits"),
            "--n-star", str(n_star), "--n-earth", str(n_earth),
            "--max-sky-ratio", str(cfg["quality"]["max_sky_ratio"] or 0),
            "--source-dir", target["directory"], "--cube", plan["cube"],
            "--corrected-dir", target["corrdir"], "--overwrite",
            "--star-group", str(target.get("star_group", 0)),
            *clip_args(cfg), *shrink_args(cfg)])


def nightly_stacked(cube):
    """Whether a cube's rows are nights, several exposures coadded, rather than
    one exposure each: its meta.fits counts them in n_exposures."""
    from astropy.table import Table
    path = os.path.join(cube, "meta.fits")
    if not os.path.exists(path):
        return False
    meta = Table.read(path)
    return "n_exposures" in meta.colnames and int(max(meta["n_exposures"])) > 1


def correct_mode(plan):
    """How the correct stage picks the files and their coefficients.

    A fit of single exposures has a coefficient row for every file: --all.
    A fit of nights does not, and the night's row is not any one exposure's:
    the sky of a night is not one sky. Then every t.fits of a fitted night is
    corrected, each with its own amplitudes solved against the fixed basis
    (reconstruct --by-night --refit), from the configuration the cube was
    built with.
    """
    if nightly_stacked(plan["cube"]):
        return ["--by-night", "--refit", "--config", plan["written_config"]]
    return ["--all"]


def clip_args(cfg):
    """The correct stage's residual clip, when correct.nsig_cut asks for one:
    panel 5 beyond that many running robust sigmas over the high pass's window."""
    nsig = (cfg.get("correct") or {}).get("nsig_cut")
    if not nsig:
        return []
    log("and setting to NaN every sample whose residual, panel 5, is beyond %.1f"
        " running robust sigmas" % float(nsig), "info")
    return ["--nsig-cut", str(float(nsig)),
            "--clip-window", str(int(cfg["highpass"]["window"]))]


def shrink_args(cfg):
    """The correct stage's own options, which the sequence figure needs too so
    that its panels 3 and 6 show what the files hold: the shrinkage of the
    observer components by their significance and its two Savitzky-Golay
    variants (pca2d.shrink), both measured in the instrument's resolution
    element, twoframe.resolution, and which samples a corrected file blanks."""
    corr = cfg.get("correct") or {}
    out = []
    mask = str(corr.get("mask") or "exposure")
    if mask != "exposure":
        out += ["--mask", mask]
        log("corrected files blank %s"
            % ("every sample any exposure left unweighted, so each of them"
               " carries the same set of lines" if mask == "common"
               else "nothing: an unweighted sample keeps its delivered flux"),
            "info")
    if corr.get("shrink"):
        out.append("--shrink")
        if corr.get("shrink_smooth"):
            out.append("--shrink-smooth")
    which = corr.get("smooth_components") or []
    if which:
        out += ["--smooth-components", ",".join(str(int(j)) for j in which)]
    tf = cfg.get("twoframe") or {}
    resolution = tf.get("resolution") or tf.get("star_resolution")
    if out and resolution:
        out += ["--resolution", str(float(resolution))]
    if out:
        log("each observer component divided out only where it is significant,"
            " all exposures together%s%s"
            % (", its significance smoothed over a resolution element"
                         if "--shrink-smooth" in out else "",
                         ", components %s smoothed first" % ",".join(map(str, which))
                         if which else ""), "info")
    return out


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
    if plan.get("members"):
        # one LBL object per star, each measured on its own; what they share is
        # the observer basis that corrected them, not their velocities
        for member in plan["members"]:
            theirs = dict(plan)
            theirs["config"] = copy.deepcopy(plan["config"])
            theirs["config"]["input"]["object"] = member["object"]
            theirs["directory"] = member["directory"]
            theirs["files"] = member["files"]
            theirs["corrdir"] = member["corrdir"]
            theirs["fitdir"] = plan.get("fitdir") or plan["outdir"]
            theirs["outdir"] = os.path.join(plan["outdir"], "lbl",
                                            member["object"])
            os.makedirs(theirs["outdir"], exist_ok=True)
            log("LBL for %s, corrected by the joint fit" % member["object"],
                "info")
            prepared = splbl.prepare(theirs)
            if block.get("run", False):
                if not prepared.get("readable", True):
                    log("not running LBL for %s: the profile above cannot read"
                        " its spectra" % member["object"], "error")
                    continue
                splbl.run(prepared["script"])
        return
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
    if args.objects:
        args.objects = [o.strip() for o in str(args.objects).split(",") if o.strip()]
        args.object = args.object or args.objects[0]
        plan = joint_plan(args, load_variant(args.config,
                                             getattr(args, "variant", None)))
    elif not args.object:
        log("--object NAME, or --objects A,B,C to fit several of them against"
            " one observer basis", "error")
        raise SystemExit(2)
    else:
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
    if plan.get("fitdir"):
        # a variant that only changes the correction: the fit, and the cube
        # and figures that go with it, are its base's
        check_reused_fit(plan)
        skipped = [s for s in ("cube", "fit", "figures") if s in wanted]
        wanted = [s for s in wanted if s not in skipped]
        if skipped:
            log("the fit is %s's, so %s not run for this variant"
                % (plan["fitdir"], ", ".join(skipped)), "info")
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
    # and the code that ran it, so the result can be traced and made again
    from .provenance import code_version, stamp as code_stamp
    plan["config"]["provenance"] = dict(code_version())
    log("code        pca2d %s" % code_stamp(plan["config"]["provenance"]), "value")
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
