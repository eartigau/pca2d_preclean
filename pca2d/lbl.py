"""Hand both sets of spectra to LBL, the uncorrected and the corrected.

The point of correcting a spectrum is the velocity that comes out of it, and
the only honest way to know whether the correction helped is to measure both.
So this stage never runs LBL on the corrected files alone. It sets up two
objects side by side in one LBL tree, from the same instrument profile and with
each building its own template:

    lbl/science/TOI-2120/          symlinks to the spectra as delivered
    lbl/science/TOI-2120_PCA2D/    symlinks to what this package wrote

and LBL's own products then sit next to each other under the same names, so the
comparison that matters is two columns of one table rather than two runs
someone has to remember to line up.

Symlinks, not copies: a few hundred t.fits are a few tens of gigabytes, they
already exist twice (delivered and corrected), and a third copy would buy
nothing. Nothing here writes into the input tree.

What it leaves beside the run's other outputs:

    lbl_config.yaml   LBL's configuration, in LBL's keys and LBL's spelling,
                      which `lbl_compute --config lbl_config.yaml` reads
    run_lbl.py        the wrap script, the runparams dict LBL users edit, with
                      both objects in it and ready to be run by hand

Running LBL takes hours, so the stage prepares by default and runs only when
asked: `lbl.run: true` in the config, or --run-lbl on the command line.
"""

from __future__ import annotations

import os
import subprocess
import sys

import yaml

from .logger import log

#: LBL's name for a spectrograph is not always its INSTRUME, and the pairing
#: with a data source is what selects the class that reads the files: NIRPS is
#: NIRPS_HA or NIRPS_HE depending on the mode it was observed in, which LBL
#: expects the caller to know rather than reading it from the header. So the
#: pairing lives in config.yaml, per instrument, and this is only the fallback
#: for a spectrograph whose block does not say.
FALLBACK = {
    "SPIROU": ("SPIROU", "APERO"),
    "NIRPS": ("NIRPS_HE", "APERO"),
}

STEPS = ("telluclean", "template", "mask", "compute", "compile")


def available():
    """(True, version) if LBL can be imported, (False, why) if it cannot."""
    try:
        import lbl
    except Exception as exc:                                  # noqa: BLE001
        return False, str(exc)
    return True, getattr(lbl, "__version__", "unknown")


def profile(config: dict) -> tuple:
    """(instrument, data_source) as LBL spells them, and where they came from.

    Read from the instrument's own block in config.yaml, because that is where
    everything else about a spectrograph is written and because choosing this
    wrongly is not an error: LBL would read the files with another instrument's
    class and return velocities.
    """
    block = config.get("lbl") or {}
    instrume = str(config["input"].get("instrument") or "").upper()
    name = block.get("instrument")
    source = block.get("data_source")
    if name and source:
        return str(name), str(source), "config.yaml"
    fallback = FALLBACK.get(instrume)
    if not fallback:
        raise SystemExit(
            "no LBL profile for INSTRUME %r. Add `lbl: {instrument: ..., "
            "data_source: ...}` to its block in config.yaml: LBL's name for a"
            " spectrograph is not always the one in the header, and NIRPS is"
            " NIRPS_HA or NIRPS_HE by the mode it was observed in."
            % (instrume or "unset"))
    return str(name or fallback[0]), str(source or fallback[1]), "the fallback table"


def object_names(config: dict, object_name: str, tag: str) -> tuple:
    """The names the two runs go under in the LBL tree, before and after."""
    suffix = str((config.get("lbl") or {}).get("suffix") or "_PCA2D_{tag}")
    return object_name, object_name + suffix.format(tag=tag)


def link_spectra(files, target: str, mode: str = "symlink") -> tuple:
    """Put `files` in `target` as symlinks (or copies).

    Returns (linked, kept, strangers): how many were made, how many were
    already right, and which .fits in the folder are none of ours.

    Idempotent, and never destructive: a name already in the folder is left
    exactly as it is, whether it is one of our symlinks from a previous run or
    a file someone else put there. LBL globs the folder, so a name that is
    already right needs nothing done to it.

    The strangers matter because LBL takes the whole folder. Anything in there
    that this run did not put there is going to be measured with everything
    else, and no error will be raised about it.
    """
    import shutil

    os.makedirs(target, exist_ok=True)
    ours = {os.path.basename(path) for path in files}
    linked = kept = 0
    for path in files:
        destination = os.path.join(target, os.path.basename(path))
        if os.path.lexists(destination):
            kept += 1
            continue
        if mode == "copy":
            shutil.copy2(path, destination)
        else:
            os.symlink(os.path.abspath(path), destination)
        linked += 1
    strangers = sorted(name for name in os.listdir(target)
                       if name.endswith(".fits") and name not in ours)
    return linked, kept, strangers


def config_document(config: dict, data_dir: str, instrument: str,
                    data_source: str) -> dict:
    """The default LBL configuration, as LBL's own keys.

    Only keys LBL knows go in: it validates every one of them against its
    parameter table and refuses a file with a name it does not recognise, which
    is a good reason to write this from the code rather than by hand.
    """
    block = config.get("lbl") or {}
    document = {
        "INSTRUMENT": instrument,
        "DATA_SOURCE": data_source,
        "DATA_DIR": os.path.abspath(data_dir),
        "DATA_TYPE": "SCIENCE",
        # LBL's per-instrument default is already '*.fits'. It is written out
        # because ours must match BOTH sets, and the corrected files are named
        # <stem>t_<M>-<N>.fits: a '*t.fits' glob would find the delivered
        # spectra and silently none of the corrected ones.
        "INPUT_FILE": str(block.get("input_file") or "*.fits"),
        # re-running must not redo the hours that are already on disk
        "SKIP_DONE": True,
    }
    if block.get("teff") is not None:
        document["OBJECT_TEFF"] = block["teff"]
    return document


def runparams(config: dict, data_dir: str, instrument: str, data_source: str,
              objects: list, config_file: str) -> dict:
    """The dict LBL's wrapper takes, with both objects in it.

    Every object is its own comparison, which is to say each builds its own
    template. Sharing one template between the corrected and uncorrected
    spectra would be measuring two things with one ruler that fits neither:
    what the correction removes is exactly what would differ between them.
    `lbl.template` overrides that for anyone who wants it anyway.
    """
    block = config.get("lbl") or {}
    steps = [str(s).strip().lower() for s in (block.get("steps") or [])]
    unknown = [s for s in steps if s not in STEPS]
    if unknown:
        raise SystemExit("unknown lbl.steps %s; known ones are %s"
                         % (", ".join(unknown), ", ".join(STEPS)))
    template = block.get("template")
    params = {
        "INSTRUMENT": instrument,
        "DATA_SOURCE": data_source,
        "DATA_DIR": os.path.abspath(data_dir),
        "CONFIG_FILE": os.path.abspath(config_file),
        "DATA_TYPES": ["SCIENCE"] * len(objects),
        "OBJECT_SCIENCE": list(objects),
        "OBJECT_COMPARISON": [template or name for name in objects],
        "OBJECT_TEFF": [block.get("teff")] * len(objects),
    }
    for step in STEPS:
        params["RUN_LBL_%s" % step.upper()] = step in steps
    # SKIP_LBL_<step> is LBL's "leave what is already done alone", which is
    # what a second run of this stage wants: the template and the mask are
    # hours, and neither depends on anything this package changed.
    for step in ("template", "mask", "compute", "compile"):
        params["SKIP_LBL_%s" % step.upper()] = True
    return params


RUNNER = '''#!/usr/bin/env python
"""Run LBL on %(n)d objects: %(objects)s.

Written by pca2d-preclean for %(object)s (%(tag)s), and left here to be read,
edited and re-run by hand: it is an ordinary LBL wrap script, and the dict
below is the one LBL users already know. What it measures is the same target
twice, as delivered and as this package corrected it, so the two velocity
series can be put side by side.

    python %(script)s

The settings that are not about which object is which live in %(config)s.
"""

from lbl.recipes import lbl_wrap

rparams = %(rparams)s

if __name__ == "__main__":
    lbl_wrap.main(rparams)
'''


def write_runner(path: str, params: dict, object_name: str, tag: str,
                 config_file: str) -> str:
    """The wrap script, with the runparams spelled out one key per line."""
    # UPPER CASE, and not a style choice: lbl_wrap reads runparams['INSTRUMENT']
    # and every other key by that exact spelling, so a lower-case dict fails on
    # the first line of its own checking with "Must define key INSTRUMENT".
    lines = ["dict("]
    for key, value in params.items():
        lines.append("    %s=%r," % (key.upper(), value))
    lines.append(")")
    body = RUNNER % {
        "n": len(params["OBJECT_SCIENCE"]),
        "objects": ", ".join(params["OBJECT_SCIENCE"]),
        "object": object_name,
        "tag": tag,
        "script": os.path.abspath(path),
        "config": os.path.basename(config_file),
        "rparams": "\n".join(lines),
    }
    with open(path, "w") as handle:
        handle.write(body)
    os.chmod(path, 0o755)
    return path


def prepare(plan) -> dict:
    """Stage both objects, write LBL's config and the script that runs it."""
    config = plan["config"]
    block = config.get("lbl") or {}
    object_name = config["input"]["object"]
    instrument, data_source, where = profile(config)
    data_dir = block.get("directory") or "lbl"
    before, after = object_names(config, object_name, plan["tag"])

    log("LBL profile: instrument %s, data source %s, from %s"
        % (instrument, data_source, where), "value")

    wanted = []
    if block.get("before", True):
        wanted.append((before, plan["directory"], plan["files"]))
    if block.get("after", True):
        corrected = sorted(
            os.path.join(plan["corrdir"], name)
            for name in (os.listdir(plan["corrdir"])
                         if os.path.isdir(plan["corrdir"]) else [])
            if name.endswith(".fits"))
        if not corrected:
            log("no corrected spectra in %s yet, so LBL is set up on the"
                " delivered ones only. Run the correct stage and this stage"
                " again to get the other half of the comparison."
                % plan["corrdir"], "warn")
        else:
            wanted.append((after, plan["corrdir"], corrected))

    objects = []
    mode = str(block.get("link") or "symlink")
    for name, source, files in wanted:
        target = os.path.join(data_dir, "science", name)
        linked, kept, strangers = link_spectra(files, target, mode)
        log("%-24s %d spectra (%d new %s, %d already there) from %s"
            % (name, linked + kept, linked, mode + "s", kept, source), "value")
        if strangers:
            log("  %d other .fits are in %s and LBL takes the whole folder, so"
                " they will be measured with the rest: %s%s"
                % (len(strangers), target, ", ".join(strangers[:3]),
                   ", ..." if len(strangers) > 3 else ""), "warn")
        objects.append(name)

    document = config_document(config, data_dir, instrument, data_source)
    config_file = os.path.join(plan["outdir"], "lbl_config.yaml")
    with open(config_file, "w") as handle:
        handle.write("# LBL's own configuration, written by pca2d-preclean.\n"
                     "# Every key here is one LBL knows; it refuses any other.\n")
        yaml.safe_dump(document, handle, sort_keys=False,
                       default_flow_style=False)

    params = runparams(config, data_dir, instrument, data_source, objects,
                       config_file)
    script = write_runner(os.path.join(plan["outdir"], "run_lbl.py"), params,
                          object_name, plan["tag"], config_file)

    if params["RUN_LBL_MASK"] and block.get("teff") is None:
        log("lbl.teff is not set. LBL picks the stellar model its mask comes"
            " from by effective temperature, so put the target's Teff in the"
            " config before running the mask step.", "warn")

    log("LBL config   %s" % config_file, "value")
    log("LBL script   %s" % script, "value")
    return {"config_file": config_file, "script": script, "objects": objects,
            "data_dir": data_dir}


def run(script: str) -> None:
    """Run the script that was just written, in this interpreter's env."""
    ok, detail = available()
    if not ok:
        raise SystemExit(
            "lbl.run is on but LBL cannot be imported here: %s. It is in"
            " environment.yml; `conda env update -f environment.yml` puts it"
            " in this environment." % detail)
    log("running LBL %s. This is hours, and LBL prints its own progress."
        % detail, "info")
    result = subprocess.run([sys.executable, script])
    if result.returncode != 0:
        raise SystemExit("LBL exited %d; the script that ran is %s"
                         % (result.returncode, script))
