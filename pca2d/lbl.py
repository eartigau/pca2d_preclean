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
    # CADC, and not APERO, even for an APERO reduction: LBL's CADC classes are
    # the ones that read the named fibre extensions (FluxAB / WaveAB / BlazeAB
    # for SPIRou, FluxA / WaveA / BlazeA for NIRPS), which is where the
    # wavelength solution of a t.fits is. Its APERO classes read the primary
    # header's wavelength polynomials, which these files do not carry there.
    "SPIROU": ("SPIROU", "CADC"),
    "NIRPS": ("NIRPS_HE", "CADC"),
}

STEPS = ("telluclean", "template", "mask", "compute", "compile")

#: Where the target's effective temperature is written, in the order they are
#: tried. APERO puts OBJTEMP on the primary header and PP_TEFF on the science
#: extension, with PP_TEFFS saying where the DRS got it. LBL needs the number
#: and does not read it from anywhere itself: it stops with "Teff is require.
#: Please add OBJECT_TEFF to config". It is in the file, so read it there.
TEFF_KEYS = ("OBJTEMP", "PP_TEFF", "OBJ_TEMP", "TEFF")


def teff_from_header(path: str):
    """(teff, keyword) from the first of TEFF_KEYS that carries a number."""
    from .io import robust_open

    try:
        with robust_open(path) as hdulist:
            for key in TEFF_KEYS:
                for hdu in hdulist:
                    value = hdu.header.get(key)
                    if value in (None, "", "None"):
                        continue
                    try:
                        number = float(value)
                    except (TypeError, ValueError):
                        continue
                    if number > 0:
                        return number, key
    except Exception:                                         # noqa: BLE001
        return None, None
    return None, None


def resolve_teff(config: dict, files) -> tuple:
    """The Teff LBL will be given, and where it came from.

    `auto`, the default, reads it from the spectra, for the same reason the
    instrument is read from them: it is written in the file, by the pipeline
    that reduced it, and a number typed into a config is a number that can be
    typed wrong. A number in the config wins over the header, and is checked
    against it out loud, because disagreeing with the file is a thing worth
    doing deliberately and not by accident.
    """
    asked = (config.get("lbl") or {}).get("teff", "auto")
    header, key = (teff_from_header(files[0]) if files else (None, None))

    if isinstance(asked, str) and asked.strip().lower() == "auto":
        if header is None:
            log("lbl.teff is `auto` and none of %s is in the header of %s."
                " LBL needs a Teff to choose the stellar model its mask comes"
                " from: put a number in lbl.teff."
                % (", ".join(TEFF_KEYS), os.path.basename(files[0]) if files
                   else "the spectra"), "warn")
            return None, "nowhere"
        return header, "%s in the header" % key
    if asked is None:
        return None, "unset"
    value = float(asked)
    if header is not None and abs(header - value) > 1.0:
        log("lbl.teff says %.0f K and %s in the header says %.0f K. The config"
            " wins, which is what it is for, but one of the two is wrong."
            % (value, key, header), "warn")
    return value, "lbl.teff in the config"


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
                    data_source: str, teff=None) -> dict:
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
    if teff is not None:
        document["OBJECT_TEFF"] = teff
    return document


def runparams(config: dict, data_dir: str, instrument: str, data_source: str,
              objects: list, config_file: str, teff=None) -> dict:
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
        "OBJECT_TEFF": [teff] * len(objects),
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


def check_profile(instrument: str, data_source: str, config_file: str,
                  sample: str) -> tuple:
    """Ask LBL to read one spectrum with the chosen profile. (ok, why not).

    Two seconds against hours. The pairing of instrument and data source picks
    the class that reads the files, and the wrong one does not fail on the
    first line: it downloads a few hundred megabytes of stellar models, sorts
    every exposure, and dies at the reference file, which is what happened the
    first time this was pointed at SPIRou with LBL's APERO class instead of its
    CADC one. What is read here is the wavelength solution, because that is the
    thing the two classes disagree about.
    """
    try:
        from lbl.instruments import select
    except Exception as exc:                                  # noqa: BLE001
        return True, "LBL is not importable here (%s), so nothing was checked" % exc
    try:
        args = select.parse_args(
            ["INSTRUMENT", "DATA_DIR", "DATA_SOURCE", "DATA_TYPE", "INPUT_FILE"],
            dict(config_file=os.path.abspath(config_file)), __name__, parse=False)
        inst = select.load_instrument(args, plogger=None)
        image, header = inst.load_science_file(sample)
        inst.get_wave_solution(sample, image, header)
    except Exception as exc:                                  # noqa: BLE001
        return False, str(exc)
    return True, "reads %s" % os.path.basename(sample)


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

    teff, whence = resolve_teff(config, plan["files"])
    if teff is not None:
        log("Teff %.0f K, from %s" % (teff, whence), "value")

    document = config_document(config, data_dir, instrument, data_source, teff)
    config_file = os.path.join(plan["outdir"], "lbl_config.yaml")
    with open(config_file, "w") as handle:
        handle.write("# LBL's own configuration, written by pca2d-preclean.\n"
                     "# Every key here is one LBL knows; it refuses any other.\n")
        yaml.safe_dump(document, handle, sort_keys=False,
                       default_flow_style=False)

    ok, why = check_profile(instrument, data_source, config_file,
                            plan["files"][0] if plan["files"] else None)
    if ok:
        log("profile check: %s/%s %s" % (instrument, data_source, why), "value")
    else:
        log("%s/%s cannot read %s: %s"
            % (instrument, data_source, os.path.basename(plan["files"][0]), why),
            "error")
        log("the instrument and data source together pick the class that reads"
            " the files, and it is set per spectrograph in config.yaml. For a"
            " t.fits with named fibre extensions the reader is the CADC one.",
            "error")

    params = runparams(config, data_dir, instrument, data_source, objects,
                       config_file, teff)
    script = write_runner(os.path.join(plan["outdir"], "run_lbl.py"), params,
                          object_name, plan["tag"], config_file)

    if params["RUN_LBL_MASK"] and teff is None:
        log("no Teff, and the mask step is on. LBL stops on that rather than"
            " guessing, so either the spectra have to carry one of %s or"
            " lbl.teff has to be a number." % ", ".join(TEFF_KEYS), "warn")

    log("LBL config   %s" % config_file, "value")
    log("LBL script   %s" % script, "value")
    return {"config_file": config_file, "script": script, "objects": objects,
            "data_dir": data_dir, "readable": ok}


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
