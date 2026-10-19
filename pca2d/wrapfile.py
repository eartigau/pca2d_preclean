"""The wrap file the window runs, now that pca2d lives inside LBL.

The window no longer runs the pipeline itself. What it was asked becomes an
ordinary LBL wrap script, the runparams dicts LBL users know, with a PCA2D
dict in the one that measures the corrected spectra; the script is run by the
python of the LBL that lbl.environment names (lbl-rapide by default), and the
correction is made there, by lbl.pca2d, called from lbl_wrap:

    DATA_DIR/science/<object>/           the delivered spectra, linked here
    DATA_DIR/science/<object>_pca2d/     the corrected ones, written by LBL
    DATA_DIR/pca2d/<object>/<M>-<N>/     the fit, the figures, the report
    DATA_DIR/run_<object>_pca2d.py       the script, to be read and re-run

A named run is <object>_pca2d_<name>, under DATA_DIR/pca2d/_<name>/, so two
of them never share an LBL object. Several objects ticked together are one
joint correction, PCA2D joint=True, against one observer basis.

What travels into PCA2D: the two counts, the settings the window changed from
the configuration, the object's own block of that configuration and the
variant, if one is named. Paths do not: LBL decides where everything goes.
"""

from __future__ import annotations

import glob
import os
import pprint
import re
import textwrap

from .logger import log

#: what the corrected object is called after the delivered one (lbl_pca2d)
SUFFIX = "_pca2d"
#: sections of the configuration that are pca2d's and may go into PCA2D;
#: `lbl` is the wrap's own business and `output` is where LBL puts things
PCA2D_SECTIONS = ("input", "domain", "target", "registration", "highpass",
                  "quality", "weights", "twoframe", "correct")
#: keys of those sections that are paths or the object's name, never sent
NOT_SENT = {("input", "directory"), ("input", "object"),
            ("input", "object_header"), ("input", "format"),
            ("input", "instrument")}
STEPS = ("template", "mask", "compute", "compile")


def corrected_name(obj: str, label: str = "") -> str:
    return obj + SUFFIX + ("_" + label if label else "")


def clean_label(text) -> str:
    return re.sub(r"[^0-9A-Za-z._+-]", "_", str(text or "").strip()).strip("_")


def sendable(block: dict) -> dict:
    """The part of a configuration block that PCA2D takes: pca2d's sections,
    less the paths and the names LBL decides."""
    out = {}
    for section, values in (block or {}).items():
        if section not in PCA2D_SECTIONS or not isinstance(values, dict):
            continue
        kept = {k: v for k, v in values.items() if (section, k) not in NOT_SENT}
        if kept:
            out[section] = kept
    return out


def merge(base: dict, more: dict) -> dict:
    out = {k: dict(v) for k, v in base.items()}
    for section, values in more.items():
        out.setdefault(section, {}).update(values)
    return out


def pca2d_dict(state, config_file, objects, defaults=None) -> dict:
    """The PCA2D dict these choices mean."""
    from .cli import load_variant
    from .config import read_yaml
    from .gui import variant_yaml

    sections = {}
    # the object's own block of the configuration, as a run of the pipeline
    # would have layered it: the first object's for a joint run
    with open(config_file) as handle:
        user = read_yaml(handle) or {}
    block = (user.get("objects") or {}).get(objects[0]) or {}
    sections = merge(sections, sendable(block))
    # the variant, the nominal plus its own lines
    variant = load_variant(config_file, state.get("variant")) \
        if state.get("variant") else None
    if variant:
        sections = merge(sections, sendable(variant))
    # what the window changed, last, as on the command line
    sections = merge(sections, sendable(variant_yaml(state, defaults)))
    for key, section in (("min_rjd", "quality"), ("max_rjd", "quality")):
        value = str(state.get(key) or "").strip()
        if value:
            sections.setdefault(section, {})[key] = float(value)
    out = {"n_star": int(state.get("n_star") or 1),
           "n_earth": int(state.get("n_earth") or 3)}
    # the two counts are switches of their own, not twoframe keys
    tf = sections.get("twoframe") or {}
    for key in ("n_star", "n_earth"):
        tf.pop(key, None)
    if "figures" in (state.get("stages") or []) or state.get("stage_figures") \
            is not None:
        out["figures"] = bool(state.get("stage_figures", True))
    label = clean_label(state.get("run_name"))
    if label:
        out["name"] = label
    if len(objects) > 1:
        out["joint"] = True
    out.update({s: v for s, v in sections.items() if v})
    return out


def lbl_steps(config: dict, state) -> list:
    """LBL's steps after the correction: none when the lbl stage is off."""
    if state.get("stage_lbl") is False:
        return []
    steps = [str(s).strip().lower() for s in
             ((config.get("lbl") or {}).get("steps") or STEPS)]
    return [s for s in STEPS if s in steps]


def base_params(config, data_dir, objects, teffs, steps, tables) -> dict:
    from .lbl import profile

    instrument, data_source, _ = profile(config)
    params = {
        "INSTRUMENT": instrument,
        "DATA_SOURCE": data_source,
        "DATA_DIR": data_dir,
        "DATA_TYPES": ["SCIENCE"] * len(objects),
        "OBJECT_SCIENCE": list(objects),
        "OBJECT_COMPARISON": list(objects),
        "OBJECT_TEFF": list(teffs),
    }
    for step in STEPS:
        params["RUN_LBL_%s" % step.upper()] = step in steps
        params["SKIP_LBL_%s" % step.upper()] = True
    if tables:
        params["RESPROJ_TABLES"] = dict(tables)
    params["REPORT"] = bool((config.get("lbl") or {}).get("report", True))
    return params


SCRIPT = '''#!%(python)s
"""%(title)s

Written by the pca2d-preclean window on %(when)s, and left here to be read,
edited and run again by hand. An ordinary LBL wrap script: the PCA2D dict is
what makes LBL correct the spectra first (lbl.recipes.lbl_pca2d), and every
key of it is a key of the pca2d configuration. Run by the LBL it was written
for:

    %(python)s %(script)s
"""

from lbl.recipes import lbl_wrap
%(blocks)s


if __name__ == "__main__":
%(main)s
'''


def _block(name, comment, params) -> str:
    lines = ["", ""] + ["# " + line for line in textwrap.wrap(comment, 76)]
    lines.append("%s = dict(" % name)
    for key, value in params.items():
        text = pprint.pformat(value, width=60, sort_dicts=False)
        text = text.replace("\n", "\n" + " " * (5 + len(key)))
        lines.append("    %s=%s," % (key, text))
    lines.append(")")
    return "\n".join(lines)


def beside_config(state, path):
    """A path as the run reads it: relative to the configuration's folder,
    which is the rule the pipeline has always followed, not to wherever the
    window was started from."""
    if not path:
        return path
    base = os.path.dirname(os.path.abspath(state.get("config") or "config.yaml"))
    return os.path.abspath(os.path.join(base, os.path.expanduser(str(path))))


def lbl_dir(state, config) -> str:
    """LBL's DATA_DIR: the LBL folder of the window, else lbl/ under the
    output root, as lbl.directory has always been resolved."""
    typed = str(state.get("lbl_dir") or "").strip() or \
        (config.get("lbl") or {}).get("directory")
    if typed:
        return beside_config(state, typed)
    root = str(state.get("out_dir") or "").strip() or \
        (config.get("output") or {}).get("directory") or "outputs"
    return os.path.join(beside_config(state, root), "lbl")


def script_path(state, config=None) -> str:
    """Where the wrap file of these choices goes, worked out without reading
    any data: the command box shows it on every keystroke."""
    objects = list(state.get("objects") or [])
    stem = "+".join(objects) if len(objects) > 1 else (objects or ["?"])[0]
    name = corrected_name(stem, clean_label(state.get("run_name")))
    return os.path.join(lbl_dir(state, config or {}), "run_%s.py" % name)


def write(state, defaults=None, dry=False) -> dict:
    """Stage the spectra and write the script. Returns {script, python, env,
    data_dir, objects, report_dirs, text}; with dry=True nothing is written
    or linked and `text` is what would be."""
    import datetime

    from .config import load_config, spectra_dir
    from .gui import variant_yaml
    from .lbl import (chosen, dtemp_table, environment_for, link_mode,
                      link_spectra, resolve_teff)

    objects = list(state.get("objects") or [])
    if not objects:
        raise SystemExit("no object ticked")
    config_file = os.path.abspath(state.get("config") or "config.yaml")
    base = load_config(config_file)
    root = beside_config(state, str(state.get("data_dir") or "").strip()
                         or (base.get("input") or {}).get("directory") or "data")
    configs = [load_config(config_file, object_name=o, data_dir=root)
               for o in objects]
    config = configs[0]
    # the LBL settings the window changed, over the configuration's
    config.setdefault("lbl", {}).update(
        variant_yaml(state, defaults).get("lbl") or {})
    data_dir = lbl_dir(state, config)
    files = {o: sorted(glob.glob(os.path.join(
        spectra_dir(c), c["input"].get("pattern", "*t.fits"))))
        for o, c in zip(objects, configs)}
    teffs = [resolve_teff(c, files[o])[0] for o, c in zip(objects, configs)]
    steps = lbl_steps(config, state)
    if steps and not (config.get("lbl") or {}).get("run", True):
        # lbl.run off: the correction is made, the velocities are not
        steps = []
    tables = dtemp_table(config.get("lbl") or {}, teffs[0])
    pca2d = pca2d_dict(state, config_file, objects, defaults)
    label = pca2d.get("name", "")
    code = chosen(config)

    runs = []
    if steps and state.get("lbl_before", (config.get("lbl") or {}).get(
            "before", True)):
        runs.append(("BEFORE", "The spectra as delivered, measured against"
                     " the template LBL builds from them.",
                     base_params(config, data_dir, objects, teffs, steps,
                                 tables)))
    after = base_params(config, data_dir, objects, teffs, steps, tables)
    after["PCA2D"] = pca2d
    runs.append(("AFTER", "The same objects corrected first by pca2d (%s),"
                 " then measured against the template LBL builds from the"
                 " corrected spectra%s." % (", ".join(
                     corrected_name(o, label) for o in objects),
                     "" if steps else "; no LBL step is asked, so only the"
                     " correction is made"), after))

    stem = "+".join(objects) if len(objects) > 1 else objects[0]
    script = os.path.join(data_dir, "run_%s.py" % corrected_name(stem, label))
    text = SCRIPT % {
        "python": code["python"] or "/usr/bin/env python",
        "title": "LBL on %s, delivered and precleaned by pca2d." % ", ".join(objects),
        "when": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "script": script,
        "blocks": "".join(_block(*run) for run in runs),
        "main": "\n".join("    lbl_wrap.main(%s)" % name for name, _, _ in runs),
    }
    pca2d_root = os.path.join(data_dir, "pca2d", "_" + label if label else "")
    tag = "%d-%d%s" % (pca2d["n_star"], pca2d["n_earth"],
                       "v" if (pca2d.get("twoframe") or {}).get("velocity_term")
                       else "")
    if len(objects) > 1:
        from .joint import joint_name
        report_dirs = [os.path.join(pca2d_root, "joint", joint_name(objects), tag)]
    else:
        report_dirs = [os.path.join(pca2d_root, objects[0], tag)]
    out = {"script": script, "python": code["python"], "code": code,
           "data_dir": data_dir, "objects": objects, "text": text,
           "report_dirs": report_dirs,
           "env": environment_for(code["python"]) if code["python"] else None}
    if dry:
        return out
    mode, why = link_mode(data_dir, (config.get("lbl") or {}).get("link"))
    if why:
        log(why, "warn")
    for obj in objects:
        target = os.path.join(data_dir, "science", obj)
        linked, kept, strangers = link_spectra(files[obj], target, mode)
        log("%-24s %d spectra (%d new %s, %d already there)"
            % (obj, linked + kept, linked, mode + "s", kept), "value")
        if strangers:
            log("  %d other .fits are in %s and LBL takes the whole folder:"
                " %s%s" % (len(strangers), target, ", ".join(strangers[:3]),
                           ", ..." if len(strangers) > 3 else ""), "warn")
    os.makedirs(data_dir, exist_ok=True)
    with open(script, "w") as handle:
        handle.write(text)
    os.chmod(script, 0o755)
    log("wrap file   %s" % script, "value")
    return out


def report_pdf(report_dirs) -> str | None:
    """The newest PDF a run left in its folder, or None."""
    found = []
    for folder in report_dirs or []:
        found += glob.glob(os.path.join(folder, "*.pdf"))
    return max(found, key=os.path.getmtime) if found else None
