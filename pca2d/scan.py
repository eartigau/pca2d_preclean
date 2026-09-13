"""What is in a data root, remembered between visits.

Reading one header of every spectrum of a campaign takes a few seconds the first
time and nothing at all afterwards: the answers are kept in ONE file per data
root, keyed by file name with its size and its modification time, so the next
visit reads only what was added or replaced. That is what lets a window show the
median SNR and the median exposure time of a target beside its name, which is
how anybody decides whether a target is worth a run, without waiting for a
campaign to be opened again.

Nothing is ever written inside the data root. A root can be a read-only archive
or a shared disk (/Volumes/irrisor here), and a pipeline that starts leaving
files in the place its inputs live is one mistake away from corrupting them. The
index goes under the user's home instead, one file per root, named after the
root so that it can be found and deleted by hand.

Only headers are read: the flux arrays are never touched, so a scan costs one
header read per file rather than the 40 MB the file holds. The keywords are the
pipeline's own (pca2d.tfits): EXTSNxxx per order for the SNR, EXPTIME for the
exposure time, INSTRUME for the instrument, which is read from the FILE and
never from a configuration.
"""
from __future__ import annotations

import datetime
import glob
import hashlib
import json
import os

import numpy as np

#: where the indexes live, outside every data root
INDEX_HOME = os.path.expanduser("~/.pca2d/scans")
#: bumped only when the SHAPE of an index changes. A new per-file field does not
#: need it: FIELDS below re-reads the files that lack it and keeps the rest.
VERSION = 1
#: every field a cached record must carry to be taken as it is. A record from
#: before a field existed is re-read rather than believed, and the campaign's
#: other numbers stay on screen while that happens: throwing the whole index
#: away instead would cost the minutes of a first scan for one new column.
FIELDS = ("snr", "exptime", "mjd", "instrument", "mag", "mag_band",
          "berv", "ecl_lat")

#: The brightness each pipeline writes, and WHICH BAND it is, which is not the
#: same for the two: NIRPS writes the J magnitude of the ESO target package,
#: SPIRou writes OBJMAG, which is H. Checked against SIMBAD on both SPIRou
#: campaigns here: TOI-2120 OBJMAG 10.45 against H 10.452 (J is 11.100) and
#: TOI-1452 10.03 against H 10.026 (J 10.604). So the band travels with the
#: number and is shown beside it; one column holding J for some targets and H
#: for others, unlabelled, would be wrong by a magnitude and say nothing about it.
MAG_KEYS = (("ESO OCS TARG JMAG", "J"), ("OBJMAG", "H"), ("HMAG", "H"),
            ("JMAG", "J"))


def index_path(root, home=None):
    """The index file for one data root, under the user's home.

    The name carries the root's last folder so a human can tell the files
    apart, and a digest of its full path so two roots ending in `science`
    cannot share one index.
    """
    full = os.path.abspath(os.path.expanduser(root or "."))
    digest = hashlib.sha1(full.encode("utf-8")).hexdigest()[:10]
    stem = os.path.basename(full.rstrip(os.sep)) or "root"
    return os.path.join(home or INDEX_HOME, "%s_%s.json" % (stem, digest))


def load(root, home=None):
    """The index of this root, or an empty one if there is none to read.

    A file written by an older version, or one that cannot be parsed, is
    treated as absent: an index is a cache, and rebuilding it costs seconds.
    """
    path = index_path(root, home)
    try:
        with open(path) as handle:
            index = json.load(handle)
    except (OSError, ValueError):
        index = {}
    if index.get("version") != VERSION or not isinstance(index.get("objects"),
                                                         dict):
        index = {"version": VERSION, "objects": {}}
    index["root"] = os.path.abspath(os.path.expanduser(root or "."))
    return index


def save(index, root, home=None):
    """Write the index, and return where it went (None if it could not be)."""
    path = index_path(root, home)
    index["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as handle:
            json.dump(index, handle, indent=1, sort_keys=True)
    except OSError:
        return None
    return path


def stamp_of(path):
    """(size, modification time) of a file: what says it is the same file."""
    info = os.stat(path)
    return [int(info.st_size), round(float(info.st_mtime), 3)]


def scan_file(path):
    """One spectrum's headers: instrument, band SNR, exposure time, date.

    Headers only. The SNR is the median of the per-order extraction SNR APERO
    wrote (EXTSNxxx), which is pca2d.tfits.band_snr over every order rather
    than over the domain's: picking the domain's orders would mean reading the
    wavelength solution, and the point of this is to be cheap. It is a number
    to sort targets by, not a number to publish.
    """
    from astropy.io import fits

    from .tfits import INSTRUMENTS

    out = {"instrument": "?", "snr": None, "exptime": None, "mjd": None,
           "mag": None, "mag_band": None, "berv": None, "ecl_lat": None}
    with fits.open(path, memmap=True) as hdulist:
        instrument = ""
        for hdu in hdulist:
            value = hdu.header.get("INSTRUME")
            if value:
                instrument = str(value).strip().upper()
                break
        out["instrument"] = instrument or "?"
        # the science extension first, PRIMARY second: APERO writes BJD and
        # EXTSNxxx on the science extension alone (pca2d.tfits)
        names = [str(hdu.name).upper() for hdu in hdulist]
        flux = str(INSTRUMENTS.get(instrument, {}).get("flux", "")).upper()
        headers = []
        if flux and flux in names:
            headers.append(hdulist[names.index(flux)].header)
        headers.append(hdulist[0].header)
        # EXTSN000..EXTSN074 are looked for by prefix rather than by order
        # number: the number of orders is in the data, and the data is what
        # this avoids reading
        snr = [_number(value) for key, value in headers[0].items()
               if key.startswith("EXTSN")]
        snr = [v for v in snr if v is not None and v > 0]
        if snr:
            out["snr"] = round(float(np.median(snr)), 2)
        for name, keys in (("exptime", ["EXPTIME"]),
                           ("mjd", ["MJDMID", "MJDATE", "BJD"]),
                           # what the star's motion through the observer's frame
                           # amounts to over a campaign, which is what decides
                           # whether the correction helps at all
                           ("berv", ["BERV"])):
            for key in keys:
                value = next((_number(h[key]) for h in headers if key in h), None)
                if value is not None:
                    out[name] = round(value, 6)
                    break
        # the ecliptic latitude, which fixes the BERV a target can EVER have:
        # |BERV| <= 29.78 cos(beta) km/s. TOI-1452 sits at +80.5 deg, so its
        # whole possible span is 9.8 km/s and no amount of observing will
        # separate the two frames for it. SPIRou writes RA_DEG/DEC_DEG, NIRPS
        # writes RA/DEC already in degrees.
        ra = next((_number(h[k]) for k in ("RA_DEG", "RA")
                   for h in headers if k in h and _number(h[k]) is not None), None)
        dec = next((_number(h[k]) for k in ("DEC_DEG", "DEC")
                    for h in headers if k in h and _number(h[k]) is not None), None)
        if ra is not None and dec is not None and abs(dec) <= 90.0:
            out["ecl_lat"] = round(_ecliptic_latitude(ra, dec), 3)
        for key, band in MAG_KEYS:
            value = next((_number(h[key]) for h in headers if key in h), None)
            # 0.0 is what NIRPS writes in the fields it did not fill, and a
            # guide magnitude of -9999.9 is a missing one, not a bright star
            if value is not None and 0.0 < value < 30.0:
                out["mag"], out["mag_band"] = round(value, 3), band
                break
    return out


#: obliquity of the ecliptic, degrees (J2000)
OBLIQUITY = 23.439291
#: the Earth's mean orbital speed, km/s
EARTH_SPEED = 29.78


def _ecliptic_latitude(ra_deg, dec_deg):
    """Ecliptic latitude from equatorial coordinates, in degrees.

    The rotation by the obliquity, done here rather than through a coordinate
    library: this is called once per spectrum during a scan, and the answer only
    has to be good to a degree to say what a target's BERV can reach.
    """
    ra = np.radians(float(ra_deg))
    dec = np.radians(float(dec_deg))
    eps = np.radians(OBLIQUITY)
    return float(np.degrees(np.arcsin(np.sin(dec) * np.cos(eps)
                                      - np.cos(dec) * np.sin(eps) * np.sin(ra))))


def berv_limit(ecl_lat):
    """The largest |BERV| this target can ever have, km/s, from its latitude."""
    if ecl_lat is None:
        return None
    return float(EARTH_SPEED * np.cos(np.radians(float(ecl_lat))))


def _number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def objects_of(root, pattern="*t.fits"):
    """{object: [file names]} for every subfolder of the root holding spectra."""
    out = {}
    if not root or not os.path.isdir(root):
        return out
    for name in sorted(os.listdir(root)):
        folder = os.path.join(root, name)
        if not os.path.isdir(folder):
            continue
        files = sorted(os.path.basename(p)
                       for p in glob.glob(os.path.join(folder, pattern)))
        if files:
            out[name] = files
    return out


def usable(record, stamp):
    """Whether a cached record can be taken as it is, rather than re-read.

    The same file, by size and modification time, AND carrying every field the
    current code shows. A file that could not be read is remembered as such so
    it is not re-read at every visit.
    """
    if not isinstance(record, dict) or record.get("stamp") != stamp:
        return False
    if record.get("unreadable"):
        return True
    return all(field in record for field in FIELDS)


def spread(items, first=10):
    """The same items, reordered so that ANY prefix samples the whole list.

    File names sort by date, so the first ten of a campaign are its first night,
    and a median taken from them is that night's weather rather than the
    campaign's: on GJ~1 the first ten spectra give a signal-to-noise of 134
    against the campaign's 163, 18% low. Reading in a comb instead, offset 0 of
    every `step`, then offset 1, and so on, costs exactly the same total reads
    and makes the estimate shown after ten files an estimate OF THE CAMPAIGN.
    """
    items = list(items)
    n = len(items)
    if n <= max(1, first):
        return items
    step = max(1, n // first)
    order = [i for offset in range(step) for i in range(offset, n, step)]
    return [items[i] for i in order]


def update(root, index=None, pattern="*t.fits", objects=None, on_file=None,
           home=None, on_listed=None, on_object=None, every=10):
    """Bring the index level with the root, reading only what changed.

    A file already in the index with the same size and modification time is
    taken as it was; one that is new or was replaced is read; one that is gone
    is dropped, as is an object folder that is. `objects` limits the work to
    some of them.

    Three callbacks, so that a window shows a campaign's properties from its
    first few spectra instead of waiting for all of them:
      `on_listed(found)`  once the folders are listed, before anything is read
      `on_file(object, done, total)`  as files are read
      `on_object(name, known, total)`  every `every` files AND when the object is
                          finished. `known` is how many of its spectra the index
                          now holds numbers for and `total` how many the folder
                          has, so known < total says "this is an estimate from
                          the first few", which is what a caller shows with a
                          tilde. The index is saved at those points too, so a
                          scan interrupted halfway keeps what it read.

    Returns (index, {"read": n, "kept": n, "gone": n}).
    """
    index = index if index is not None else load(root, home)
    table = index.setdefault("objects", {})
    found = objects_of(root, pattern)
    if on_listed is not None:
        on_listed({name: len(files) for name, files in found.items()})
    tally = {"read": 0, "kept": 0, "gone": 0}
    for name in list(table):
        if name not in found and (objects is None or name in objects):
            del table[name]
    # What each object needs read, before any of it is read: the rows can then
    # show what is already known, all of them at once.
    plans = {}
    for name, files in found.items():
        if objects is not None and name not in objects:
            continue
        entry = table.setdefault(name, {})
        cached = entry.get("files") or {}
        fresh = {}
        todo = []
        for filename in files:
            path = os.path.join(root, name, filename)
            try:
                stamp = stamp_of(path)
            except OSError:
                continue
            old = cached.get(filename)
            if usable(old, stamp):
                fresh[filename] = old
                tally["kept"] += 1
            else:
                # The SAME file, only missing a field this version shows: it
                # stays on screen, as it was, while it is read again. Dropped
                # instead, adding one column emptied every row of everything it
                # already knew, and the object being re-read lost its instrument,
                # its SNR and its count until the last of its spectra was back
                # (2026-09-13, GJ 1 and GL699_SPIROU showing no instrument).
                if isinstance(old, dict) and old.get("stamp") == stamp:
                    fresh[filename] = old
                todo.append((filename, path, stamp))
        tally["gone"] += len([f for f in cached if f not in fresh])
        # what the row can already show, before a single file of it is read
        entry["files"] = fresh
        if on_object is not None and fresh:
            on_object(name, len(fresh), len(files))
        plans[name] = {"entry": entry, "fresh": fresh, "total": len(files),
                       # read in a comb, so the numbers shown after ten files
                       # describe the campaign and not its first night
                       "todo": spread(todo, every or 10), "done": 0,
                       # fixed here: counting what REMAINS would shrink the
                       # total as the batch is consumed
                       "to_read": len(todo)}

    # A ROUND at a time, `every` files of each object, rather than one object
    # from end to end. With a dozen campaigns the last of them used to wait for
    # the eleven before it; now every row carries an estimate within the first
    # round, and they all sharpen together (the user, 2026-09-13: "ça permet de
    # voir où on va avant de finir la première").
    step = max(1, int(every or 10))
    while any(plan["todo"] for plan in plans.values()):
        for name, plan in plans.items():
            batch, plan["todo"] = plan["todo"][:step], plan["todo"][step:]
            if not batch:
                continue
            for filename, path, stamp in batch:
                try:
                    record = scan_file(path)
                except Exception:                                 # noqa: BLE001
                    # a file that cannot be read is recorded as unreadable
                    # rather than re-read at every visit, fields left empty
                    record = {"instrument": "?", "snr": None, "exptime": None,
                              "mjd": None, "mag": None, "mag_band": None,
                              "berv": None, "ecl_lat": None, "unreadable": True}
                record["stamp"] = stamp
                plan["fresh"][filename] = record
                plan["done"] += 1
                tally["read"] += 1
                if on_file is not None:
                    on_file(name, plan["done"], plan["to_read"])
            plan["entry"]["files"] = plan["fresh"]
            if on_object is not None:
                on_object(name, len(plan["fresh"]), plan["total"])

    for name, plan in plans.items():
        plan["entry"]["scanned"] = datetime.datetime.now().isoformat(
            timespec="seconds")
    return index, tally


def summary(index, name):
    """What a window shows on one row: counts, medians, instrument, dates.

    `files` is how many spectra the index holds numbers for, which during a scan
    is fewer than the folder has: the caller knows the folder's count and says
    so, and marks the medians as estimates.
    """
    files = ((index.get("objects") or {}).get(name) or {}).get("files") or {}
    out = {"object": name, "files": len(files), "snr": None, "exptime": None,
           "instrument": "?", "first": None, "last": None,
           "mag": None, "mag_band": None}
    bands = [f.get("mag_band") for f in files.values()
             if isinstance(f, dict) and f.get("mag_band")]
    if bands:
        # the band is a property of the pipeline, so one campaign has one band
        out["mag_band"] = max(set(bands), key=bands.count)
    for field in ("snr", "exptime", "mag"):
        values = [f[field] for f in files.values()
                  if isinstance(f, dict) and f.get(field) is not None]
        if values:
            out[field] = float(np.median(values))
    instruments = sorted({f.get("instrument") for f in files.values()
                          if isinstance(f, dict) and f.get("instrument")
                          and f.get("instrument") != "?"})
    # a folder holding two instruments is a folder to look at, so say both
    # rather than the first one read
    out["instrument"] = "+".join(instruments) if instruments else "?"
    dates = [f["mjd"] for f in files.values()
             if isinstance(f, dict) and f.get("mjd") is not None]
    if dates:
        out["first"], out["last"] = min(dates), max(dates)
    return out


def summaries(index):
    """One summary per object, by name."""
    return [summary(index, name) for name in sorted(index.get("objects") or {})]


#: the width of a BERV bin, in km/s
BERV_BIN = 3.0


def _ecl_lat(index, name):
    """The object's ecliptic latitude, from the first spectrum that has one."""
    files = ((index.get("objects") or {}).get(name) or {}).get("files") or {}
    for record in files.values():
        if isinstance(record, dict) and record.get("ecl_lat") is not None:
            return float(record["ecl_lat"])
    return None


def bervs_of(index, name):
    """Every barycentric velocity this object was observed at, in km/s."""
    files = ((index.get("objects") or {}).get(name) or {}).get("files") or {}
    return np.array([f["berv"] for f in files.values()
                     if isinstance(f, dict) and f.get("berv") is not None],
                    dtype=float)


def berv_coverage(index, names, width=BERV_BIN):
    """What the chosen stars' barycentric coverage actually is.

    Returns (edges, {object: counts per bin}, summary), the bins `width` km/s
    wide and shared by every object so the bars stack.

    `summary` carries, for the union and for each object, the span (max minus
    min) and the EFFECTIVE coverage, which is the number of bins that hold at
    least one exposure times the bin width. The two differ exactly where it
    matters: a campaign observed at two extremes and nowhere between has a wide
    span and little coverage, and it is the coverage that says how well the two
    frames can be told apart (docs, the barycentric span decides).
    """
    per_object = {n: bervs_of(index, n) for n in names}
    per_object = {n: v for n, v in per_object.items() if v.size}
    if not per_object:
        return np.array([]), {}, {"span": 0.0, "effective": 0.0, "objects": {}}
    every = np.concatenate(list(per_object.values()))
    lo = float(np.floor(every.min() / width) * width)
    hi = float(np.ceil(every.max() / width) * width)
    edges = np.arange(lo, hi + width, width)
    counts = {n: np.histogram(v, bins=edges)[0] for n, v in per_object.items()}
    filled = np.zeros(len(edges) - 1, dtype=bool)
    for c in counts.values():
        filled |= c > 0
    limits = [berv_limit(_ecl_lat(index, n)) for n in per_object]
    limits = [v for v in limits if v is not None]
    summary = {"span": float(every.max() - every.min()),
               "effective": float(filled.sum() * width),
               "n": int(every.size),
               # twice the largest amplitude any of them can reach: what the
               # sky allows, against what the campaign got
               "possible": 2.0 * max(limits) if limits else None,
               "objects": {n: {"span": float(v.max() - v.min()),
                               "effective": float((counts[n] > 0).sum() * width),
                               "n": int(v.size)}
                           for n, v in per_object.items()}}
    return edges, counts, summary


def nights_of(index, name):
    """The nights this object was observed on, as integer MJDs.

    The night is what the atmosphere belongs to, so this is what says whether
    two campaigns can share an observer basis at all (docs/joint_fit.md). A UT
    date boundary falls in the middle of a night at both sites, but a campaign's
    nights are counted the same way for every object, so the comparison holds.
    """
    files = ((index.get("objects") or {}).get(name) or {}).get("files") or {}
    return {int(np.floor(f["mjd"])) for f in files.values()
            if isinstance(f, dict) and f.get("mjd") is not None}


def shared_nights(index, names):
    """(nights common to all of them, [each one's own count]).

    Empty or nearly empty is the warning a joint fit needs: a basis fitted on
    stars that were not observed together is an average over conditions none of
    them met. Measured on PROXIMA+GJ1+GJ3090: two nights of 258, 147 and 99.
    """
    sets = [nights_of(index, name) for name in names]
    sets = [s for s in sets if s]
    if len(sets) < 2:
        return set(), [len(s) for s in sets]
    return set.intersection(*sets), [len(s) for s in sets]
