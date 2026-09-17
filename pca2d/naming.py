"""What to call a run, from the command it is.

Two things have to be read off a folder name in an output root months later:
WHAT was reduced, and whether it was reduced the same way as the one beside
it. The targets give the first, and no list of settings short enough to be a
folder name gives the second, so the command itself is hashed: same command,
same six characters; one flag different anywhere, a different six. Unique in
the only sense that matters here, which is that a clash is improbable rather
than impossible.

The command, and not a handful of settings: the roots a run reads and writes,
the stages it runs and every LBL setting change what it produces as surely as
the component counts do, and hashing ten keys left two different runs under
one name. Two things are taken out of it first, so that the name means the
same on both sides:

    the name itself   or it would hash the name it is computing
    the targets' order   the same set in another order is the same reduction

The window hashes the command it shows (`gui.build_command`); the pipeline
hashes its own `sys.argv`, which is what `--name auto` is.
"""

from __future__ import annotations

import hashlib
import re

#: the flags whose value is a list of targets
TARGETS = ("--objects",)
#: the flag left out of the hash
NAME = "--name"


def normalised(argv):
    """The command line a name is hashed from, as a list."""
    out = []
    skip = False
    for i, word in enumerate(list(argv or [])):
        if skip:
            skip = False
            continue
        if word == NAME:
            skip = True                      # its value goes with it
            continue
        if word.startswith(NAME + "="):
            continue
        flag = word.split("=", 1)[0]
        if flag in TARGETS:
            names = word.split("=", 1)[1] if "=" in word else None
            if names is None and i + 1 < len(argv):
                names, skip = argv[i + 1], True
            out += [flag, ",".join(sorted((names or "").split(",")))]
            continue
        out.append(word)
    return out


def targets_head(names, most=3):
    """The part of a name that says what was reduced."""
    names = sorted(str(n).strip() for n in (names or []) if str(n).strip())
    if not names:
        head = "run"
    elif len(names) <= most:
        head = "+".join(names)
    else:
        # three names is already a long folder name; past that, say how many
        head = "%s+%d" % ("+".join(names[:2]), len(names) - 2)
    return re.sub(r"[^0-9A-Za-z._+-]", "_", head)[:40].strip("_+") or "run"


def run_name(argv, names=None, digits=6):
    """'<targets>_<hash of the command>' for the command `argv`.

    `names` are the targets when the caller knows them; otherwise they are
    read from the command's own --object and --objects.
    """
    argv = list(argv or [])
    if names is None:
        names = []
        for i, word in enumerate(argv):
            flag, _, value = word.partition("=")
            if flag in TARGETS + ("--object",):
                value = value or (argv[i + 1] if i + 1 < len(argv) else "")
                names += [n for n in value.split(",") if n]
    payload = " ".join(normalised(argv))
    short = hashlib.blake2b(payload.encode("utf-8"),
                            digest_size=8).hexdigest()[:digits]
    return "%s_%s" % (targets_head(names), short)
