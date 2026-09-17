# The optimisation megarun: how to run it, on rali

For the Claude instance driving this campaign on the workstation, and for EA
reading along from a laptop. The plan itself, what is to be understood and
what it is judged by, is `optimisation_megarun.md` at the root of this
repository. This file is the mechanics: where everything is, what to launch,
in what order, and what to write down after each run.

**The one rule: one run at a time.** A run is a fit and then LBL, it takes
between forty minutes and three hours, and two at once take longer than two
in a row and can fill the memory. Never launch a second while one is going.

## What is already there, verified on 2026-09-17

| what | where |
| --- | --- |
| project root | `/home/artigau/pca2d_optimisation` (45 TB free) |
| the code | `/data/spirou/pca2d_preclean`, a clone of this repository |
| the environment | `conda activate pca2d-preclean` (it holds LBL 0.67.008) |
| pdflatex | `/usr/bin/pdflatex`, so the reports compile |
| NIRPS spectra | `/cosmos99/nirps/apero-data/nirps_he_offline/objects/<OBJECT>` |
| SPIRou spectra | `/cosmos99/spirou/apero-data/spirou_offline/objects/<OBJECT>` |
| the station | rali, 40 cores, 376 GB of memory |

Where each target actually is, counted in `t.fits`. The plan asks for three
of them "with SPIRou"; they exist only in NIRPS, and GL406 exists in both:

| target | NIRPS | SPIRou |
| --- | --- | --- |
| GL406 | 411 | 679 |
| GJ1 | 292 | - |
| GJ707 | 175 | - |
| TOI4552 | 119 | - |
| TOIM4508 | 174 | - |
| TOI782 | 201 | - |
| GL725B | - | 971 |
| GL251 | - | 748 |
| GL48 | - | 989 |
| TOI2120 | - | 326 |
| TOI6091 | - | 245 |

## Step 0: the code, the environment, a screen that survives

```bash
git -C /data/spirou/pca2d_preclean pull
conda activate pca2d-preclean
cd /data/spirou/pca2d_preclean && python -m pytest -q      # a minute, all green
tmux new -s megarun                                        # or screen
```

Everything below is run inside that tmux session, so a dropped connection
does not stop a run.

## Step 1: stage the spectra, once

Per instrument, per object, symlinks only. Nothing is copied.

```bash
ROOT=/home/artigau/pca2d_optimisation
stage () {                       # stage <instrument> <source dir> <object>
  mkdir -p "$ROOT/data/$1/$3"
  find "$2/$3" -maxdepth 1 -name '*t.fits' ! -name '*_pp_*' \
       -exec ln -sf {} "$ROOT/data/$1/$3/" \;
  echo "$1 $3: $(ls "$ROOT/data/$1/$3" | wc -l) spectra"
}
NIRPS=/cosmos99/nirps/apero-data/nirps_he_offline/objects
SPIROU=/cosmos99/spirou/apero-data/spirou_offline/objects
for o in GL406 GJ1 GJ707 TOI4552 TOIM4508 TOI782; do stage NIRPS "$NIRPS" $o; done
for o in GL725B GL251 GL48 TOI2120 TOI6091 GL406; do stage SPIROU "$SPIROU" $o; done
```

The counts printed must match the table above. A count of 0 means the object
folder is spelled differently there: check with `ls $NIRPS | grep -i <name>`.

## Step 2: one configuration for the whole campaign

```bash
cd /home/artigau/pca2d_optimisation
cp /data/spirou/pca2d_preclean/config.yaml config.yaml
```

Then edit `config.yaml` so that the campaign writes where it should. The
`objects:` blocks of the repository's configuration are kept: they carry each
target's Teff and its published planets, which LBL and the report need.

```yaml
input:
  directory: /home/artigau/pca2d_optimisation/data/NIRPS   # per instrument
output:
  directory: /home/artigau/pca2d_optimisation/outputs
  cache_directory: /home/artigau/pca2d_optimisation/cache
  fits_directory: null            # everything under the output root
lbl:
  directory: /home/artigau/pca2d_optimisation/lbl
  run: true
```

`--data-dir` on the command line chooses the instrument for a run, so the
value above is only the default.

## Step 3: what one run looks like

```bash
cd /home/artigau/pca2d_optimisation
conda activate pca2d-preclean
mkdir -p logs

RUN=nominal                        # the name of this scenario
OBJ=GL406
INST=NIRPS
pca2d-preclean --object $OBJ \
    --config  /home/artigau/pca2d_optimisation/config.yaml \
    --data-dir /home/artigau/pca2d_optimisation/data/$INST \
    --out-dir  /home/artigau/pca2d_optimisation/outputs \
    --lbl-dir  /home/artigau/pca2d_optimisation/lbl \
    --no-fits-dir --name $RUN \
    --n-star 0 --n-earth 7 --weight velocity --high-pass 100 \
    --velocity-term false \
    2>&1 | tee logs/${RUN}_${INST}_${OBJ}.log
```

Several objects fitted together against one observer basis: replace
`--object GL406` with `--objects TOI4552,TOIM4508,TOI782`. The run then
writes into `outputs/_<name>/joint/<A+B+C>/<tag>/`.

A dry run prints everything it resolved and touches nothing:

```bash
pca2d-preclean ... --dry-run
```

**Check before every launch**: `df -h /home/artigau` (a run writes its
corrected spectra, as big as the target's own data) and that no other run is
going (`pgrep -af pca2d`).

## Step 4: the sweep, and the order it must be run in

Changing the two component counts, the velocity term, the correction's metric
or the shrinkage **reuses the cube**. Changing the high pass, the grid step or
the nightly coadding **builds a new one** (tens of minutes each). So: group
the runs by high pass, and inside a group by object set.

Phase A, on one small NIRPS target (GL406, 411 spectra) and one small SPIRou
target (TOI6091, 245), one parameter at a time around the nominal
(`--n-star 0 --n-earth 7 --weight velocity --high-pass 100 --velocity-term
false`). Name each run after what it changes, since the name is the folder:

| name | what changes |
| --- | --- |
| `nominal` | - |
| `w-flux` | `--weight flux` |
| `star1`, `star2` | `--n-star 1`, `--n-star 2` |
| `earth1`, `earth3`, `earth11` | `--n-earth 1`, `3`, `11` |
| `velterm` | `--velocity-term true` |
| `hp50`, `hp200` | `--high-pass 50`, `--high-pass 200` (new cube) |

Phase B, the sets of the plan, each alone and then together, at the best
settings Phase A found:

- NIRPS: `TOI4552,TOIM4508,TOI782` together, then each alone
- NIRPS: `GJ1,GJ707` together, then each alone
- SPIRou: `TOI2120,TOI6091` together, then each alone
- SPIRou: `GL725B,GL251,GL48` together, then each alone. GL725B is the one
  with almost no barycentric coverage (BERV within ±4.6 km/s), which is where
  the method is expected to have the least to work with

Phase C: the two or three settings that came out best in Phase A, run on
every target, so that the conclusion rests on more than one star.

## Step 5: after every run, score it and push

```bash
cd /data/spirou/pca2d_preclean
python megarun/score.py /home/artigau/pca2d_optimisation/outputs \
    --csv megarun/results.csv --status megarun/status.md
git add megarun/results.csv megarun/status.md && \
git commit -m "megarun: <what was run>, robust sigma <before> -> <after>" && \
git push
```

`megarun/status.md` is the table EA reads from the laptop after a `git pull`:
one row per star per run, the settings it was given, the robust sigma before
and after, and the gain. `results.csv` holds every number, for the figures at
the end.

Say in each commit message what was run and what it was worth. Those messages
are the campaign's log.

## Step 6: what to look at when something is surprising

- the run's own report, `outputs/.../<object>_<tag>.pdf`: the velocities, the
  V_tot bias with its ΔBIC, the phase-folded known planets, the periodograms
- a scan across component counts of one object, in one PDF:
  ```bash
  python -m pca2d.lblscan --object GL406 --tags 0-1 0-3 0-7 0-11 \
      --suffix '_PCA2D_{tag}' --outputs /home/artigau/pca2d_optimisation/outputs \
      --lbl-dir /home/artigau/pca2d_optimisation/lbl \
      --out /home/artigau/pca2d_optimisation/outputs/gl406_scan.pdf
  ```
- what is on the disks:
  ```bash
  python -m pca2d.housekeeping --config /home/artigau/pca2d_optimisation/config.yaml \
      --out-dir /home/artigau/pca2d_optimisation/outputs
  ```
- the runs already made, with their start times:
  `python -m pca2d.runs /home/artigau/pca2d_optimisation/outputs`

## Step 7: the paper, at the end

The plan asks for a full paper with references, and for the literature to be
read rather than remembered: use the Consensus tool for YARARA (Cretignier et
al.) and what else has been done on telluric and instrumental residuals in
near-infrared radial velocities, and cite what is found. Figures and tables
come from `megarun/results.csv` and from the runs' own reports. Write it in
`paper/`, beside the existing draft.

## What to tell EA, and when

- after each run: the commit of `megarun/status.md` is enough
- when a phase ends: a short summary in the commit message, with the numbers
- when something is stuck (a stage failing twice, the disk filling, a fit
  refusing to start): stop, write what happened in `megarun/status.md`, push,
  and say so rather than working around it
