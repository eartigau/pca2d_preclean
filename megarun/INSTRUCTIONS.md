# The optimisation megarun: how to run it, on rali

For the Claude instance driving this campaign on the workstation, and for EA
reading along from a laptop. The plan itself, what is to be understood and
what it is judged by, is `optimisation_megarun.md` at the root of this
repository. This file is the mechanics: where everything is, what to launch,
in what order, and what to write down after each run.

**Two rules.**

1. **One run at a time.** A run is a fit and then LBL, it takes between forty
   minutes and three hours, and two at once take longer than two in a row and
   can fill the memory. Never launch a second while one is going.
2. **Every run is named**, `--name <scenario>` or `--name auto`. A named run
   gets its own folder AND its own LBL object
   (`<target>_PCA2D_<counts>_<name>`); an unnamed one does not, so two
   scenarios of the same target with the same component counts would be one
   LBL object and the second would overwrite the first's velocities.

   Since 2026-09-18 (commit after `0807e4c`) the name is followed by the hash
   of the command line: `--name earth3` makes `_earth3_a1b2c3` and
   `GL406_PCA2D_0-7_earth3_a1b2c3`, and the compilation PDF is
   `GL406_0-7_a1b2c3.pdf`. `--name auto` is the targets and the hash, nothing
   typed. The hash covers everything the command says except how much of the
   run happens this time: `--stages`, `--dry-run`, `--rebuild-cube`,
   `--clean-cache` and `--lbl-before` are left out, so **a run resumed with
   `--stages lbl` keeps its name and lands in the folder it left.** The
   instrument is in the line through `--data-dir`, `--out-dir` and
   `--lbl-dir`, so the same target on two instruments is two runs.

   **Runs made before that stay where they are.** A folder `_<name>` already
   on disk is that run, and is used as it is, when it was made by the same
   command or before any hash was recorded; only one made by another command
   sends the new run to `_<name>_<hash>`. `pca2d.runs` and `megarun/score.py`
   read both, and give an old run the hash its own recorded command gives.
   The one change a resume of an old run sees: its PDF, written again, takes
   the new name with the hash beside the old one.

   `reuse_fit` (the correction-only variants) now reaches named runs: it takes
   a name (`nominal`, `earth3`, matching `_earth3` or `_earth3_<hash>`), a
   hash alone, or a folder, and the newest fit of what matches.

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

Where each target actually is, counted in `t.fits`. GL406 is the only one in
both, and it is run on each instrument separately: two instruments are two
runs, never one fit.

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
  directory: /home/artigau/pca2d_optimisation/outputs/NIRPS
  cache_directory: /home/artigau/pca2d_optimisation/cache
  fits_directory: null            # everything under the output root
lbl:
  directory: /home/artigau/pca2d_optimisation/lbl/NIRPS
  run: true
```

The three directories that carry the instrument are given on every command
line (`--data-dir`, `--out-dir`, `--lbl-dir`), so the values above are only
defaults. The cube cache is shared: its key already hashes the input
directory, so two instruments never meet in it.

## Step 3: what one run looks like

```bash
cd /home/artigau/pca2d_optimisation
conda activate pca2d-preclean
mkdir -p logs

ROOT=/home/artigau/pca2d_optimisation
RUN=nominal                        # the name of this scenario
OBJ=GL406
INST=NIRPS                         # or SPIROU
pca2d-preclean --object $OBJ \
    --config   $ROOT/config.yaml \
    --data-dir $ROOT/data/$INST \
    --out-dir  $ROOT/outputs/$INST \
    --lbl-dir  $ROOT/lbl/$INST \
    --no-fits-dir --name $RUN \
    --n-star 0 --n-earth 7 --weight velocity --high-pass 100 \
    --velocity-term false \
    2>&1 | tee logs/${RUN}_${INST}_${OBJ}.log
```

Several objects fitted together against one observer basis: replace
`--object GL406` with `--objects TOI4552,TOIM4508,TOI782`. The run then
writes into `outputs/<INST>/_<name>/joint/<A+B+C>/<tag>/`. Only objects of
the same instrument are ever fitted together: one observer basis is one
spectrograph's sky.

The output root and the LBL tree carry the instrument, and the name carries
the scenario, in the folder and in the LBL object alike. That is what keeps
GL406's two runs, one per instrument, and the sweep's scenarios apart: same
target, same counts, different velocities.

`--name auto` instead of a scenario name hashes the command line, which is
the same identifier the window proposes. Use it when a run is a one-off; for
the sweep, a name that says what it changes (`hp50`, `earth3`) is easier to
read in `megarun/status.md` months later.

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
- GL406 alone on NIRPS and alone on SPIRou, the same settings both times: the
  same star through two spectrographs, which says what belongs to the method
  and what belongs to the instrument

Phase C: the two or three settings that came out best in Phase A, run on
every target, so that the conclusion rests on more than one star.

## Step 5: after every run, score it and push

```bash
cd /data/spirou/pca2d_preclean
python megarun/score.py /home/artigau/pca2d_optimisation/outputs/NIRPS \
                       /home/artigau/pca2d_optimisation/outputs/SPIROU \
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
- several scenarios of one object against the delivered velocities, in one
  PDF. `--compare` takes the LBL object names themselves, which is what the
  per-scenario suffixes make them:
  ```bash
  python -m pca2d.lblscan --object GL406 \
      --compare GL406_PCA2D_0-1_earth1="1 observer" \
                GL406_PCA2D_0-3_earth3="3 observers" \
                GL406_PCA2D_0-7_nominal="7 observers" \
                GL406_PCA2D_0-11_earth11="11 observers" \
      --lbl-dir /home/artigau/pca2d_optimisation/lbl/NIRPS \
      --out /home/artigau/pca2d_optimisation/outputs/NIRPS/gl406_earth.pdf
  ```
- what is on the disks:
  ```bash
  python -m pca2d.housekeeping --config /home/artigau/pca2d_optimisation/config.yaml \
      --out-dir /home/artigau/pca2d_optimisation/outputs
  ```
- the runs already made, with their start times:
  `python -m pca2d.runs /home/artigau/pca2d_optimisation/outputs/NIRPS`

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
