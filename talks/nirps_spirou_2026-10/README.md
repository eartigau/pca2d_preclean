# Two-frame precleaning: a talk for the NIRPS and SPIRou teams

`two_frame_precleaning.pdf`, 17 slides, October 2026.

## Rebuilding it

```sh
pdflatex talk.tex && pdflatex talk.tex
```

The figures in `figs/` come from two places:

- **The pipeline's own figure generator**, so the talk cannot show something
  the pipeline does not produce:

  ```sh
  python docs/make_figures.py --run outputs/TOI2120/1-3 \
      --cube cache/cube_tfits_1b1971abce62 --format pdf --out <figs>
  ```

- **Two scripts written for the talk**, both of which re-measure their numbers
  from the `.rdb` files at every run rather than hard-coding them, so a figure
  cannot quietly go stale:
  - `gain_vs_snr.py`: the gain against SNR over five campaigns, with paired
    bootstrap error bars.
  - `berv_talk.py`: the BERV-bias figure reshaped for a 16:9 slide. It calls
    `pca2d.texreport`'s own functions, so it cannot diverge from the report's
    version; only the titles are shortened, because the numbers are on the
    next slide.

## The numbers quoted

Everything comes from `outputs/TOI2120/1-3` (run of 2026-10-16) and the LBL
`.rdb` products in `lbl/lblrdb` and `outputs/lbl/lblrdb`. The five-campaign
table is measured on epochs the two series share, paired.

| Star | Instrument | SNR | Epochs | rms in | rms out | Gain |
|---|---|---|---|---|---|---|
| TOI-2120 | SPIRou | 33 | 80 | 47.98 | 14.20 | 3.38x |
| TOI-1452 | NIRPS | 38 | 51 | 101.43 | 69.37 | 1.46x |
| GJ 328 | NIRPS | 181 | 46 | 9.34 | 9.03 | 1.03x |
| GJ 1 | NIRPS | 233 | 140 | 2.54 | 2.56 | 0.99x |
| Proxima | NIRPS | 287 | 250 | 2.96 | 3.29 | 0.90x |

Spearman rho = -1.000 over the five (p = 0.000 on the rank test, but five
points only); the power-law fit gives gain ~ SNR^-0.45.

**Proxima is to be redone** with the HARPS updates; that row is provisional.
