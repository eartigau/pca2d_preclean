"""Les exposants telluriques de LBL: l'eau, le sec, et le sec contre la masse d'air.

`lbl_telluclean` ajuste deux exposants par pose et les ecrit dans l'en-tete du
spectre nettoye: TLPEH2O pour l'eau, TLPEOTR pour tout le reste (O2, CO2, CH4
ensemble, et non par espece chimique). Cette page les montre.

    page 1  les deux series temporelles, l'eau et le sec, et la masse d'air
            sous elles pour comparaison
    page 2  le sec CONTRE la masse d'air, et le residu en dessous

CE QUE LA PAGE 2 VERIFIE. Les especes seches sont bien melangees dans
l'atmosphere, donc leur colonne suit la masse d'air et le nuage doit tomber sur
une droite. MAIS PAS FORCEMENT SUR LA PREMIERE BISSECTRICE: TAPAS est calcule
pour une colonne de reference, et un observatoire en altitude en a moins (La
Silla est a 2400 m), ce qui deplace l'ordonnee a l'origine sans toucher a la
pente. La droite ajustee est donc tracee avec SA pente et SON offset, les deux
libres, et la pente 1 passant par zero est rappelee a cote pour l'oeil. Le
residu est pris contre la droite ajustee, qui est le modele qu'on teste.

Si le sec a ete FIXE a la masse d'air a l'ajustement, ce que la configuration
livree fait, le nuage est exactement sur la droite et le residu est nul: la
page le dit alors plutot que de laisser croire a une mesure.

    python telluric_exponents.py --science DIR --out PAGE.pdf
"""
from __future__ import annotations

import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

RJD0 = 2400000.0
#: l'ecart sous lequel on considere que le sec n'a pas ete ajuste du tout
FIXED_TOL = 1e-9


def read_exponents(directory, modes=None):
    """(rjd, eau, sec, masse d'air, mode) de chaque spectre nettoye du dossier.

    `modes` filtre sur HIERARCH ESO INS MODE: deux modes d'un meme
    spectrographe sont deux instruments, et les melanger sur une figure fait
    croire a une dispersion qui n'est que la difference des deux.
    """
    rows = []
    for path in sorted(glob.glob(os.path.join(directory, "*.fits"))):
        try:
            header = fits.getheader(path)
        except Exception:                                     # noqa: BLE001
            continue
        if "TLPEH2O" not in header:
            continue
        mode = str(header.get("HIERARCH ESO INS MODE", "")).strip()
        if modes and mode not in modes:
            continue
        bjd = header.get("HIERARCH ESO QC BJD")
        rows.append((
            float(bjd) - RJD0 if bjd is not None else np.nan,
            float(header["TLPEH2O"]),
            float(header.get("TLPEOTR", np.nan)),
            float(header.get("HIERARCH ESO TEL AIRM START", np.nan)),
            mode,
        ))
    if not rows:
        return None
    return dict(
        rjd=np.array([r[0] for r in rows]),
        water=np.array([r[1] for r in rows]),
        dry=np.array([r[2] for r in rows]),
        airmass=np.array([r[3] for r in rows]),
        mode=np.array([r[4] for r in rows]),
    )


def dry_is_fixed(data):
    """Le sec a-t-il ete fixe a la masse d'air plutot qu'ajuste?"""
    good = np.isfinite(data["dry"]) & np.isfinite(data["airmass"])
    if not good.any():
        return False
    return bool(np.max(np.abs(data["dry"][good] - data["airmass"][good]))
                < FIXED_TOL)


def fit_line(x, y):
    """(pente, ordonnee a l'origine) des moindres carres, les DEUX libres.

    L'ordonnee est libre parce qu'un observatoire en altitude voit moins de
    colonne seche que la reference de TAPAS: cet offset est physique et doit
    pouvoir apparaitre au lieu d'etre absorbe par la pente.
    """
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 3:
        return np.nan, np.nan
    design = np.vstack([x[good], np.ones(int(good.sum()))]).T
    slope, intercept = np.linalg.lstsq(design, y[good], rcond=None)[0]
    return float(slope), float(intercept)


def draw(data, out):
    fixed = dry_is_fixed(data)
    with PdfPages(out) as pdf:
        # ---- page 1: les series temporelles
        fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        axes[0].plot(data["rjd"], data["water"], "o", ms=4, color="C0")
        axes[0].set_ylabel("exposant H2O")
        axes[0].set_title("la colonne d'eau, ajustee pose par pose"
                          " (TLPEH2O)", fontsize=10)
        axes[1].plot(data["rjd"], data["dry"], "o", ms=4, color="C1")
        axes[1].set_ylabel("exposant sec")
        axes[1].set_title("tout le reste d'un bloc, non par espece chimique"
                          " (TLPEOTR)%s"
                          % ("   FIXE a la masse d'air" if fixed else ""),
                          fontsize=10)
        axes[2].plot(data["rjd"], data["airmass"], "o", ms=4, color="0.4")
        axes[2].set_ylabel("masse d'air")
        axes[2].set_xlabel("rjd")
        for ax in axes:
            ax.grid(alpha=0.3)
        fig.suptitle("Les exposants telluriques de lbl_telluclean, %d poses"
                     % len(data["rjd"]))
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # ---- page 2: le sec contre la masse d'air, et le residu
        fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1]))
        x, y = data["airmass"], data["dry"]
        axes[0].plot(x, y, "o", ms=5, color="C1", label="les poses")
        span = np.linspace(np.nanmin(x) * 0.98, np.nanmax(x) * 1.02, 10)
        axes[0].plot(span, span, "k--", lw=1.0,
                     label="pente 1 par l'origine")
        slope, intercept = fit_line(x, y)
        if np.isfinite(slope):
            axes[0].plot(span, slope * span + intercept, "-", lw=1.4,
                         color="C3",
                         label="ajuste: %.3f x airmass %+.3f"
                               % (slope, intercept))
            resid = y - (slope * x + intercept)
        else:
            resid = np.full_like(y, np.nan)
        axes[0].set_ylabel("exposant sec")
        axes[0].legend(fontsize=8)
        axes[0].grid(alpha=0.3)
        axes[1].plot(x, resid, "o", ms=5, color="C3")
        axes[1].axhline(0.0, color="k", lw=1.0, ls=":")
        axes[1].set_xlabel("masse d'air")
        axes[1].set_ylabel("residu")
        axes[1].grid(alpha=0.3)
        if fixed:
            note = ("LE SEC EST FIXE a la masse d'air a l'ajustement: le nuage"
                    " est sur la droite par CONSTRUCTION\net le residu est nul."
                    " Cette page verifie la contrainte, elle ne la mesure pas.")
        else:
            rms = float(np.nanstd(resid, ddof=1))
            note = ("le sec est AJUSTE librement; residu contre la droite"
                    " ajustee: rms %.4f\nune ordonnee a l'origine non nulle est"
                    " attendue, l'observatoire etant en altitude" % rms)
        fig.suptitle("Le sec contre la masse d'air\n%s" % note, fontsize=10)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--science", required=True,
                        help="le dossier des spectres NETTOYES (<OBJET>_tc)")
    parser.add_argument("--out", required=True)
    parser.add_argument("--modes", default="HARPS",
                        help="les modes gardes, separes par des virgules;"
                             " vide pour tout garder")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    modes = [m for m in args.modes.split(",") if m] if args.modes else None
    data = read_exponents(args.science, modes)
    if data is None:
        raise SystemExit("aucun spectre portant TLPEH2O dans %s" % args.science)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    draw(data, args.out)
    print("%d poses -> %s" % (len(data["rjd"]), args.out))


if __name__ == "__main__":
    main()
