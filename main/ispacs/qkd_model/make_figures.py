"""Talk figures for the BB84 key-rate model.

    fig_qkd_main.png   (a) SKR vs distance: no decoy / decoy finite key / ideal decoy bound
                       (b) QBER vs distance for the same three cases
    fig_qkd_range.png  maximum secure distance per scenario

Sweeps are cached in results/sweeps.npz; delete it (or pass --rerun) after changing the model.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import qkd_model as q

HERE = Path(__file__).resolve().parent
CACHE = HERE / "results" / "sweeps.npz"
FIG = HERE / "figures"
N_MAIN = 1e11  # 1000 s of pulses at 100 MHz (~17 min)
PAPER_RANGE_KM = 137.2  # ISPACS paper estimate (infinite decoys, Alice loss included, no dead time)

# --- palette (dataviz reference, light mode) ---
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
BOUND = "#898781"

plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 11,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "legend.frameon": False, "legend.labelcolor": INK2,
})


def compute():
    dev = q.load_device()
    L = np.arange(0, 201, 2.0)
    res = q.sweep(L, dev, q.SPAD, N_list=[1e9, N_MAIN, 1e13])
    key = f"fk_{N_MAIN:.0e}"
    res["chip_qber"] = dev.e_det[0]  # key (Z) basis
    sn = q.sweep(np.arange(0, 301, 2.0), dev, q.SNSPD, N_list=[N_MAIN])
    res["snspd_L"], res["snspd_fk"] = sn["L"], sn[key]
    CACHE.parent.mkdir(exist_ok=True)
    np.savez(CACHE, **res)
    return dict(np.load(CACHE))


def skr_text(r):
    """Key rate as a short label: 118 kbps, 7.2 kbps, 990 bps."""
    if r >= 1e4:
        return f"{r / 1e3:.0f} kbps"
    if r >= 1e3:
        return f"{r / 1e3:.1f} kbps"
    return f"{r:.0f} bps"


def end_label(ax, L, r, color, text, left=False, floor=1.0):
    """Dot at the last point still on the plot (rate >= floor) + direct label."""
    i = np.nonzero(r >= floor)[0][-1]
    ax.plot(L[i], r[i], "o", ms=6, color=color, mec=SURFACE, mew=1.5, zorder=5)
    ax.annotate(text, (L[i], r[i]), xytext=(-8 if left else 6, 0), textcoords="offset points",
                ha="right" if left else "left", va="center", color=INK2, fontsize=10)


def fig_main(d):
    L = d["L"]
    key = f"fk_{N_MAIN:.0e}"
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"wspace": 0.28})

    # (a) key rate
    series = [("inf_decoy", BOUND, "--", "Ideal decoy (asymptotic bound)"),
              (key, BLUE, "-", "Decoy BB84, finite key (17 min block)"),
              ("no_decoy", ORANGE, "-", "No decoy (current chip)")]
    for k, c, ls, lab in series:
        r = np.where(d[k] > 0, d[k], np.nan)
        a.semilogy(L, r, ls, color=c, lw=2, label=lab)
        end_label(a, L, d[k], c, f"{q.max_range(L, d[k]):.0f} km", left=(k == key))

    # SKR values at selected distances (decoy below its curve, no-decoy to its left)
    labels = [(key, BLUE, 0, (5, 17), "left"), (key, BLUE, 50, (4, -14), "left"),
              (key, BLUE, 100, (4, -14), "left"),
              ("no_decoy", ORANGE, 0, (8, -17), "left"), ("no_decoy", ORANGE, 50, (-6, -12), "right")]
    for k, c, Lp, offset, ha in labels:
        r = d[k][int(np.argmin(abs(L - Lp)))]
        if r <= 0:
            continue
        a.plot(Lp, r, "o", ms=5.5, color=c, mec=SURFACE, mew=1.2, zorder=5)
        a.annotate(skr_text(r), (Lp, r), xytext=offset, textcoords="offset points",
                   ha=ha, va="center", color=INK2, fontsize=9.5)
    a.set_ylim(1, 1e6)
    a.set_xlim(0, 200)
    a.set_xlabel("Fiber distance (km)")
    a.set_ylabel("Secret key rate (bit/s)")
    a.set_title("(a) Secret key rate", loc="left", fontsize=12)
    a.legend(loc="upper right", fontsize=9.5)

    # (b) QBER, drawn only where the key rate is positive; same colors as (a), labelled at the end
    for k, c, ls, lab in series:
        e = np.where(d[k] > 0, d[k + "_qber"] * 100, np.nan)
        b.plot(L, e, ls, color=c, lw=2, label=lab)
        i = np.nonzero(d[k] > 0)[0][-1]
        b.plot(L[i], e[i], "o", ms=6, color=c, mec=SURFACE, mew=1.5, zorder=5)
        b.annotate(f"{e[i]:.1f} %", (L[i], e[i]), xytext=(0, 9), textcoords="offset points",
                   ha="center", color=INK2, fontsize=10)
    b.axhline(11, color=INK2, lw=1, ls=":")
    b.text(198, 11.2, "BB84 threshold (11 %)", ha="right", va="bottom", color=INK2, fontsize=9.5)
    b.text(4, 9.6, "Colors as in (a). Curves end where the key rate reaches zero.\n"
           f"Chip polarization error alone: {100 * d['chip_qber']:.3f} % (circuit simulation)",
           va="top", color=INK2, fontsize=9.5, linespacing=1.5)
    b.set_xlim(0, 200)
    b.set_ylim(0, 13)
    b.set_xlabel("Fiber distance (km)")
    b.set_ylabel("QBER (%)")
    b.set_title("(b) Quantum bit error rate", loc="left", fontsize=12)

    fig.text(0.01, -0.04, "InGaAs SPAD: η = 20 %, dark 10⁻⁶/gate, afterpulse 1 %, dead time 10 µs · 100 MHz · "
             "0.2 dB/km · f_EC = 1.16 · ε_sec = 10⁻¹⁰ · Bob chip loss 4.5 dB (Alice's loss is pre-attenuation)",
             fontsize=8.5, color=MUTED)
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "fig_qkd_main.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_range(d):
    L = d["L"]
    rows = [("No decoy (current chip)", q.max_range(L, d["no_decoy"]), ORANGE),
            ("Decoy, 10 s block", q.max_range(L, d["fk_1e+09"]), BLUE),
            ("Decoy, 17 min block", q.max_range(L, d[f"fk_{N_MAIN:.0e}"]), BLUE),
            ("Decoy, ideal (bound)", q.max_range(L, d["inf_decoy"]), BOUND),
            ("Decoy, 17 min + SNSPD", q.max_range(d["snspd_L"], d["snspd_fk"]), BLUE)]
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    ax.set_axisbelow(True)
    y = np.arange(len(rows))[::-1]
    ax.axvline(PAPER_RANGE_KM, color=INK2, lw=1, ls=":", zorder=1)
    for yi, (lab, v, c) in zip(y, rows):
        ax.barh(yi, v, height=0.55, color=c, zorder=2)
        ax.text(v + 4, yi, f"{v:.0f} km", va="center", color=INK, fontsize=10.5, zorder=3,
                bbox=dict(boxstyle="square,pad=0.15", fc=SURFACE, ec="none"))
    ax.text(PAPER_RANGE_KM, len(rows) - 0.35, f"paper estimate ({PAPER_RANGE_KM:.0f} km)",
            ha="center", va="bottom", color=INK2, fontsize=9)
    ax.set_yticks(y, [r[0] for r in rows], color=INK2)
    ax.set_xlim(0, 260)
    ax.set_xlabel("Maximum secure distance (km)")
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(FIG / "fig_qkd_range.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    data = compute() if ("--rerun" in sys.argv or not CACHE.exists()) else dict(np.load(CACHE))
    fig_main(data)
    fig_range(data)
    print("saved figures to", FIG)
