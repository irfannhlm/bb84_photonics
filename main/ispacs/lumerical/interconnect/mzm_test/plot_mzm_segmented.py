"""Segmented MZM slide figure: truth table of heater sections -> phase -> polarization.

The 3pi/2 heater (L = 300 um) is split into L/3 (basis, pi/2) and 2L/3 (bit, pi), both driven
NRZ 0 V / 1.2 V. Each row draws the top arm to scale with the driven sections filled.
Ideal phases are shown. The polarization per (basis, bit) comes from the Alice circuit model
(bob_test/bob_model.py, Stokes S1/S2): 00 -> H, 01 -> V, 10 -> A, 11 -> D.
Output: fig_mzm_segmented.png next to this script.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

from plot_mzm_sweep import AXIS, GRID, HERE, INK, INK2, MUTED, ORANGE, SURFACE

V_HIGH = 1.2
# (basis, bit, delta phi, MZM output, polarization, arrow angle in deg)
ROWS = [(0, 0, "0", "Through", "H", 0), (0, 1, "π", "Drop", "V", 90),
        (1, 0, "π/2", "50:50", "A", -45), (1, 1, "3π/2", "50:50", "D", 45)]

# column x positions (data units); heater bar: L/3 then 2L/3, drawn to scale
X_BASIS, X_BIT = 0.5, 1.35
X_H0, UNIT = 2.2, 1.1  # bar start, width of L/3
X_DPHI, X_OUT, X_POL = 6.35, 7.75, 9.35
ROW_H, BAR_H = 1.0, 0.52


def heater_bar(ax, y, basis_on, bit_on):
    ax.plot([X_H0 - 0.25, X_H0 + 3 * UNIT + 0.25], [y, y], color=MUTED, lw=2.2, zorder=1, solid_capstyle="butt")
    for x0, w, on in ((X_H0, UNIT, basis_on), (X_H0 + UNIT, 2 * UNIT, bit_on)):
        ax.add_patch(Rectangle((x0 + 0.03, y - BAR_H / 2), w - 0.06, BAR_H, lw=1.4,
                               ec=ORANGE if on else INK2, fc=ORANGE if on else SURFACE, zorder=2))
        ax.text(x0 + w / 2, y, f"{V_HIGH} V" if on else "0 V", ha="center", va="center", fontsize=10.5,
                color=SURFACE if on else MUTED, fontweight="bold" if on else "normal", zorder=3)


def pol_arrow(ax, x, y, angle):
    r = 0.3
    dx, dy = r * np.cos(np.radians(angle)), r * np.sin(np.radians(angle))
    ax.add_patch(FancyArrowPatch((x - dx, y - dy), (x + dx, y + dy), arrowstyle="<|-|>", mutation_scale=11,
                                 lw=1.8, color=INK, zorder=3))


def main():
    fig, ax = plt.subplots(figsize=(10.2, 4.3))
    ax.set_xlim(0, 10.2)
    ax.set_ylim(-0.7, 5.35)
    ax.set_aspect("equal")
    ax.axis("off")

    # headers
    yh = 4.7
    hdr = dict(ha="center", va="center", fontsize=11, color=INK2)
    ax.text(X_BASIS, yh, "Basis", **hdr)
    ax.text(X_BIT, yh, "Bit", **hdr)
    ax.text(X_H0 + UNIT / 2, yh + 0.1, "basis heater\nL/3  →  π/2", **{**hdr, "fontsize": 10.5})
    ax.text(X_H0 + 2 * UNIT, yh + 0.1, "bit heater\n2L/3  →  π", **{**hdr, "fontsize": 10.5})
    ax.text(X_DPHI, yh, "Δφ", **hdr)
    ax.text(X_OUT, yh, "MZM output", **hdr)
    ax.text(X_POL, yh, "Polarization", **hdr)
    ax.plot([0.1, 10.1], [yh - 0.45, yh - 0.45], color=AXIS, lw=1)

    for i, (basis, bit, dphi, out, pol, ang) in enumerate(ROWS):
        y = 3.6 - i * ROW_H - (0.2 if i >= 2 else 0)  # small gap between the two bases
        cell = dict(ha="center", va="center", fontsize=12, color=INK)
        ax.text(X_BASIS, y, str(basis), **cell)
        ax.text(X_BIT, y, str(bit), **cell)
        heater_bar(ax, y, basis, bit)
        ax.text(X_DPHI, y, dphi, **cell)
        ax.text(X_OUT, y, out, **{**cell, "color": INK2, "fontsize": 11})
        ax.text(X_POL - 0.3, y, pol, **{**cell, "fontweight": "bold"})
        pol_arrow(ax, X_POL + 0.3, y, ang)
    ax.plot([0.1, 10.1], [2.0, 2.0], color=GRID, lw=1)  # between the two bases

    ax.text(X_H0 + 1.5 * UNIT, -0.35, f"both sections driven NRZ: 0 V / {V_HIGH} V", ha="center",
            va="center", fontsize=10, color=MUTED, style="italic")

    out = HERE / "fig_mzm_segmented.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
