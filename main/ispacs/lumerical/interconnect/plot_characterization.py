"""Alice / Bob frequency-domain characterization summary figure.

    (top)    Alice transmission vs wavelength for the four BB84 states, projected on the
             intended polarization and its orthogonal one (gap = polarization extinction).
    (bottom) metrics at 1550 nm: Alice IL / ER / azimuth error, Bob IL / ER per input.

Data: bob_test/interconnect_reference.npz (saved sweep results exported by
bob_test/export_reference.py from alice_test_freq.icp::sweep_basis and
bob_test_freq.icp::sweep_base). No Lumerical needed.

Output: fig_alice_bob_characterization.png next to this file.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REF = np.load(HERE / "bob_test" / "interconnect_reference.npz")
OUT = HERE / "fig_alice_bob_characterization.png"
WL0 = 1550e-9

# --- palette (dataviz reference, light mode; same as qkd_model/make_figures.py) ---
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 10,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "legend.frameon": False, "legend.labelcolor": INK2,
})

# ============================================================
# ALICE: complex fiber field per state -> intended / orthogonal power
# ============================================================
wl = REF["alice/gain_X/wavelength"].ravel()
i0 = np.argmin(abs(wl - WL0))
gX, gY = REF["alice/gain_X/X gain (dB)"], REF["alice/gain_Y/Y gain (dB)"]  # (wl, bit, basis)
aX, aY = REF["alice/angle_X/X angle (rad)"], REF["alice/angle_Y/Y angle (rad)"]
Ex = np.sqrt(10 ** (gX / 10)) * np.exp(1j * aX)
Ey = np.sqrt(10 ** (gY / 10)) * np.exp(1j * aY)

# (label, basis idx, bit idx, target azimuth deg) -- assignment from the 1550 nm azimuths
ALICE_STATES = [("H", 0, 0, 0.0), ("V", 0, 1, 90.0), ("-45°", 1, 0, -45.0), ("+45°", 1, 1, 45.0)]


def project(ex, ey, az_deg):
    """Power in the linear polarization at az_deg and in its orthogonal (linear, dB)."""
    t = np.radians(az_deg)
    co = np.abs(np.cos(t) * ex + np.sin(t) * ey) ** 2
    cross = np.abs(-np.sin(t) * ex + np.cos(t) * ey) ** 2
    return co, cross


def azimuth_ellipticity(ex, ey):
    s0 = abs(ex) ** 2 + abs(ey) ** 2
    s1 = abs(ex) ** 2 - abs(ey) ** 2
    s2 = 2 * np.real(ex * np.conj(ey))
    s3 = -2 * np.imag(ex * np.conj(ey))
    return np.degrees(0.5 * np.arctan2(s2, s1)), np.degrees(0.5 * np.arcsin(s3 / s0))


alice_rows = []
alice_curves = []
for name, bas, bit, target in ALICE_STATES:
    ex, ey = Ex[:, bit, bas], Ey[:, bit, bas]
    co, cross = project(ex, ey, target)
    alice_curves.append((name, bas, bit, 10 * np.log10(co), 10 * np.log10(cross)))
    az, ell = azimuth_ellipticity(ex[i0], ey[i0])
    il = -10 * np.log10(co[i0] + cross[i0])
    er = 10 * np.log10(co[i0] / cross[i0])
    alice_rows.append([name, f"{bas} / {bit}", f"{il:.2f}", f"{er:.1f}",
                       f"{az:+.2f}", f"{az - target:+.2f}", f"{ell:+.2f}"])

# ============================================================
# BOB: two detector outputs vs input polarization (rotator) and basis
# ============================================================
rot = np.degrees(REF["bob/gain1/rotation"].ravel())
g1, g2 = REF["bob/gain1/TE gain (dB)"], REF["bob/gain2/TE gain (dB)"]  # (wl, rotation, base)
bob_rows = []
for base, inputs in ((0, (0.0, 90.0)), (1, (-45.0, 45.0))):
    for r in inputs:
        k = int(np.argmin(abs(rot - r)))
        a, b = g1[i0, k, base], g2[i0, k, base]
        il = -10 * np.log10(10 ** (a / 10) + 10 ** (b / 10))
        bob_rows.append([f"{r:+.0f}°" if r else "0°", str(base), f"{a:.2f}", f"{b:.2f}",
                         f"{il:.2f}", f"{abs(a - b):.1f}"])
# mismatched basis: should split 50/50 -> report the worst imbalance
mism = [abs(g1[i0, k, base] - g2[i0, k, base])
        for base, inputs in ((0, (-45.0, 45.0)), (1, (0.0, 90.0)))
        for k in [int(np.argmin(abs(rot - r))) for r in inputs]]

# ============================================================
# FIGURE
# ============================================================
fig = plt.figure(figsize=(12.5, 8.2))
gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1], hspace=0.45, wspace=0.08,
                      left=0.06, right=0.985, top=0.83, bottom=0.03)
fig.suptitle("BB84 chip characterization (INTERCONNECT, frequency domain)",
             x=0.06, ha="left", fontsize=14, color=INK, fontweight="semibold")
fig.text(0.06, 0.915, "Alice: fiber-output transmission per state, projected on the intended polarization "
         "and its orthogonal; the gap is the polarization extinction ratio",
         ha="left", color=INK2, fontsize=10)

wl_nm = wl * 1e9
axes = []
for j, (name, bas, bit, co_db, cross_db) in enumerate(alice_curves):
    ax = fig.add_subplot(gs[0, j], sharey=axes[0] if axes else None)
    axes.append(ax)
    ax.plot(wl_nm, co_db, color=BLUE, lw=2, label="intended pol.")
    ax.plot(wl_nm, cross_db, color=ORANGE, lw=1.2, label="orthogonal pol.")
    ax.axvline(1550, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=0)
    ax.plot(1550, co_db[i0], "o", ms=6, color=BLUE, mec=SURFACE, mew=1.5, zorder=5)
    ax.plot(1550, cross_db[i0], "o", ms=6, color=ORANGE, mec=SURFACE, mew=1.5, zorder=5)
    box = dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.9)
    ax.annotate(f"{co_db[i0]:.2f} dB", (1550, co_db[i0]), xytext=(6, 6),
                textcoords="offset points", color=INK2, fontsize=9, bbox=box)
    ax.annotate(f"ER {co_db[i0] - cross_db[i0]:.1f} dB", (1550, cross_db[i0]), xytext=(0, 34),
                textcoords="offset points", ha="center", color=INK, fontsize=9, bbox=box, zorder=6)
    ax.set_title(f"{name}   (basis {bas}, bit {bit})", loc="left", fontsize=11)
    ax.set_xlim(wl_nm.min(), wl_nm.max())
    ax.set_xticks([1535, 1550, 1565])
    ax.set_xlabel("Wavelength (nm)")
    if j == 0:
        ax.set_ylabel("Transmission (dB)")
    else:
        ax.tick_params(labelleft=False)
axes[0].set_ylim(-75, 0)
axes[0].legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=2, fontsize=10, handlelength=1.5,
               borderaxespad=0)


def table(ax, title, header, rows, note, worst):
    """worst: {column index: row index} cells shown bold (worst case per metric)."""
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=12, pad=4)
    t = ax.table(cellText=rows, colLabels=header, loc="upper left", cellLoc="center",
                 bbox=[0, 0.16, 1, 0.80])
    t.auto_set_font_size(False)
    t.set_fontsize(10)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.set_linewidth(0.6)
        cell.visible_edges = "B"
        cell.set_facecolor(SURFACE)
        cell.get_text().set_color(INK if r else INK2)
        if r == 0:
            cell.get_text().set_fontweight("semibold")
            cell.set_edgecolor(AXIS)
        elif worst.get(c) == r - 1:
            cell.get_text().set_fontweight("bold")
            cell.get_text().set_fontsize(11)
    ax.text(0, 0.07, note, transform=ax.transAxes, ha="left", va="top", color=INK2, fontsize=9.5)


col = lambda rows, j: np.array([float(r[j]) for r in rows])  # noqa: E731
table(fig.add_subplot(gs[1, :2]), "Alice (transmitter) @ 1550 nm",
      ["State", "Basis / bit", "IL (dB)", "ER (dB)", "Azimuth (°)", "Angle err. (°)", "Ellipticity (°)"],
      alice_rows,
      "Bold = worst case.  ER = intended / orthogonal power, so it already includes the\n"
      "angle error and ellipticity.",
      {2: col(alice_rows, 2).argmax(), 3: col(alice_rows, 3).argmin(),
       5: abs(col(alice_rows, 5)).argmax()})

table(fig.add_subplot(gs[1, 2:]), "Bob (receiver) @ 1550 nm, matched basis",
      ["Input (rotator)", "Basis", "Det. 1 (dB)", "Det. 2 (dB)", "IL (dB)", "ER (dB)"],
      bob_rows,
      f"Bold = worst case.  Mismatched basis splits 50/50 within {max(mism):.2f} dB.\n"
      "Det. 1 = gain1 (MMI_1.opt_2), Det. 2 = gain2 (MMI_1.opt_1).",
      {4: col(bob_rows, 4).argmax(), 5: col(bob_rows, 5).argmin()})

fig.savefig(OUT, dpi=200)
print(f"saved {OUT}")
for r in alice_rows:
    print("Alice", r)
for r in bob_rows:
    print("Bob  ", r)
