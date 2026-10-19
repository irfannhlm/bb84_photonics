"""Alice -> Bob time-domain functional verification figure.

    (top)    first N_SHOW bit slots: Alice basis/bit drive, Bob basis drive, Bob's two optical
             outputs, and the decoded result per slot (sifted bit + check, or discarded).
    (bottom) per-bit extinction ratio for every matched-basis slot + run summary.

Data: waveforms from alice_bob.ich (INTERCONNECT result buffer, HDF5, read with h5py - no
Lumerical needed); per-slot results from BB84_transcript_<rate>Mbps.csv (extract_bit.lsf) when it
exists, otherwise recomputed here. Overlapping slots are cross-checked. The .ich can hold a partial
buffer (optical monitors shorter than the drives), so every signal keeps its own time axis.
Signal map (checked against the Transkrip CSVs): NRZ_1 = Alice basis, NRZ_2 = Alice bit,
NRZ_3 = Bob basis; Bob::opt_1 = detector 0, Bob::opt_2 = detector 1.

Detection uses ONE gate per slot, centred on the main pulse for both detectors. Gating each
detector on its own peak (as extract_bit.lsf originally did) makes the dark detector integrate
the ~5.4 ns echo of the S-parameter filters instead of the empty slot, which drops the ER to
10-20 dB.

Output: fig_alice_bob_timing.png next to this file.
"""
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ICH = HERE / "alice_bob.ich"
OUT = HERE / "fig_alice_bob_timing.png"
GATE = 1e-9  # full gate width (s), SPAD-like
N_SHOW = 12  # slots in the waveform strip
V_TH = 0.6  # logic threshold on the 0 / 1.2 V drives

# --- palette (dataviz reference, light mode; same as qkd_model/make_figures.py) ---
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GOOD = "#008300"
plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 10,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "legend.frameon": False, "legend.labelcolor": INK2,
})

# ============================================================
# READ THE RESULT BUFFER
# ============================================================
f = h5py.File(ICH, "r")
ROOT = "::Root Element::"


def electrical(name):
    d = f[ROOT + name + "::output::monitor/data"][:]
    return d["time"], d["value"]


def optical_power(name):
    """Total power summed over the modes of an optical monitor."""
    d = f[ROOT + name + "::monitor/data"]
    t = d["time"][:]
    p = np.array([sum(m["value"][0]["real"] ** 2 + m["value"][0]["imag"] ** 2 for m in s)
                  for s in d["signal"][:]])
    return t, p


drives = {k: electrical(n) for k, n in
          (("Alice basis", "NRZ_1"), ("Alice bit", "NRZ_2"), ("Bob basis", "NRZ_3"))}
t, p0 = optical_power("Bob::opt_1")
t1, p1 = optical_power("Bob::opt_2")
p1 = np.interp(t, t1, p1)  # monitors can differ by a few samples
T = 20e-9  # 50 Mbps (Root Element "bitrate")
CSV = HERE / f"BB84_transcript_{1e-6 / T:.0f}Mbps.csv"
dt = t[1] - t[0]

# ============================================================
# BIT DECISIONS
# ============================================================
n_bits = int(min(tt[-1] for tt, _ in drives.values()) / T)
mid = (np.arange(n_bits) + 0.5) * T
basis_a, bit_a, basis_b = ((np.interp(mid, *drives[k]) > V_TH).astype(int)
                           for k in ("Alice basis", "Alice bit", "Bob basis"))

total = p0 + p1
t_first = t[np.argmax(np.where(t < T, total, 0))]  # first pulse at Bob = timing reference
centres = t_first + np.arange(n_bits) * T
keep = centres + GATE / 2 <= t[-1]
centres, basis_a, bit_a, basis_b = centres[keep], basis_a[keep], bit_a[keep], basis_b[keep]


def gated(p, c):
    s = np.abs(t - c) <= GATE / 2
    return p[s].sum() * dt


e0 = np.array([gated(p0, c) for c in centres])
e1 = np.array([gated(p1, c) for c in centres])
bit_b = (e1 > e0).astype(int)
contrast = 10 * np.log10(np.maximum(e0, e1) / np.minimum(e0, e1))
matched = basis_a == basis_b
errors = matched & (bit_b != bit_a)
STATE = {(0, 0): "H", (0, 1): "V", (1, 0): "−45°", (1, 1): "+45°"}
strip = dict(basis_a=basis_a, bit_a=bit_a, basis_b=basis_b, bit_b=bit_b, matched=matched)

if CSV.exists():  # authoritative per-slot results from INTERCONNECT (full run)
    tr = pd.read_csv(CSV)
    k = len(centres)
    same = ((tr.alice_basis[:k] == basis_a) & (tr.alice_bit[:k] == bit_a) & (tr.bob_basis[:k] == basis_b)
            & (tr.bob_bit[:k] == bit_b)).all()
    dev = np.abs(tr.contrast_dB[:k].to_numpy() - contrast)[matched].max()
    print(f"cross-check vs {CSV.name}: first {k} slots identical bits = {same}, "
          f"max matched ER difference {dev:.2f} dB")
    basis_a, bit_a, basis_b = tr.alice_basis.to_numpy(), tr.alice_bit.to_numpy(), tr.bob_basis.to_numpy()
    bit_b, contrast = tr.bob_bit.to_numpy(), tr.contrast_dB.to_numpy()
    matched = tr.basis_match.eq("MATCHED").to_numpy()
    errors = matched & (bit_b != bit_a)

n = len(contrast)
summary = (f"{n} slots · {matched.sum()} matched basis ({matched.mean():.0%}) · "
           f"{errors.sum()} bit errors · matched ER {contrast[matched].min():.1f}–"
           f"{contrast[matched].max():.1f} dB · mismatched split within {contrast[~matched].max():.2f} dB")
print(summary)

# ============================================================
# FIGURE
# ============================================================
fig = plt.figure(figsize=(12.5, 9.2))
gs = fig.add_gridspec(7, 1, height_ratios=[0.55, 0.55, 0.55, 1.6, 0.55, 0.55, 1.6], hspace=0.18,
                      left=0.11, right=0.985, top=0.885, bottom=0.06)
fig.suptitle("BB84 Alice → Bob functional verification (INTERCONNECT, time domain, 50 Mbps)",
             x=0.11, ha="left", fontsize=14, color=INK, fontweight="semibold")
fig.text(0.11, 0.915, "Classical Gaussian pulses (0.38 ns FWHM) through both chips; one 1 ns gate per "
         "slot on Bob's two outputs. Thermal heaters modelled as ideal fast shifters.",
         ha="left", color=INK2, fontsize=10)

x_end = N_SHOW * T * 1e9
shift = (t_first - 0.5 * T)  # Bob output lags the drive; align slot k of both on the same axis
ax_d = [fig.add_subplot(gs[i]) for i in range(3)]
for ax, (lab, (tt, v)) in zip(ax_d, drives.items()):
    ax.plot(tt * 1e9, v, color=INK2, lw=1.5)
    ax.set_ylim(-0.3, 1.6)
    ax.set_yticks([0, 1.2], ["0", "1"])
    ax.set_ylabel(lab, rotation=0, ha="right", va="center", color=INK2)
    ax.grid(axis="y", visible=False)

ax_o = fig.add_subplot(gs[3], sharex=ax_d[0])
pn = 1 / total.max()
ax_o.plot((t - shift) * 1e9, p0 * pn, color=BLUE, lw=1.5, label="detector 0")
ax_o.plot((t - shift) * 1e9, p1 * pn, color=ORANGE, lw=1.5, label="detector 1")
ax_o.set_ylabel("Bob output\n(norm. power)", rotation=0, ha="right", va="center", color=INK2)
ax_o.set_ylim(-0.05, 1.35)
ax_o.legend(loc="upper left", ncol=2, fontsize=9)

ax_r = fig.add_subplot(gs[4], sharex=ax_d[0])
ax_r.set_ylim(0, 1)
ax_r.set_yticks([])
ax_r.grid(False)
ax_r.spines["left"].set_visible(False)
ax_r.set_ylabel("Result", rotation=0, ha="right", va="center", color=INK2)
for k in range(N_SHOW):
    xc = (k + 0.5) * T * 1e9
    ax_r.text(xc, 0.78, f"{STATE[(strip["basis_a"][k], strip["bit_a"][k])]} → {'+×'[strip["basis_b"][k]]}",
              ha="center", va="center", fontsize=9, color=INK2)
    if strip["matched"][k]:
        ok = strip["bit_b"][k] == strip["bit_a"][k]
        ax_r.text(xc, 0.3, f"bit {strip["bit_b"][k]} {'✓' if ok else '✗'}", ha="center", va="center",
                  fontsize=10, color=GOOD if ok else INK, fontweight="semibold")
    else:
        ax_r.text(xc, 0.3, "discarded", ha="center", va="center", fontsize=9, color=MUTED)

for ax in ax_d + [ax_o, ax_r]:
    for k in range(N_SHOW + 1):
        ax.axvline(k * T * 1e9, color=GRID, lw=0.8, zorder=0)
    ax.grid(axis="x", visible=False)
    ax.set_xlim(0, x_end)
    ax.tick_params(labelbottom=False, bottom=False)
ax_r.set_xticks((np.arange(N_SHOW) + 0.5) * T * 1e9, [str(k + 1) for k in range(N_SHOW)])
ax_r.tick_params(labelbottom=True)
ax_r.set_xlabel("Bit slot (20 ns each; Bob outputs shifted by the Alice → Bob delay)")

# --- per-bit ER, all matched slots ---
ax_e = fig.add_subplot(gs[6])
idx = np.nonzero(matched)[0]
for (b, v), col in zip(sorted(STATE), (BLUE, ORANGE, AQUA, YELLOW)):
    s = idx[(basis_a[idx] == b) & (bit_a[idx] == v)]
    ax_e.plot(s + 1, contrast[s], "o", ms=6, color=col, mec=SURFACE, mew=1, label=STATE[(b, v)])
ax_e.set_xlim(0, n + 1)
ax_e.set_ylim(25, 55)
ax_e.set_ylabel("Matched-basis\nER (dB)", rotation=0, ha="right", va="center", color=INK2)
ax_e.set_xlabel("Bit slot")
ax_e.legend(loc="upper right", ncol=4, fontsize=9, handletextpad=0.2, columnspacing=1.2)
ax_e.set_title(summary, loc="left", fontsize=10.5, pad=6)

fig.savefig(OUT, dpi=200)
print(f"saved {OUT}")
