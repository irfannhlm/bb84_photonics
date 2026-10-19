"""MZM transmission vs heater voltage, with the four BB84 operating points marked.

Model: the mzm_test circuit rebuilt with bob_test/bob_model.py parts --
    MMI_1 -> wg_heater_1 (300 um, swept)  -> MMI_2      input at MMI_1.opt_2
          -> wg_heater_2 (300 um, 0 V)    ->            outputs MMI_2.opt_3 / opt_4
using the calibrated heater (phase ~ V^2 L, loss ~ L) and the MMI S-parameter file.

Check: the saved INTERCONNECT sweep in mzm_test_sweep/ is read directly (no Lumerical needed):
    sweep_N.ich  HDF5 impulse response at MMI_2 opt_3 / opt_4 -> port energy = sum |E|^2
    sweep_N.icp  DC_1 "amplitude" (big-endian double after the property name) -> voltage
It is scaled to the model's total at 0 V (the impulse energy has no absolute reference).
Output: fig_mzm_voltage_sweep.png next to this script.
"""
import struct
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bob_test"))
import bob_model as bm  # noqa: E402

SWEEP = HERE / "mzm_test_sweep"
MON = "::Root Element::MMI_2::{}::monitor/data"
L_HEATER = 300e-6

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 11,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "legend.frameon": False, "legend.labelcolor": INK2,
})


# ============================================================
# MODEL
# ============================================================
def mzm_model(f0, k_phase, loss):
    """Return T(V) -> (through, drop) power transmission for the bare MZM."""
    conn = [("MMI1.opt_4", "hA.opt_1"), ("hA.opt_2", "MMI2.opt_2"),
            ("MMI1.opt_3", "hB.opt_1"), ("hB.opt_2", "MMI2.opt_1")]
    ext = ["MMI1.opt_1", "MMI1.opt_2", "MMI2.opt_3", "MMI2.opt_4"]

    def T(V):
        comps = {"MMI1": bm.mmi(f0), "MMI2": bm.mmi(f0),
                 "hA": bm.heater(L_HEATER, V, k_phase, loss), "hB": bm.heater(L_HEATER, 0.0, k_phase, loss)}
        S = bm.solve_network(comps, conn, ext)
        return abs(S[2, 1]) ** 2, abs(S[3, 1]) ** 2

    return T


def operating_points(T, k_phase):
    """0 V, and the voltages where the through fraction hits 1/2, 0, 1/2 (dphi = pi/2, pi, 3pi/2)."""
    v_nom = lambda phi: np.sqrt(phi / (k_phase * L_HEATER))  # noqa: E731
    frac = lambda V: T(V)[0] / sum(T(V))  # noqa: E731
    v_half1 = brentq(lambda V: frac(V) - 0.5, v_nom(0.25 * np.pi), v_nom(0.75 * np.pi))
    v_half2 = brentq(lambda V: frac(V) - 0.5, v_nom(1.25 * np.pi), v_nom(1.75 * np.pi))
    vv = np.linspace(v_nom(0.8 * np.pi), v_nom(1.2 * np.pi), 4001)
    v_pi = vv[np.argmin([frac(V) for V in vv])]
    return [(0.0, "Fully through", "0"), (v_half1, "50:50", "π/2"),
            (v_pi, "Fully drop", "π"), (v_half2, "50:50", "3π/2")]


# ============================================================
# INTERCONNECT REFERENCE
# ============================================================
def dc1_voltage(icp):
    b = icp.read_bytes()
    p = b.find("DC_1".encode("utf-16le"))
    q = b.find("amplitude".encode("utf-16le"), p) + len("amplitude".encode("utf-16le"))
    return struct.unpack(">d", b[q + 4:q + 12])[0]


def port_energy(f, port):
    return sum(abs(complex(r["signal"][0]["value"]["real"].item(),
                           r["signal"][0]["value"]["imag"].item())) ** 2 for r in f[MON.format(port)][:])


def load_interconnect():
    V, E3, E4 = [], [], []
    for n in range(1, len(list(SWEEP.glob("sweep_*.ich"))) + 1):
        V.append(dc1_voltage(SWEEP / f"sweep_{n}.icp"))
        with h5py.File(SWEEP / f"sweep_{n}.ich", "r") as f:
            E3.append(port_energy(f, "opt_3"))
            E4.append(port_energy(f, "opt_4"))
    V, E3, E4 = map(np.array, (V, E3, E4))
    # the output stops changing above ~1.2 V in INTERCONNECT; keep only the live points
    stuck = np.nonzero(np.abs(np.diff(E3)) < 1e-9 * E3.max())[0]
    n = stuck[0] + 1 if len(stuck) else len(V)
    return V[:n], E3[:n], E4[:n]


# ============================================================
# FIGURE
# ============================================================
def db(x):
    return 10 * np.log10(np.maximum(x, 1e-7))


def main():
    f0, k_phase, loss, _ = bm.calibrate()
    T = mzm_model(f0, k_phase, loss)
    pts = operating_points(T, k_phase)

    Vf = np.linspace(0, 1.5, 751)
    Tf = np.array([T(V) for V in Vf])
    Vi, E3, E4 = load_interconnect()
    scale = sum(T(0.0)) / (E3[0] + E4[0])
    Ti = np.c_[E3, E4] * scale

    # top axis: phase difference, ticked where each phase occurs -> evenly spaced phases sit at
    # unevenly spaced voltages (phase ~ V^2)
    v_of = lambda phi_pi: np.sqrt(phi_pi * np.pi / (k_phase * L_HEATER))  # noqa: E731

    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ports = ((BLUE, "Through"), (ORANGE, "Drop"))
    for v, _, _ in pts:
        ax.axvline(v, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)

    for j, (c, lab) in enumerate(ports):
        ax.plot(Vf, Tf[:, j], color=c, lw=2, zorder=3)
        ax.plot(Vi, Ti[:, j], "x", ms=4.5, mew=1.1, color=c, alpha=0.7, zorder=4)
        ax.annotate(lab, (Vf[-1], Tf[-1, j]), xytext=(6, 0), textcoords="offset points",
                    va="center", color=INK2, fontsize=10.5, annotation_clip=False)
    for v, _, _ in pts:
        for t, (c, _) in zip(T(v), ports):
            ax.plot(v, t, "o", ms=8, color=c, mec=SURFACE, mew=1.8, zorder=5)
    ax.plot([], [], "x", ms=4.5, mew=1.1, color=MUTED, label="INTERCONNECT sweep")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.5), fontsize=9.5, handletextpad=0.3)

    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Transmission")
    ax.set_xlim(-0.03, 1.5)
    ax.set_xticks([v for v, _, _ in pts] + [1.5], [f"{v:.2f}" for v, _, _ in pts] + ["1.50"])
    ax.set_xlabel("Heater voltage (V)")

    top = ax.secondary_xaxis("top")
    top.set_xticks([v_of(p) for p in (0, 0.5, 1, 1.5, 2)], ["0", "π/2", "π", "3π/2", "2π"])
    top.set_xlabel("Phase difference Δφ", color=INK2, labelpad=6)
    top.tick_params(colors=INK2)
    top.spines["top"].set_visible(True)
    top.spines["top"].set_color(AXIS)

    out = HERE / "fig_mzm_voltage_sweep.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"heater: {k_phase * L_HEATER / np.pi:.4f} pi/V^2 over 300 um, loss {loss * L_HEATER:.3f} dB")
    for v, name, phase in pts:
        t3, t4 = T(v)
        print(f"{name:14s} {phase:5s} V = {v:.3f} V   through {db(t3):7.2f} dB  drop {db(t4):7.2f} dB"
              f"   split {t3 / (t3 + t4):.3f}")
    ti = np.array([T(V) for V in Vi])
    print(f"INTERCONNECT points used: {len(Vi)} (0 - {Vi.max():.3f} V); "
          f"max |model - IC| = {abs(ti - Ti).max():.4f} (linear)")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
