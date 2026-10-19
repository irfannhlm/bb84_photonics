"""Independent cross-check of the BB84 model with Perceval (Quandela).

The optical and photon-statistics part of qkd_model.py is recomputed with Perceval's
polarization + loss simulator:

  Alice chip   rails 0/1 = MMI_1.opt_1/opt_2 -> simulated Jones vector (exact Unitary)
  PSGC         rails (fiber X, fiber Y) -> H/V of ONE fiber mode (HWP + PBS)
  Fiber        LC (0.2 dB/km) and optional PR (residual polarization rotation)
  Bob chip     PBS + HWP back to rails -> simulated 2x2 map (Unitary + LC + Unitary, i.e. SVD)
               -> LC (detector efficiency)
  Detection    threshold detectors on MMI_1.opt_1 (bit 0) and MMI_1.opt_2 (bit 1)
  Source       phase-randomized weak coherent pulse = Poisson mixture of Fock states |n>

Perceval replaces the analytic optics (qkd_model.OPTICS); detector noise (dark counts,
afterpulsing, dead time) and the key-rate formulas are then applied unchanged, using the
same optimized protocol parameters. Perceval has no dark-count model, so that part is
shared by construction; everything optical and quantum comes from Perceval.

Outputs
  fig_perceval_check.png  SKR and QBER vs distance: model curves with Perceval points
  console                 deviations, plus textbook BB84 checks (mismatched bases -> 50 %,
                          intercept-resend attack -> 25 %, fiber rotation)
"""
from contextlib import contextmanager
from math import exp, factorial

import matplotlib.pyplot as plt
import numpy as np
import perceval as pcvl
from matplotlib.lines import Line2D
from perceval.components import HWP, LC, PBS, PR, PS, Unitary

import qkd_model as q
from make_figures import BLUE, BOUND, CACHE, FIG, INK, INK2, N_MAIN, ORANGE, SURFACE

bm = q.bob_model
F0, K_PHASE, LOSS, _ = bm.calibrate()
V = bm.V_HIGH
N_MAX = 8  # photon-number cutoff for the Poisson mixture
DET = q.SPAD
STATES = [(0, 0), (0, 1), (1, 0), (1, 1)]  # (basis, bit); bit 0 -> MMI_1.opt_1


# ==========================================
# CIRCUIT BUILDING BLOCKS
# ==========================================
def add_lossy(p, T):
    """Add a non-unitary 2x2 map T on rails (0, 1) as Unitary - LC - Unitary (SVD)."""
    U, s, Vh = np.linalg.svd(T)
    p.add(0, Unitary(pcvl.Matrix(Vh)))
    for m in (0, 1):
        p.add(m, LC(float(np.clip(1 - s[m] ** 2, 0, 1))))
    p.add(0, Unitary(pcvl.Matrix(U)))


def rails_to_fiber(p):
    """PSGC as an ideal combiner: rail 0 -> H, rail 1 -> V of fiber mode 1."""
    p.add(1, PS(-np.pi / 2))  # cancel the i from the HWP
    p.add(1, HWP(np.pi / 4))  # H -> iV
    p.add(0, PBS())  # H crosses to mode 1, V stays: both end up in mode 1


def fiber_to_rails(p):
    """Inverse: H of fiber mode 1 -> rail 0, V -> rail 1 (both H afterwards)."""
    p.add(0, PBS())
    p.add(1, HWP(np.pi / 4))
    p.add(1, PS(-np.pi / 2))


def add_alice(p, basis, bit):
    """Alice's output state on the rails (fiber X, Y) for a photon entering MMI_1.opt_1.

    mu is defined at Alice's fiber output, so the state is normalized; it is prepared with
    the exact 2x2 unitary whose first column is Alice's simulated Jones vector.
    """
    E = bm.alice_jones(bm.circuit(F0, V * bit, V * basis, K_PHASE, LOSS))
    a, b = E / np.linalg.norm(E)
    p.add(0, Unitary(pcvl.Matrix(np.array([[a, -np.conj(b)], [b, np.conj(a)]]))))


def bob_T(basis):
    return bm.circuit(F0, 0.0, V * basis, K_PHASE, LOSS)[:2, 2:]


def link(a_basis, bit, b_basis, L, rot=0.0):
    """Processor for Alice -> fiber -> Bob -> detectors."""
    p = pcvl.Processor("SLOS", 2)
    add_alice(p, a_basis, bit)
    rails_to_fiber(p)
    if L > 0:
        p.add(1, LC(1 - 10 ** (-q.ALPHA * L / 10)))
    if rot:
        p.add(1, PR(rot))
    fiber_to_rails(p)
    add_lossy(p, bob_T(b_basis))
    for m in (0, 1):
        p.add(m, LC(1 - DET.eta))
    p.min_detected_photons_filter(0)
    return p


def click_probs(p, n):
    """P(click on rail 0 only, rail 1 only, both) for n photons entering MMI_1.opt_1."""
    if n == 0:
        return np.zeros(3)
    p.with_input(pcvl.BasicState(f"|{n}{{P:H}},0>"))
    out = np.zeros(3)
    for s, pr in p.probs()["results"].items():
        c0, c1 = s[0] > 0, s[1] > 0
        if c0 and c1:
            out[2] += pr
        elif c0:
            out[0] += pr
        elif c1:
            out[1] += pr
    return out


class PercevalOptics:
    """Optical click probabilities from Perceval, cached per (distance, basis).

    table[n] = (correct only, wrong only, both) for an n-photon Fock state, averaged over
    the two bits; a weak coherent pulse is the Poisson mixture of these.
    """

    def __init__(self):
        self.cache = {}

    def table(self, L, basis):
        key = (float(L), basis)
        if key not in self.cache:
            t = np.zeros((N_MAX + 1, 3))
            for b, bit in STATES:
                if b != basis:
                    continue
                p = link(b, bit, b, L)
                for n in range(1, N_MAX + 1):
                    c = click_probs(p, n)  # (rail 0 only, rail 1 only, both)
                    t[n] += np.array([c[bit], c[1 - bit], c[2]]) / 2
            self.cache[key] = t
        return self.cache[key]

    def wcp(self, mu, L, dev, det, basis):
        w = np.array([exp(-mu) * mu**n / factorial(n) for n in range(N_MAX + 1)])
        return w @ self.table(L, basis)

    def single(self, L, dev, det, basis):
        return self.table(L, basis)[1]


@contextmanager
def use_optics(optics):
    saved = dict(q.OPTICS)
    q.OPTICS.update(wcp=optics.wcp, single=optics.single)
    try:
        yield
    finally:
        q.OPTICS.update(saved)


# ==========================================
# CHECKS
# ==========================================
def check_rates(dev, distances):
    """SKR and QBER for the three protocols, model vs Perceval, same protocol parameters.

    Returns {protocol: array[distance, (SKR model, SKR Perceval, QBER model, QBER Perceval)]}.
    """
    cache = dict(np.load(CACHE))
    key = f"fk_{N_MAIN:.0e}"
    optics = PercevalOptics()
    out = {k: np.full((len(distances), 4), np.nan) for k in ("no_decoy", key, "inf_decoy")}
    for i, L in enumerate(distances):
        for k, fn in [("no_decoy", q.no_decoy_rate), ("inf_decoy", q.infinite_decoy_rate)]:
            r_m, e_m, mu = q.optimize_mu(fn, L, dev, DET)
            with use_optics(optics):
                r_p, e_p = fn(L, dev, DET, mu)
            out[k][i] = r_m, r_p, e_m, e_p
        x = cache[key + "_params"][int(np.argmin(abs(cache["L"] - L)))]
        r_m, e_m = q.finite_key_rate(L, dev, DET, N_MAIN, x)
        with use_optics(optics):
            r_p, e_p = q.finite_key_rate(L, dev, DET, N_MAIN, x)
        out[key][i] = r_m, r_p, e_m, e_p

    print("=== SKR and QBER: Python model vs Perceval (same protocol parameters) ===")
    for k, v in out.items():
        ok = (v[:, 0] > 0) & (v[:, 1] > 0)
        dr = np.max(np.abs(v[ok, 1] / v[ok, 0] - 1))
        de = np.max(np.abs(v[ok, 3] / v[ok, 2] - 1))
        print(f" {k:10s}: max relative difference  SKR {dr:.1e}   QBER {de:.1e}   ({ok.sum()} distances)")
    return out


def check_textbook():
    print("\n=== Textbook BB84 behaviour (single photons, 0 km) ===")
    res = {}
    # matched / mismatched bases
    for match in (True, False):
        e = []
        for b, bit in STATES:
            bb = b if match else 1 - b
            P = click_probs(link(b, bit, bb, 0), 1)
            e.append(P[1 - bit] / P[:2].sum())
        res["matched" if match else "mismatched"] = np.mean(e)
    # intercept-resend: Eve measures in a random basis with an ideal analyzer and resends
    e_ir = []
    for b, bit in STATES:
        for eve in (0, 1):
            rot = np.pi / 4 if eve else 0.0
            pe = pcvl.Processor("SLOS", 2)  # Eve's analyzer on Alice's fiber output
            add_alice(pe, b, bit)
            rails_to_fiber(pe)
            if rot:
                pe.add(1, PR(rot))
            pe.add(0, PBS())  # H -> mode 0, V -> mode 1
            pe.min_detected_photons_filter(1)
            pe.with_input(pcvl.BasicState("|1{P:H},0>"))
            for s, pr in pe.probs()["results"].items():
                pol = "H" if s[0] else "V"
                pb = pcvl.Processor("SLOS", 2)  # Eve resends her result into the fiber mode
                if rot:
                    pb.add(1, PR(-rot))
                fiber_to_rails(pb)
                add_lossy(pb, bob_T(b))
                pb.min_detected_photons_filter(0)
                pb.with_input(pcvl.BasicState(f"|0,1{{P:{pol}}}>"))
                P = np.zeros(2)
                for so, x in pb.probs()["results"].items():
                    P += np.array([so[0] > 0, so[1] > 0]) * x
                e_ir.append(pr * P[1 - bit] / P.sum() / 2)
    res["intercept-resend"] = np.sum(e_ir) / len(STATES)
    # residual polarization rotation in the fiber
    for deg in (2, 5):
        e = []
        for b, bit in STATES:
            P = click_probs(link(b, bit, b, 0, rot=np.radians(deg)), 1)
            e.append(P[1 - bit] / P[:2].sum())
        res[f"rotation {deg} deg"] = np.mean(e)
    model = lambda rot: np.mean(bm.device_summary(rot)["e_det"])  # noqa: E731
    expected = {"matched": f"model {model(0) * 100:.3f} %", "mismatched": "50 %", "intercept-resend": "25 %"}
    for deg in (2, 5):
        expected[f"rotation {deg} deg"] = f"model {model(np.radians(deg)) * 100:.3f} %"
    for k, v in res.items():
        print(f" {k:18s}: {v * 100:7.3f} %   (expected {expected[k]})")
    return res


def figure(distances, out, res):
    cache = dict(np.load(CACHE))
    L = cache["L"]
    key = f"fk_{N_MAIN:.0e}"
    series = [("inf_decoy", BOUND, "--", "Ideal decoy (asymptotic bound)"),
              (key, BLUE, "-", "Decoy BB84, finite key (17 min block)"),
              ("no_decoy", ORANGE, "-", "No decoy (current chip)")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"wspace": 0.28})
    dot = dict(ms=7, mfc=SURFACE, mew=1.6, ls="none", zorder=5)
    for k, c, ls, _ in series:
        a.semilogy(L, np.where(cache[k] > 0, cache[k], np.nan), ls, color=c, lw=2)
        b.plot(L, np.where(cache[k] > 0, cache[k + "_qber"] * 100, np.nan), ls, color=c, lw=2)
        v = out[k]
        ok = v[:, 1] >= 1  # points still on the SKR axis
        a.plot(distances[ok], v[ok, 1], "o", mec=c, **dot)
        b.plot(distances[ok], v[ok, 3] * 100, "o", mec=c, **dot)

    handles = [Line2D([], [], color=c, ls=ls, lw=2, label=lab) for _, c, ls, lab in series]
    handles += [Line2D([], [], color=INK2, lw=2, label="Python model (lines)"),
                Line2D([], [], marker="o", mec=INK, mfc=SURFACE, mew=1.6, ls="none", ms=7,
                       label="Perceval (points)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9.5, bbox_to_anchor=(0.5, -0.06))
    a.set_ylim(1, 1e6)
    a.set_xlim(0, 200)
    a.set_xlabel("Fiber distance (km)")
    a.set_ylabel("Secret key rate (bit/s)")
    a.set_title("(a) Secret key rate", loc="left", fontsize=12)

    b.axhline(11, color=INK2, lw=1, ls=":")
    b.text(198, 11.2, "BB84 threshold (11 %)", ha="right", va="bottom", color=INK2, fontsize=9.5)
    b.text(4, 9.6, "Perceval also reproduces textbook BB84:\n"
           f"mismatched bases {res['mismatched'] * 100:.1f} % (expected 50 %)\n"
           f"intercept-resend attack {res['intercept-resend'] * 100:.1f} % (expected 25 %)",
           va="top", color=INK2, fontsize=9.5, linespacing=1.5)
    b.set_xlim(0, 200)
    b.set_ylim(0, 13)
    b.set_xlabel("Fiber distance (km)")
    b.set_ylabel("QBER (%)")
    b.set_title("(b) Quantum bit error rate", loc="left", fontsize=12)
    fig.savefig(FIG / "fig_perceval_check.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    dev = q.load_device()
    distances = np.arange(0, 201, 10.0)
    out = check_rates(dev, distances)
    res = check_textbook()
    figure(distances, out, res)
    print("\nsaved", FIG / "fig_perceval_check.png")
