"""Benchmark the key-rate formulas in qkd_model.py against published results.

The security formulas used for the talk figures (no_decoy_rate, infinite_decoy_rate,
key_length + maximize_protocol) are run with each paper's own parameters and compared
with curves digitized from the papers' figures (benchmark_data/*.csv, see digitize_eps.py)
and with the numbers stated in their text. PDFs are in references/.

  [A] Ma, Qi, Zhao, Lo, PRA 72, 012326 (2005), Fig. 3 -- GYS experiment, asymptotic:
      infinite decoys (mu = 0.48) and vacuum + weak decoy (mu = 0.48, nu = 0.05, 0),
      standard BB84 (q = 1/2), f = 1.22. Text: 142.05 km and 140.55 km.
  [B] Lo, Ma, Chen, PRL 94, 230504 (2005), Fig. 1 -- GYS, GLLP without decoys. Text:
      "only about 30 km". The paper does not state mu for this curve, so this is only a
      consistency check (console output, not in the figure).
  [C] Lim, Curty, Walenta, Xu, Zbinden, PRA 89, 022307 (2014), Fig. 1 -- finite-key
      3-intensity decoy BB84, fixed key-basis block n_X = 10^4..10^9, eps_sec = kappa*ell.
      Text: n_X = 10^4 reaches 135 km; at 100 km R(10^9) / R(10^7) ~ 1.75.

[A] and [B] use the detector model of qkd_model.py (OPTICS + add_noise) with the paper's
parameters. [C] uses the detection and error model written in Lim et al. (their
"dedicated fiber" equations) so that only the security bound is being compared.
"""
import csv
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import brentq

import qkd_model as q
from make_figures import BLUE, BOUND, FIG, INK, INK2, SURFACE

HERE = Path(__file__).resolve().parent
DATA = HERE / "benchmark_data"


@contextmanager
def settings(**kw):
    """Temporarily override module-level constants of qkd_model (ALPHA, F_EC, MU3, ...)."""
    saved = {k: getattr(q, k) for k in kw}
    for k, v in kw.items():
        setattr(q, k, v)
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(q, k, v)


def load_digitized(name):
    """{style: [segments (km, rate) sorted by distance]} merged per line style."""
    rows = [r for r in csv.reader(open(DATA / f"{name}_digitized.csv")) if not r[0].startswith("#")][1:]
    segs = {}
    for style, i, km, rate in rows:
        segs.setdefault((style, int(i)), []).append((float(km), float(rate)))
    out = {}
    for (style, _), pts in segs.items():
        out.setdefault(style, []).append(np.array(sorted(pts)))
    return out


def merged(segments):
    pts = np.vstack(segments)
    return pts[np.argsort(pts[:, 0])]


def interp_log(curve, L):
    return 10 ** np.interp(L, curve[:, 0], np.log10(curve[:, 1]), left=np.nan, right=np.nan)


def zero_crossing(f, lo, hi):
    """Largest distance with positive rate (f decreasing), by bisection."""
    return brentq(lambda L: f(L) - 1e-300, lo, hi, xtol=1e-3) if f(lo) > 0 >= f(hi) - 1e-300 else np.nan


# ==========================================
# [A] + [B]: GYS experiment, asymptotic formulas
# ==========================================
GYS = dict(alpha=0.21, e_det=0.033, Y0=1.7e-6, eta_bob=0.045)  # Ma et al. 2005, Table I


def gys():
    dev = q.Device(np.array([1.0, 1.0]), np.array([GYS["e_det"]] * 2))
    det = q.Detector("GYS", eta=GYS["eta_bob"], p_dark=GYS["Y0"] / 2, p_ap=0.0, t_dead=0.0)
    return dev, det


def ma_infinite(L, mu=0.48):
    dev, det = gys()
    with settings(ALPHA=GYS["alpha"], F_EC=1.22):
        return 0.5 * q.infinite_decoy_rate(L, dev, det, mu)[0] / q.F_REP  # q = 1/2, per pulse


def ma_vacuum_weak(L, mu=0.48, nu=0.05):
    """Asymptotic vacuum + weak decoy (Ma 2005) = Lim's bound with N -> inf and mu3 = 0;
    key from signal pulses only and no vacuum term, as in Ma's Eq. for this protocol."""
    dev, det = gys()
    N, eps = 1e40, 1e-12
    mus, ps = np.array([mu, nu, 0.0]), np.array([1 - 2e-9, 1e-9, 1e-9])
    with settings(ALPHA=GYS["alpha"], F_EC=1.22):
        A, y0 = q.detector_state(mus, ps, L, dev, det)
        Q, E = np.array([q.gains(m, L, dev, det, 0, A, y0) for m in mus]).T
        n = N * ps * Q
        ell = q.key_length(mus, ps, n, n * E, n, n * E, eps_sec=eps, eps_cor=eps, vacuum_key=False)
    return max(0.5 * ell / N, 0.0) if ell is not None else 0.0


def lmc_no_decoy(L):
    dev, det = gys()
    with settings(ALPHA=GYS["alpha"], F_EC=1.0):
        return q.optimize_mu(q.no_decoy_rate, L, dev, det)[0] / q.F_REP  # q = 1, optimal mu


def optimal_mu_check():
    """Ma 2005: mu_opt ~ 0.48 (f = 1.22) and ~0.54 (f = 1) for GYS, where Y0 is negligible."""
    dev, det = gys()
    out = {}
    for f in (1.22, 1.0):
        with settings(ALPHA=GYS["alpha"], F_EC=f):
            out[f] = q.optimize_mu(q.infinite_decoy_rate, 20.0, dev, det)[2]
    return out


# ==========================================
# [C]: Lim et al. 2014, finite key
# ==========================================
LIM = dict(eta_bob=0.10, p_dc=6e-7, p_ap=0.04, e_mis=5e-3, f_ec=1.16, kappa=1e-15, eps_cor=1e-15)


def lim_counts(L, mus):
    """Lim et al. 'dedicated fiber' model: detection rate R_k and error probability e_k."""
    eta_sys = 10 ** (-0.2 * L / 10) * LIM["eta_bob"]
    D = 1 - (1 - 2 * LIM["p_dc"]) * np.exp(-eta_sys * mus)
    R = D * (1 + LIM["p_ap"])
    e = LIM["p_dc"] + LIM["e_mis"] * (1 - np.exp(-eta_sys * mus)) + LIM["p_ap"] * D / 2
    return R, e


def lim_rate(L, n_x, x):
    """Secret key per pulse for a fixed key-basis block n_x, eps_sec = kappa * ell."""
    if not q.protocol_params_ok(x):
        return -np.inf
    mu1, mu2, p1, p2, px = x
    mus, ps = np.array([mu1, mu2, q.MU3]), np.array([p1, p2, 1 - p1 - p2])
    R, e = lim_counts(L, mus)
    N = n_x / (px**2 * np.sum(ps * R))
    nK, mK = N * px**2 * ps * R, N * px**2 * ps * e
    nT, mT = N * (1 - px) ** 2 * ps * R, N * (1 - px) ** 2 * ps * e
    eps, ell = 1e-10, None
    for _ in range(8):  # fixed point of eps_sec = kappa * ell
        ell = q.key_length(mus, ps, nK, mK, nT, mT, eps_sec=eps, eps_cor=LIM["eps_cor"], f_ec=LIM["f_ec"])
        if ell is None:
            return -np.inf  # decoy bounds infeasible
        if ell <= 0:
            break  # no key: keep the (negative) value so the optimizer can climb out
        eps = min(LIM["kappa"] * ell, 1e-3)
    return ell / N


def lim_curve(n_x, distances):
    x = [0.5, 0.1, 0.6, 0.3, 0.9]
    out = []
    for L in distances:
        xn = q.maximize_protocol(lambda x: lim_rate(L, n_x, x) * 1e8, [x] + q.DEFAULT_STARTS)
        r = lim_rate(L, n_x, xn)
        if r > 0:
            x = xn
        out.append(max(r, 0.0))
    return np.array(out)


def crossing(L, r, level=1e-7):
    """Distance where the rate falls to `level` (Lim's figure stops at its 1e-7 axis)."""
    ok = r > 0
    L, lr = np.asarray(L)[ok], np.log10(r[ok])
    i = np.nonzero(lr < np.log10(level))[0]
    if not len(i):  # drops from above the level straight to zero within one grid step
        return L[-1]
    j = i[0]
    return np.interp(np.log10(level), [lr[j], lr[j - 1]], [L[j], L[j - 1]])


# ==========================================
# MAIN
# ==========================================
def main():
    ma, lmc, lim = load_digitized("ma2005"), load_digitized("lmc2005"), load_digitized("lim2014")
    results = {}

    print("=== [A] Ma et al. 2005, GYS, asymptotic (rates per pulse, q = 1/2) ===")
    for name, fn, style, text in [("infinite decoys", ma_infinite, "DA", 142.05),
                                  ("vacuum + weak decoy", ma_vacuum_weak, "SO", 140.55)]:
        paper = merged(ma[style])
        Lmax = zero_crossing(fn, 1.0, 200.0)
        Ls = np.arange(0, 141, 10.0)
        mine = np.array([fn(L) for L in Ls])
        ratio = mine / interp_log(paper, Ls)
        results[name] = (Ls, mine, paper)
        print(f" {name:20s}: max distance {Lmax:6.2f} km (paper text {text} km, figure {paper[-1, 0]:.1f} km); "
              f"rate / paper curve {np.nanmin(ratio):.3f}-{np.nanmax(ratio):.3f}")
    mus = optimal_mu_check()
    print(f" optimal mu at 20 km: {mus[1.22]:.3f} for f = 1.22 (paper 0.48), {mus[1.0]:.3f} for f = 1 (paper 0.54)")

    print("\n=== [B] Lo, Ma, Chen 2005, GYS, GLLP without decoys (f = 1) ===")
    paper = merged(lmc["DA"])
    print(f" paper figure ends at {paper[-1, 0]:.1f} km (text: 'only about 30 km'); mu is not stated in the paper")
    print(f" optimal mu at every distance : max distance {zero_crossing(lmc_no_decoy, 1.0, 80.0):6.2f} km")
    dev, det = gys()
    with settings(ALPHA=GYS["alpha"], F_EC=1.0):
        fixed = lambda L: q.no_decoy_rate(L, dev, det, 0.012)[0] / q.F_REP  # noqa: E731
        Ls = np.array([5, 10, 15, 20, 25.0])
        ratio = 0.5 * np.array([fixed(L) for L in Ls]) / interp_log(paper, Ls)
        print(f" fixed mu = 0.012, q = 1/2    : max distance {zero_crossing(fixed, 1.0, 80.0):6.2f} km, "
              f"rate / paper curve {ratio.min():.2f}-{ratio.max():.2f} (consistent, but settings not stated)")

    print("\n=== [C] Lim et al. 2014, finite key (rates per pulse) ===")
    solid = sorted(lim["SO"], key=lambda c: c[-1, 0])  # n_X = 1e4 ... 1e9 by reach
    Ls = np.arange(0, 201, 5.0)
    lim_res = {}
    for s, paper in zip(range(4, 10), solid):
        r = lim_curve(10.0**s, Ls)
        lim_res[s] = (r, paper)
        ratio = r / interp_log(paper, Ls)
        print(f" n_X = 1e{s}: reach to 1e-7 {crossing(Ls, r):6.1f} km (paper figure {paper[-1, 0]:.1f} km); "
              f"rate / paper curve {np.nanmin(ratio[Ls <= paper[-1, 0] - 5]):.3f}-{np.nanmax(ratio[Ls <= paper[-1, 0] - 5]):.3f}")
    r9, r7 = lim_res[9][0][Ls == 100][0], lim_res[7][0][Ls == 100][0]
    print(f" R(1e9) / R(1e7) at 100 km: {r9 / r7:.2f} (paper text ~1.75)")
    r_inf = lim_curve(1e25, Ls)
    paper_inf = merged(lim["DD"])
    ratio = r_inf / interp_log(paper_inf, Ls)
    print(f" asymptotic: reach to 1e-7 {crossing(Ls, r_inf):6.1f} km (paper figure {paper_inf[-1, 0]:.1f} km); "
          f"rate / paper curve {np.nanmin(ratio[Ls <= 185]):.3f}-{np.nanmax(ratio[Ls <= 185]):.3f}")

    figure(results, Ls, lim_res, r_inf, paper_inf)


def figure(results, Ls, lim_res, r_inf, paper_inf):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"wspace": 0.25})
    dot = dict(ms=6.5, mfc=SURFACE, mew=1.5, ls="none", zorder=5)

    for name, c, ls in [("infinite decoys", BOUND, "--"), ("vacuum + weak decoy", BLUE, "-")]:
        L, mine, paper = results[name]
        a.semilogy(paper[:, 0], paper[:, 1], ls, color=c, lw=2)
        ok = mine > 1e-8
        a.plot(L[ok], mine[ok], "o", mec=c, **dot)
    a.set_xlim(0, 160)
    a.set_ylim(1e-7, 1e-2)
    a.set_xlabel("Fiber distance (km)")
    a.set_ylabel("Secret key rate (per pulse)")
    a.set_title("(a) Asymptotic decoy: Ma et al. 2005, Fig. 3 (GYS)", loc="left", fontsize=11.5)
    a.legend(handles=[Line2D([], [], color=BOUND, ls="--", lw=2, label="Ideal decoy (Ma Fig. 3)"),
                      Line2D([], [], color=BLUE, lw=2, label="Vacuum + weak decoy (Ma Fig. 3)"),
                      Line2D([], [], marker="o", mec=INK, mfc=SURFACE, mew=1.5, ls="none", ms=6.5,
                             label="This model, paper's parameters")],
             loc="upper right", fontsize=9)

    ramp = ["#86b6ef", "#5598e7", "#3987e5", "#2a78d6", "#1c5cab", "#104281"]  # n_X = 1e4 ... 1e9
    for (s, (r, paper)), c in zip(lim_res.items(), ramp):
        b.semilogy(paper[:, 0], paper[:, 1], "-", color=c, lw=1.8)
        ok = r >= 1e-7
        b.plot(Ls[ok], r[ok], "o", mec=c, **dot)
    b.semilogy(paper_inf[:, 0], paper_inf[:, 1], "--", color=BOUND, lw=1.8)
    ok = r_inf >= 1e-7
    b.plot(Ls[ok], r_inf[ok], "o", mec=BOUND, **dot)
    b.set_xlim(0, 200)
    b.set_ylim(1e-7, 2e-2)
    b.set_xlabel("Fiber distance (km)")
    b.set_title("(b) Finite key: Lim et al. 2014, Fig. 1 (block size $n_X$)", loc="left", fontsize=11.5)
    sizes = b.legend(handles=[Line2D([], [], color=c, lw=2, label=f"$n_X = 10^{s}$") for s, c in zip(lim_res, ramp)],
                     loc="lower left", fontsize=9, title="Block size", title_fontsize=9)
    b.add_artist(sizes)
    b.legend(handles=[Line2D([], [], color=INK2, lw=2, label="Published curve (digitized)"),
                      Line2D([], [], marker="o", mec=INK, mfc=SURFACE, mew=1.5, ls="none", ms=6.5,
                             label="This model, paper's parameters"),
                      Line2D([], [], color=BOUND, ls="--", lw=2, label="Asymptotic")],
             loc="upper right", fontsize=9)
    fig.savefig(FIG / "fig_benchmark_published.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("\nsaved", FIG / "fig_benchmark_published.png")


if __name__ == "__main__":
    main()
