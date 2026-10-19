"""BB84 key-rate model driven by the simulated circuit (bob_model.py).

Chain:  Alice chip -> mu at Alice's fiber output -> fiber (alpha dB/km) -> Bob chip -> 2 SPADs

* Device: Bob's transmission and the matched-basis polarization error per basis come
  from the INTERCONNECT-validated circuit model. Z = rectilinear (key basis, lowest
  error), X = diagonal (test basis). Alice's insertion loss does not enter: mu is set
  at her output, so her loss only costs laser power.
* Detector: gated InGaAs SPAD with efficiency, dark counts, afterpulsing and
  (non-paralyzable) dead time; double clicks get a random bit (squashing).
* Protocols:
    no_decoy_rate        asymptotic GLLP, no decoys (all multi-photon pulses tagged)
    infinite_decoy_rate  asymptotic GLLP with perfectly known single-photon yield (bound)
    finite_key_rate      3-intensity decoy BB84, finite key, Lim et al. PRA 89, 022307 (2014)
"""
import sys  
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lumerical/interconnect/bob_test"))
import bob_model  # noqa: E402

# ==========================================
# PARAMETERS
# ==========================================
F_REP = 100e6  # pulse rate (Hz)
ALPHA = 0.2  # fiber loss (dB/km)
F_EC = 1.16  # error-correction inefficiency
EPS_SEC, EPS_COR = 1e-10, 1e-15
MU3 = 2e-4  # weakest decoy ("vacuum") intensity


@dataclass(frozen=True)
class Detector:
    """Per-detector parameters. Values are typical, not a specific product."""
    name: str
    eta: float  # detection efficiency (gate/jitter losses folded in)
    p_dark: float  # dark-count probability per gate
    p_ap: float  # afterpulse probability per click
    t_dead: float  # dead time after a click (s)


SPAD = Detector("InGaAs SPAD (gated)", eta=0.20, p_dark=1e-6, p_ap=0.01, t_dead=10e-6)
SNSPD = Detector("SNSPD", eta=0.85, p_dark=1e-7, p_ap=0.0, t_dead=50e-9)


@dataclass(frozen=True)
class Device:
    eta_bob: np.ndarray  # Bob transmission [Z, X]
    e_det: np.ndarray  # matched-basis polarization error [Z, X]


def load_device(misalign_deg=0.0):
    """Device from the circuit model, plus an optional residual fiber misalignment
    (polarization rotation left after compensation). The rotation adds to the chip's
    error as a field; its sign is unknown in practice, so +/- rotations are averaged."""
    th = np.radians(misalign_deg)
    ds = [bob_model.device_summary(r) for r in ({th, -th} if th else {0.0})]
    return Device(np.mean([d["eta_bob"] for d in ds], axis=0), np.mean([d["e_det"] for d in ds], axis=0))


def h2(x):
    x = np.clip(x, 1e-15, 0.5)
    return -x * np.log2(x) - (1 - x) * np.log2(1 - x)


# ==========================================
# DETECTION
# ==========================================
def optical_wcp(mu, L, dev, det, basis):
    """Optical click probabilities (correct only, wrong only, both) for a phase-randomized
    weak coherent pulse of mean mu, before detector noise. Poisson light splits into
    independent Poisson streams, so the two detectors fire independently."""
    n = mu * transmittance(L, dev, det, basis)
    e = dev.e_det[basis]
    c, w = 1 - np.exp(-n * (1 - e)), 1 - np.exp(-n * e)
    return np.array([c * (1 - w), w * (1 - c), c * w])


def optical_single(L, dev, det, basis):
    """Same for exactly one photon: it reaches the correct detector, the wrong one, or is lost."""
    eta, e = transmittance(L, dev, det, basis), dev.e_det[basis]
    return np.array([eta * (1 - e), eta * e, 0.0])


# Swappable source of optical click probabilities (perceval_check.py plugs in Perceval here).
OPTICS = {"wcp": optical_wcp, "single": optical_single}


def add_noise(p, y0):
    """Add independent noise clicks (prob y0 per detector) to optical click probabilities p.
    Returns (click probability, error probability); double clicks get a random bit."""
    pc, pw, pb = p
    none = 1 - pc - pw - pb
    Pc = pc * (1 - y0) + none * y0 * (1 - y0)
    Pw = pw * (1 - y0) + none * y0 * (1 - y0)
    Pb = pb + (pc + pw) * y0 + none * y0**2
    Q = Pc + Pw + Pb
    return Q, (Pw + 0.5 * Pb) / Q


def detector_state(mus, probs, L, dev, det):
    """Self-consistent dead-time availability A and noise click probability y0.

    The click rate (averaged over the intensity mix, key basis) sets how often the
    detectors are dead and how many afterpulses they produce.
    """
    opt = [OPTICS["wcp"](m, L, dev, det, 0) for m in mus]
    A, y_ap = 1.0, 0.0
    for _ in range(30):
        y0 = det.p_dark + y_ap
        q_u = sum(p * add_noise(o, y0)[0] for o, p in zip(opt, probs))
        A = 1 / (1 + F_REP * q_u / 2 * det.t_dead)  # per detector, half the clicks each
        y_ap = det.p_ap * A * q_u / 2
    return A, det.p_dark + y_ap


def gains(mu, L, dev, det, basis, A, y0):
    """Observed gain (clicks per pulse, after dead time) and QBER."""
    Q, E = add_noise(OPTICS["wcp"](mu, L, dev, det, basis), y0)
    return A * Q, E


def single_photon(L, dev, det, basis, y0):
    """Single-photon yield and error (dead time applied by the caller)."""
    return add_noise(OPTICS["single"](L, dev, det, basis), y0)


def transmittance(L, dev, det, basis):
    return 10 ** (-ALPHA * L / 10) * dev.eta_bob[basis] * det.eta


# ==========================================
# ASYMPTOTIC PROTOCOLS (efficient BB84, sifting -> 1)
# ==========================================
def no_decoy_rate(L, dev, det, mu):
    A, y0 = detector_state([mu], [1.0], L, dev, det)
    QZ, EZ = gains(mu, L, dev, det, 0, A, y0)
    QX, EX = gains(mu, L, dev, det, 1, A, y0)
    delta = (1 - np.exp(-mu) * (1 + mu)) / QZ  # multi-photon fraction of detections
    if delta >= 1:
        return 0.0, EZ
    r = QZ * ((1 - delta) * (1 - h2(EX / (1 - delta))) - F_EC * h2(EZ))
    return max(r, 0.0) * F_REP, EZ


def infinite_decoy_rate(L, dev, det, mu):
    A, y0 = detector_state([mu], [1.0], L, dev, det)
    QZ, EZ = gains(mu, L, dev, det, 0, A, y0)
    Y1, _ = single_photon(L, dev, det, 0, y0)
    _, e1X = single_photon(L, dev, det, 1, y0)
    Q1 = A * mu * np.exp(-mu) * Y1
    r = Q1 * (1 - h2(e1X)) - F_EC * QZ * h2(EZ)
    return max(r, 0.0) * F_REP, EZ


def optimize_mu(rate_fn, L, dev, det):
    res = minimize_scalar(lambda lm: -rate_fn(L, dev, det, 10**lm)[0], bounds=(-4, 0.3), method="bounded")
    mu = 10**res.x
    r, E = rate_fn(L, dev, det, mu)
    return r, E, mu


# ==========================================
# FINITE-KEY DECOY BB84 (Lim et al. 2014)
# ==========================================
def key_length(mus, ps, nK, mK, nT, mT, eps_sec=EPS_SEC, eps_cor=EPS_COR, f_ec=None, vacuum_key=True):
    """Secret key length, Lim et al. PRA 89, 022307 (2014), Eqs. (1)-(5).

    mus, ps : the three intensities (mu1 > mu2 + mu3, mu2 > mu3) and their probabilities
    nK, mK  : key-basis detections and errors per intensity (Lim's basis X; Z here)
    nT, mT  : test-basis detections and errors per intensity (Lim's basis Z; X here)
    vacuum_key: count key from vacuum events (Lim et al.); Ma et al. 2005 leave it out.
    Returns ell (bits), or None if the decoy bounds are infeasible.
    """
    f_ec = F_EC if f_ec is None else f_ec
    mu1, mu2, mu3 = mus
    log_term = np.log(21 / eps_sec)  # Hoeffding with eps_1 = eps_sec / 21

    def hoeff(n_k):
        d = np.sqrt(n_k.sum() / 2 * log_term)
        return np.exp(mus) / ps * (n_k - d), np.exp(mus) / ps * (n_k + d)

    tau = lambda n: np.sum(ps * np.exp(-mus) * mus**n / np.prod(np.arange(1, n + 1)))  # noqa: E731
    t0, t1 = tau(0), tau(1)
    den = mu1 * (mu2 - mu3) - mu2**2 + mu3**2

    def s01(n):
        lo, hi = hoeff(n)
        s0 = max(t0 * (mu2 * lo[2] - mu3 * hi[1]) / (mu2 - mu3), 0.0)  # Eq. (2)
        s1 = t1 * mu1 * (lo[1] - hi[2] - (mu2**2 - mu3**2) / mu1**2 * (hi[0] - s0 / t0)) / den  # Eq. (3)
        return s0, s1

    sK0, sK1 = s01(nK)
    _, sT1 = s01(nT)
    mlo, mhi = hoeff(mT)
    vT1 = t1 * (mhi[1] - mlo[2]) / (mu2 - mu3)  # Eq. (4)
    if sK1 <= 0 or sT1 <= 0 or vT1 <= 0:
        return None
    b = min(vT1 / sT1, 0.5)
    # statistical penalty on the phase error; the log turns negative only for astronomically
    # large blocks, where the penalty vanishes (gamma -> 0)
    arg = (sK1 + sT1) / (sK1 * sT1 * (1 - b) * b) * 21**2 / eps_sec**2
    gamma = np.sqrt((sK1 + sT1) * (1 - b) * b / (sK1 * sT1 * np.log(2)) * np.log2(max(arg, 1.0)))
    phi = min(b + gamma, 0.5)  # Eq. (5)
    leak = f_ec * nK.sum() * h2(mK.sum() / nK.sum())
    s0 = sK0 if vacuum_key else 0.0
    return s0 + sK1 * (1 - h2(phi)) - leak - 6 * np.log2(21 / eps_sec) - np.log2(2 / eps_cor)  # Eq. (1)


def protocol_params_ok(x):
    mu1, mu2, p1, p2, pZ = x
    return mu1 > mu2 + MU3 and mu2 > MU3 and p1 > 0 and p2 > 0 and 1 - p1 - p2 > 0 and 0 < pZ < 1


def finite_key_rate(L, dev, det, N, x):
    """x = (mu1, mu2, p_mu1, p_mu2, p_Z). Returns (key rate in bit/s, Z-basis QBER).
    The rate is -inf when the decoy bounds are infeasible and negative when no key is left."""
    if not protocol_params_ok(x):
        return -np.inf, np.nan
    mu1, mu2, p1, p2, pZ = x
    mus, ps = np.array([mu1, mu2, MU3]), np.array([p1, p2, 1 - p1 - p2])
    A, y0 = detector_state(mus, ps, L, dev, det)
    QZ, EZ = np.array([gains(m, L, dev, det, 0, A, y0) for m in mus]).T
    QX, EX = np.array([gains(m, L, dev, det, 1, A, y0) for m in mus]).T
    nZ = N * ps * pZ**2 * QZ  # expected sifted counts per intensity
    nX = N * ps * (1 - pZ) ** 2 * QX
    ell = key_length(mus, ps, nZ, nZ * EZ, nX, nX * EX)
    if ell is None:
        return -np.inf, np.nan  # decoy bounds infeasible
    return ell * F_REP / N, np.sum(nZ * EZ) / nZ.sum()  # may be negative: no key at this point


def maximize_protocol(objective, starts, n_random=400, n_polish=3, seed=0):
    """Maximize objective(x) over x = (mu1, mu2, p_mu1, p_mu2, p_Z).

    1. Random search (fixed seed, so results are reproducible) plus the given starts: the
       optimum can sit far from intuition (e.g. small blocks favour sending mostly decoys).
    2. Nelder-Mead polishing from the best few candidates, in an unconstrained
       parametrization that keeps mu3 < mu2 < mu1, probabilities in (0, 1).
    objective(x) returns the value to maximize (-inf when the decoy bounds are infeasible).
    """
    def unpack(z):
        mu1 = 10 ** z[0]
        mu2 = MU3 + (mu1 - MU3) / (1 + np.exp(-z[1]))  # MU3 < mu2 < mu1
        p = np.exp([z[2], z[3], 0.0]); p /= p.sum()
        return np.array([mu1, mu2, p[0], p[1], 1 / (1 + np.exp(-z[4]))])

    def pack(x):
        mu1, mu2, p1, p2, pZ = x
        f = (mu2 - MU3) / (mu1 - MU3)
        p3 = 1 - p1 - p2
        return np.array([np.log10(mu1), np.log(f / (1 - f)), np.log(p1 / p3), np.log(p2 / p3), np.log(pZ / (1 - pZ))])

    def cost(z):
        with np.errstate(all="ignore"):  # extreme random samples overflow harmlessly
            v = objective(unpack(z))
        return -v if np.isfinite(v) else 1e30  # finite penalty keeps Nelder-Mead well defined

    rng = np.random.default_rng(seed)
    cands = [pack(x) for x in starts if protocol_params_ok(x)]
    cands += [np.array([rng.uniform(-1.5, 0.2), rng.normal(0, 2), rng.normal(0, 2), rng.normal(0, 2), rng.normal(1, 2)])
              for _ in range(n_random)]
    scored = sorted(((cost(z), i) for i, z in enumerate(cands)))
    best = None
    for c0, i in scored[:n_polish]:
        res = minimize(cost, cands[i], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 4000})
        if best is None or res.fun < best.fun:
            best = res
    return unpack(best.x)


DEFAULT_STARTS = [[0.5, 0.1, 0.6, 0.3, 0.9], [0.5, 0.1, 0.6, 0.3, 0.7], [0.4, 0.1, 0.5, 0.3, 0.55]]


def optimize_finite(L, dev, det, N, x0):
    """Maximize the finite-key rate over intensities, their probabilities and p_Z."""
    x = maximize_protocol(lambda x: finite_key_rate(L, dev, det, N, x)[0] / F_REP * N, [x0] + DEFAULT_STARTS)
    r, E = finite_key_rate(L, dev, det, N, x)
    return max(r, 0.0), E, x


# ==========================================
# SWEEPS
# ==========================================
def sweep(distances, dev, det, N_list=()):
    """SKR (bit/s) and Z-QBER vs distance for all protocols."""
    out = {"L": np.asarray(distances)}
    for name, fn in [("no_decoy", no_decoy_rate), ("inf_decoy", infinite_decoy_rate)]:
        r, e = np.array([optimize_mu(fn, L, dev, det)[:2] for L in distances]).T
        out[name], out[name + "_qber"] = r, e
    for N in N_list:
        x = np.array([0.5, 0.1, 0.6, 0.3, 0.9])
        rs, es, xs = [], [], []
        for L in distances:
            r, e, xn = optimize_finite(L, dev, det, N, x)
            if r > 0:
                x = xn  # warm start the next distance
            rs.append(r); es.append(e); xs.append(xn)
        key = f"fk_{N:.0e}"
        out[key], out[key + "_qber"], out[key + "_params"] = np.array(rs), np.array(es), np.array(xs)
    return out


def max_range(L, r):
    pos = np.asarray(L)[np.asarray(r) > 0]
    return pos.max() if len(pos) else 0.0


if __name__ == "__main__":
    dev = load_device()
    print("device: eta_bob", dev.eta_bob.round(4), " e_det", dev.e_det)
    for L in [0, 50, 100]:
        print(L, "km  no-decoy", optimize_mu(no_decoy_rate, L, dev, SPAD)[0],
              " inf-decoy", optimize_mu(infinite_decoy_rate, L, dev, SPAD)[0],
              " finite 1e11", optimize_finite(L, dev, SPAD, 1e11, [0.5, 0.1, 0.6, 0.3, 0.9])[0])
