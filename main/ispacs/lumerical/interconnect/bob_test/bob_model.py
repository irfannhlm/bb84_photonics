"""Python rebuild of the BB84 Alice/Bob circuit (segmented MZM + 2D PSGC) at one wavelength.

The circuit is the same compound element used for both Alice and Bob in INTERCONNECT:

    MMI_1.opt_4 -> wg_heater_1 (100 um, basis) -> wg_heater_5 (200 um, bit) -> MMI_2.opt_2   (top arm)
    MMI_1.opt_3 -> wg_heater_4 (100 um, 0 V)   -> wg_heater_6 (200 um, 0 V) -> MMI_2.opt_1   (bottom arm)
    MMI_2.opt_3 -> wg_heater_3 (300 um, 0 V) -> PSGC.WG_X
    MMI_2.opt_4 -> wg_heater_2 (300 um, 0 V) -> PSGC.WG_Y
    PSGC.Fiber  (mode 1, mode 2)

External ports of the solved network: [MMI_1.opt_1, MMI_1.opt_2, Fiber mode 1, Fiber mode 2].
The full S-matrix network (including reflections) is solved, like INTERCONNECT does.

Steps:
  1. Calibrate the heater (phase per V^2*L and loss per length) on Alice's sweep_basis.
  2. Validate, with no further fitting, against Bob's sweep_base (8 points).
  3. Predict Bob's response to elliptical/circular light and the Alice -> Bob error matrix.

Data: run export_reference.py once (needs Lumerical) to create interconnect_reference.npz.
"""
import re
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
PSGC_DAT = REPO / "tests/lumerical/interconnect/psgc_lumerical.dat"  # as referenced by SPAR_1
MMI_DAT = REPO / "tests/lumerical/CML_Build/BB84_Custom_CML/source/MMI_2x2_TE_c/MMI_2x2_TE_c.dat"
REF_NPZ = HERE / "interconnect_reference.npz"

C0 = 299792458.0
V_HIGH = 1.2  # logic-high drive voltage (V)
L_REF = 300e-6  # PDK reference: 3pi/2 at 1.2 V for 300 um


# ============================================================
# S-PARAMETER FILES
# ============================================================
def load_sparam(path):
    """Read a Lumerical .dat S-parameter file.

    Each block header is (out_port, mode_label, out_mode, in_port, in_mode, ...); the
    out/in order was confirmed against the Alice (WG->fiber) and Bob (fiber->WG) losses,
    because the PSGC data is not reciprocal.
    """
    lines = Path(path).read_text().splitlines()
    blocks, i = {}, 0
    while i < len(lines):
        m = re.match(r"\s*\((.*)\)\s*$", lines[i])
        if m and "transmission" in lines[i]:
            t = [x.strip().strip('"') for x in m.group(1).split(",")]
            n = int(re.match(r"\s*\((\d+)", lines[i + 1]).group(1))
            arr = np.array([[float(v) for v in lines[i + 2 + k].split()] for k in range(n)])
            blocks[(t[0], int(t[2]), t[3], int(t[4]))] = arr
            i += 2 + n
        else:
            i += 1
    return blocks


def sparam_matrix(blocks, ports, f0):
    """Complex S-matrix at f0 (linear interpolation of |S| and unwrapped phase)."""
    S = np.zeros((len(ports), len(ports)), complex)
    for (po, mo, pi, mi), a in blocks.items():
        mag = np.interp(f0, a[:, 0], a[:, 1])
        ph = np.interp(f0, a[:, 0], np.unwrap(a[:, 2]))
        S[ports.index((po, mo)), ports.index((pi, mi))] = mag * np.exp(1j * ph)
    return S


# ============================================================
# COMPONENTS
# ============================================================
PSGC_BLOCKS = load_sparam(PSGC_DAT)
MMI_BLOCKS = load_sparam(MMI_DAT)


def psgc(f0):
    ports = [("WG_X", 1), ("WG_Y", 1), ("Fiber", 1), ("Fiber", 2)]
    return ["WG_X", "WG_Y", "F1", "F2"], sparam_matrix(PSGC_BLOCKS, ports, f0)


def mmi(f0):
    ports = [(f"port {i}", 1) for i in range(1, 5)]
    return ["opt_1", "opt_2", "opt_3", "opt_4"], sparam_matrix(MMI_BLOCKS, ports, f0)


def heater(length, voltage, k_phase, loss_db_per_m):
    """Thermo-optic shifter: phase = k * V^2 * L (heater power ~ V^2), loss ~ L.

    The static propagation phase is dropped: both MZM arms (100+200 um) and both
    output heaters (300 um) have equal length, so it cancels.
    """
    t = 10 ** (-loss_db_per_m * length / 20) * np.exp(1j * k_phase * voltage**2 * length)
    return ["opt_1", "opt_2"], np.array([[0, t], [t, 0]])


# ============================================================
# NETWORK SOLVER
# ============================================================
def solve_network(components, connections, external):
    """Reduce connected S-matrices to the S-matrix seen at the external ports."""
    index, blocks = {}, []
    for name, (ports, S) in components.items():
        for p in ports:
            index[f"{name}.{p}"] = len(index)
        blocks.append(S)
    n = len(index)
    S = np.zeros((n, n), complex)
    o = 0
    for b in blocks:
        S[o:o + len(b), o:o + len(b)] = b
        o += len(b)

    internal = [index[p] for pair in connections for p in pair]
    ext = [index[p] for p in external]
    assert len(set(internal) | set(ext)) == n, "every port must be connected or external"
    pos = {g: k for k, g in enumerate(internal)}
    P = np.zeros((len(internal), len(internal)))
    for a, b in connections:
        P[pos[index[a]], pos[index[b]]] = P[pos[index[b]], pos[index[a]]] = 1

    S_ii, S_ie = S[np.ix_(internal, internal)], S[np.ix_(internal, ext)]
    S_ei, S_ee = S[np.ix_(ext, internal)], S[np.ix_(ext, ext)]
    return S_ee + S_ei @ P @ np.linalg.solve(np.eye(len(internal)) - S_ii @ P, S_ie)


CONNECTIONS = [
    ("MMI1.opt_4", "h1.opt_1"), ("h1.opt_2", "h5.opt_1"), ("h5.opt_2", "MMI2.opt_2"),
    ("MMI1.opt_3", "h4.opt_1"), ("h4.opt_2", "h6.opt_1"), ("h6.opt_2", "MMI2.opt_1"),
    ("MMI2.opt_3", "h3.opt_1"), ("h3.opt_2", "SPAR.WG_X"),
    ("MMI2.opt_4", "h2.opt_1"), ("h2.opt_2", "SPAR.WG_Y"),
]
EXTERNAL = ["MMI1.opt_1", "MMI1.opt_2", "SPAR.F1", "SPAR.F2"]


def circuit(f0, v_bit, v_basis, k_phase, loss_db_per_m):
    """4x4 S-matrix, ports [MMI_1.opt_1, MMI_1.opt_2, Fiber mode 1, Fiber mode 2]."""
    h = lambda L, V: heater(L, V, k_phase, loss_db_per_m)  # noqa: E731
    components = {
        "SPAR": psgc(f0), "MMI1": mmi(f0), "MMI2": mmi(f0),
        "h1": h(100e-6, v_basis), "h5": h(200e-6, v_bit),
        "h4": h(100e-6, 0.0), "h6": h(200e-6, 0.0),
        "h2": h(300e-6, 0.0), "h3": h(300e-6, 0.0),
    }
    return solve_network(components, CONNECTIONS, EXTERNAL)


def alice_jones(S):
    """Fiber Jones vector (modes 1, 2) for light launched into MMI_1.opt_1."""
    return S[2:, 0]


def bob_matrix(S):
    """2x2 map: fiber Jones (modes 1, 2) -> amplitudes at [MMI_1.opt_1, MMI_1.opt_2]."""
    return S[:2, 2:]


# ============================================================
# POLARIZATION HELPERS (fiber-mode frame: mode 1 = 'X', mode 2 = 'Y')
# ============================================================
STATES = {
    "H (mode 1)": np.array([1, 0]),
    "V (mode 2)": np.array([0, 1]),
    "D (+45)": np.array([1, 1]) / np.sqrt(2),
    "A (-45)": np.array([1, -1]) / np.sqrt(2),
    "R (+90 deg)": np.array([1, 1j]) / np.sqrt(2),
    "L (-90 deg)": np.array([1, -1j]) / np.sqrt(2),
}


def stokes(E):
    """Normalized Stokes (S1, S2, S3); S3 > 0 when mode 2 leads mode 1 by +90 deg."""
    e1, e2 = E
    s0 = abs(e1) ** 2 + abs(e2) ** 2
    return np.array([abs(e1) ** 2 - abs(e2) ** 2, 2 * np.real(e1 * np.conj(e2)),
                     -2 * np.imag(e1 * np.conj(e2))]) / s0


def rotator(theta):
    """INTERCONNECT polarization rotator in the PSGC fiber-mode frame.

    Fitted from bob_test: it rotates opposite to the mode-1/mode-2 frame, so its
    '-45 deg' setting is this frame's D (+45) and vice versa.
    """
    return np.array([np.cos(theta), -np.sin(theta)])


def db(x):
    return 10 * np.log10(np.maximum(np.abs(x) ** 2, 1e-30))


# ============================================================
# CALIBRATION & DEVICE SUMMARY
# ============================================================
COMBOS = [(0, 0), (0, 1), (1, 0), (1, 1)]  # (basis, bit)


def calibrate(ref=None):
    """Fit heater phase coefficient and loss on Alice's sweep_basis.

    Returns (f0, k_phase, loss_db_per_m, alice_ref); alice_ref is [state, (mode 1, mode 2)] in dB.
    """
    ref = np.load(REF_NPZ) if ref is None else ref
    lam = ref["alice/gain_X/wavelength"].ravel()
    i = np.argmin(abs(lam - 1550e-9))
    f0 = C0 / lam[i]
    gx, gy = ref["alice/gain_X/X gain (dB)"][i], ref["alice/gain_Y/Y gain (dB)"][i]  # [bit, basis]
    alice_ref = np.array([[gx[bit, b], gy[bit, b]] for b, bit in COMBOS])

    k0 = 1.5 * np.pi / (V_HIGH**2 * L_REF)  # nominal: 3pi/2 at 1.2 V, 300 um
    a0 = 0.23 / L_REF  # paper: 0.23 dB per 300 um

    def residual(x):
        model = np.array([db(alice_jones(circuit(f0, V_HIGH * bit, V_HIGH * b, x[0] * k0, x[1] * a0)))
                          for b, bit in COMBOS])
        weight = np.where(alice_ref > -30, 1.0, 0.1)  # deep nulls sit near INTERCONNECT's numerical floor
        return ((model - alice_ref) * weight).ravel()

    fit = least_squares(residual, [1.0, 1.0])
    return f0, fit.x[0] * k0, fit.x[1] * a0, alice_ref


def fiber_rotation(theta):
    """Polarization rotation by theta in the fiber-mode frame (same matrix as Perceval's PR)."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s], [-s, c]])


def device_summary(rot=0.0):
    """Per-basis quantities for the QKD model (basis 0 = rectilinear, 1 = diagonal).

    Alice's states are normalized at her fiber output (mu is defined there), so only
    Bob's transmission and the end-to-end polarization error enter the key rate.
    rot: residual polarization rotation in the fiber (rad), applied as a field, so it
    adds coherently with the chip's own polarization error.
    Returns dict with eta_bob[basis], e_det[basis] (matched-basis error probability,
    averaged over the two bits) and il_alice_db (for reference only).
    """
    f0, k_phase, loss, _ = calibrate()
    R = fiber_rotation(rot)
    eta_bob, e_det, il_alice = np.zeros(2), np.zeros(2), []
    for b, bit in COMBOS:
        E = alice_jones(circuit(f0, V_HIGH * bit, V_HIGH * b, k_phase, loss))
        il_alice.append(-10 * np.log10(np.sum(np.abs(E) ** 2)))
        T = bob_matrix(circuit(f0, 0.0, V_HIGH * b, k_phase, loss))
        out = np.abs(T @ R @ (E / np.linalg.norm(E))) ** 2
        eta_bob[b] += out.sum() / 2
        e_det[b] += out[1 - bit] / out.sum() / 2
    return {"eta_bob": eta_bob, "e_det": e_det, "il_alice_db": np.array(il_alice)}


# ============================================================
# MAIN
# ============================================================
def main():
    ref = np.load(REF_NPZ)
    f0, k_phase, loss, alice_ref = calibrate(ref)
    lam = ref["alice/gain_X/wavelength"].ravel()
    i = np.argmin(abs(lam - 1550e-9))
    ax, ay = ref["alice/angle_X/X angle (rad)"][i], ref["alice/angle_Y/Y angle (rad)"][i]
    combos = COMBOS

    print(f"=== Heater calibration (fitted on Alice, {lam[i]*1e9:.2f} nm) ===")
    print(f"phase at 1.2 V, 300 um : {k_phase * V_HIGH**2 * L_REF / np.pi:.4f} pi  (nominal 1.5 pi)")
    print(f"heater loss            : {loss * L_REF:.3f} dB / 300 um  (paper 0.23)")

    print("\n=== Alice: model vs INTERCONNECT (fiber mode 1 / mode 2, dB; X-Y phase, deg) ===")
    for (b, bit), r in zip(combos, alice_ref):
        E = alice_jones(circuit(f0, V_HIGH * bit, V_HIGH * b, k_phase, loss))
        ph_m = np.degrees(np.angle(E[0] / E[1]))
        ph_r = np.degrees(np.angle(np.exp(1j * (ax[bit, b] - ay[bit, b]))))
        m = db(E)
        print(f"basis {b} bit {bit}: model {m[0]:7.2f} {m[1]:7.2f} {ph_m:7.1f} | ref {r[0]:7.2f} {r[1]:7.2f} {ph_r:7.1f}")

    # ---------- 2. Bob validation (no fitting) ----------
    lb = ref["bob/gain1/wavelength"].ravel()
    j = np.argmin(abs(lb - lam[i]))
    rot, base = ref["bob/gain1/rotation"].ravel(), ref["bob/gain1/base"].ravel()
    g1, g2 = ref["bob/gain1/TE gain (dB)"][j], ref["bob/gain2/TE gain (dB)"][j]
    print("\n=== Bob: model vs INTERCONNECT sweep_base (dB) ===")
    print("(gain1 = MMI_1.opt_2, gain2 = MMI_1.opt_1)")
    print("base  rot   | model opt_1  opt_2 | ref gain2  gain1")
    for nb, vb in enumerate(base):
        T = bob_matrix(circuit(f0, 0.0, vb, k_phase, loss))
        for nr, th in enumerate(rot):
            out = db(T @ rotator(th))
            print(f"{vb:3.1f} {np.degrees(th):6.1f} | {out[0]:11.2f} {out[1]:6.2f} | {g2[nr, nb]:9.2f} {g1[nr, nb]:6.2f}")

    # ---------- 3a. Bob with elliptical / circular light ----------
    print("\n=== Bob prediction for all six cardinal states (fiber-mode frame) ===")
    print("basis  input        | opt_1 (dB) opt_2 (dB) | P(opt_1)")
    for vb, name in [(0.0, "rect"), (V_HIGH, "diag")]:
        T = bob_matrix(circuit(f0, 0.0, vb, k_phase, loss))
        for s, E in STATES.items():
            out = np.abs(T @ E) ** 2
            print(f"{name:5s}  {s:12s} | {db(np.sqrt(out[0])):10.2f} {db(np.sqrt(out[1])):10.2f} | {out[0] / out.sum():.4f}")

    print("\n=== Bob measurement axes (state that fully exits opt_1) ===")
    for vb, name, ideal in [(0.0, "rect", (1, 0, 0)), (V_HIGH, "diag", (0, -1, 0))]:
        T = bob_matrix(circuit(f0, 0.0, vb, k_phase, loss))
        M1 = T[0].conj()[:, None] * T[0][None, :]  # POVM element of opt_1 (unnormalized)
        w, v = np.linalg.eigh(M1)
        s = stokes(v[:, -1])
        tilt = np.degrees(np.arccos(np.clip(np.dot(s, ideal), -1, 1)))
        print(f"{name}: S = ({s[0]:+.4f}, {s[1]:+.4f}, {s[2]:+.4f})  "
              f"off ideal axis {tilt:.2f} deg on Poincare sphere (S3 part {np.degrees(np.arcsin(s[2])):+.2f} deg)")

    # ---------- 3b. Alice -> Bob error matrix (direct fiber link) ----------
    print("\n=== Alice -> Bob, matched bases (direct link, unit fiber) ===")
    print("basis bit | Bob opt_1   opt_2  (dB) | P(error), bit 0 -> opt_1")
    for b, bit in combos:
        E = alice_jones(circuit(f0, V_HIGH * bit, V_HIGH * b, k_phase, loss))
        T = bob_matrix(circuit(f0, 0.0, V_HIGH * b, k_phase, loss))
        out = np.abs(T @ E) ** 2
        p_err = out[1 - bit] / out.sum()
        print(f"  {b}    {bit}  | {db(np.sqrt(out[0])):8.2f} {db(np.sqrt(out[1])):8.2f} | {p_err:.2e}")


if __name__ == "__main__":
    main()
