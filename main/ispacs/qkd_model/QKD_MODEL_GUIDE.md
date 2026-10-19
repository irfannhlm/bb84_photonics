# How the BB84 key-rate model works

A guide for photonics engineers: what the model computes, why, and how each step
was checked. If you are used to INTERCONNECT/Tidy3D (powers, S-parameters, ER, IL),
the good news is that **most of the model is classical photonics**. The "quantum" part
comes down to two ideas: light arrives as discrete photons, and the security proof
decides how much of the raw key Eve could know.

---

## 1. The big picture

```
 Your chip (classical)          Photons & detectors            Security proof
 ─────────────────────          ───────────────────            ──────────────
 Alice circuit  ┐                                              
 PSGC S-matrix  ├─► Jones vectors ─► click probabilities ─► QBER, gains ─► secret key rate (SKR)
 Bob circuit    ┘   (bob_model.py)   (qkd_model.py)          (qkd_model.py, formulas from papers)
                    validated vs      checked with Perceval    checked against published curves
                    INTERCONNECT
```

The question the model answers: **for each fiber length, how many secret bits per
second can Alice and Bob produce, given this chip, typical detectors and a security
proof?** The answer is the SKR-vs-distance curve.

---

## 2. Glossary

Grouped from "optics you know" to "security you don't (yet)".

### Light and polarization

| Term | Meaning |
|---|---|
| **Jones vector** | Two complex numbers giving the field amplitude and phase in two orthogonal polarizations (here: the PSGC's two fiber modes, "X" and "Y"). Exactly what INTERCONNECT's ONA gain + angle give you. |
| **Azimuth, ellipticity** | The angle of the polarization ellipse, and how "round" it is. Ellipticity 0 means linear polarization. |
| **Stokes vector, Poincaré sphere** | Another way to write a polarization: a point on a sphere. H and V are opposite points, D and A are opposite points 90° away, circular states sit at the poles. A polarization rotation in fiber is a rotation of this sphere. |
| **H, V, D, A** | Horizontal (0°), vertical (90°), diagonal (+45°), anti-diagonal (−45°). The four BB84 states. |
| **Basis** | A pair of orthogonal states used together: rectilinear {H, V} (called **Z** in the code) or diagonal {D, A} (called **X**). |
| **Transmittance η** | Fraction of power (= fraction of photons) that survives a stage: `10^(−loss_dB/10)`. |

### Photons

| Term | Meaning |
|---|---|
| **Photon / Fock state \|n⟩** | A pulse containing exactly n photons. |
| **Weak coherent pulse (WCP)** | An attenuated laser pulse. It does not contain a fixed number of photons: the number is random with a **Poisson distribution**. |
| **μ (mu)** | Mean photon number per pulse at Alice's output (e.g. 0.5). With μ = 0.5: 61% of pulses are empty, 30% have 1 photon, 9% have 2 or more. |
| **Phase randomization** | Each laser pulse has a random optical phase. Then the pulse is mathematically a *mixture* of Fock states with Poisson weights. Gain-switched lasers do this naturally. Every security proof here assumes it. |
| **Multi-photon pulse** | A pulse with ≥ 2 photons. Dangerous: see PNS below. |

### Detection

| Term | Meaning |
|---|---|
| **Threshold detector** | Can only say "click" (≥ 1 photon) or "no click". It cannot count photons. Single-photon detectors work this way. |
| **SPAD** | Single-photon avalanche diode (InGaAs for 1550 nm). **Gated**: it is armed only in a short window around each expected pulse. |
| **Detection efficiency η_d** | Probability that a photon arriving at the detector causes a click (20% for the SPAD here). |
| **Dark count** | A click with no photon (thermal noise). Probability per gate, here 10⁻⁶. |
| **Afterpulse** | A false click caused by charge trapped from a previous click. Here 1% per click. |
| **Dead time** | After a click the detector is switched off for a while (10 µs here) to suppress afterpulses. Limits the maximum click rate. |
| **Double click** | Both detectors click at once. The protocol assigns a random bit (required by the security proofs). |
| **SNSPD** | Superconducting nanowire detector: ~85% efficiency, very low noise, but needs cryogenic cooling. Used only as a comparison. |

### Protocol and security

| Term | Meaning |
|---|---|
| **BB84** | Alice sends one of four states (random basis, random bit). Bob measures in a random basis. They then publicly compare *bases* (not bits). |
| **Sifting** | Keeping only the events where Alice's and Bob's bases matched. With mismatched bases Bob's result is random (50%), so those bits are discarded. |
| **Efficient BB84 / basis bias p_Z** | Choosing one basis most of the time (e.g. 90%) so most events survive sifting; the other basis is used to test for Eve. |
| **Gain Q** | Probability that a sent pulse produces a click at Bob. |
| **QBER** | Quantum bit error rate: the fraction of sifted bits where Bob's bit ≠ Alice's bit. |
| **Error correction (EC)** | Alice and Bob fix their differing bits by exchanging parity information over a public channel. That information is visible to Eve and costs `f_EC × h(QBER)` bits per bit; f_EC = 1.16 means 16% worse than the theoretical minimum. |
| **Binary entropy h(x)** | `h(x) = −x log₂x − (1−x) log₂(1−x)`. "How many bits of uncertainty" a coin with probability x carries. h(0) = 0, h(0.11) ≈ 0.5, h(0.5) = 1. |
| **Privacy amplification (PA)** | Compressing the corrected key (with a hash function) so that what Eve knows is erased. The compression needed is set by the **phase error rate**. |
| **Phase error rate e₁ (φ)** | The error rate Eve would have caused on single photons, measured in the *other* basis. It bounds Eve's information. For ideal single photons, PA costs `h(e₁)`. |
| **Photon-number-splitting (PNS) attack** | On a multi-photon pulse, Eve keeps one photon and forwards the rest. She learns the bit without causing any error. So multi-photon detections must be treated as fully known to Eve. |
| **Single-photon yield Y₁** | Probability that a 1-photon pulse is detected. Only those detections can produce secret key. |
| **Decoy states** | Alice randomly varies μ (signal, decoy, near-vacuum). Eve cannot tell the pulses apart, so comparing click rates at different μ lets Alice and Bob *measure* Y₁ and e₁ instead of assuming the worst. Needs an intensity modulator. |
| **GLLP** | Gottesman–Lo–Lütkenhaus–Preskill (2004): the key-rate formula for imperfect sources. Key comes only from single photons; multi-photon events are tagged as known to Eve. |
| **Asymptotic** | The idealization of infinitely many pulses, so every rate is known exactly. |
| **Finite key** | A real key is made from a finite **block** of N pulses, so every rate is estimated with statistical uncertainty; the proof must assume the worst value inside the confidence interval. |
| **ε_sec, ε_cor** | Allowed failure probabilities: ε_sec = 10⁻¹⁰ (Eve learns something), ε_cor = 10⁻¹⁵ (Alice's and Bob's keys differ). Smaller ε means a more conservative key. |
| **SKR** | Secret key rate: secret bits per second (or per pulse) after sifting, EC and PA. |

---

## 3. Step by step

### Step 1: Classical circuit → polarization states  (`bob_model.py`)

This is your territory. The circuit (MMIs, heaters, PSGC) is a linear optical network,
so I rebuilt it in Python from the **same S-parameter files INTERCONNECT uses**
(PSGC from Tidy3D, MMI from your CML) and solved the full network including
reflections.

- Only two numbers had to be fitted: the heater phase per V²·length (fit: 3π/2 × 1.006
  at 1.2 V, 300 µm) and the heater loss (0.15 dB per 300 µm). They were fitted on
  Alice's sweep only.
- **Validation:** with no further fitting, all 8 of Bob's INTERCONNECT sweep points
  match (bright outputs within 0.1 dB, ~40 dB nulls within ~1 dB).

**The key link between classical and quantum:** for a linear circuit, a single photon
follows the *same* transfer matrix as the classical field. The probability that the
photon exits a port equals the fraction of classical power that exits that port. So
from your circuit we get, for each basis:

| Quantity | Meaning | Value |
|---|---|---|
| η_Bob | Fraction of photons that pass Bob's chip | 0.356 (−4.48 dB) |
| e_det (rectilinear, Z) | Probability a photon lands in the *wrong* detector | 1.2 × 10⁻⁴ |
| e_det (diagonal, X) | same, diagonal basis | 7.2 × 10⁻⁴ |

e_det is computed from the full Jones vectors, Alice's actual output state projected
onto Bob's actual measurement. That is why the paper's ER and polarization-angle error
should not be added together: for rectilinear states they are the *same* imperfection.

Alice's chip loss (5.5 dB) does **not** enter: μ is defined at Alice's output, so she
simply attenuates less. Her loss only costs laser power.

### Step 2: Photons through the link  (`optical_wcp`, `optical_single`)

Total transmittance from Alice's output to a click:

  η = 10^(−0.2·L/10) × η_Bob × η_d

A pulse with mean μ photons delivers on average μη photons to Bob's detectors. Because
the photon number is Poisson, the number reaching each detector is again Poisson, so
the probability that the correct detector fires is `1 − exp(−μη(1−e_det))` and that
the wrong one fires is `1 − exp(−μη·e_det)`. These two are independent for Poisson
light (not for a single photon, which can only go one way; that case has its own
formula).

### Step 3: Detector noise  (`add_noise`, `detector_state`)

On top of the photon clicks, each detector clicks by itself with probability y₀ per gate:

- y₀ = dark count + afterpulses, where the afterpulse part is proportional to the
  recent click rate;
- dead time removes a fraction of gates: availability `A = 1 / (1 + click rate × dead time)`.
  At 0 km the click rate is so high that the detectors are mostly "dead": this is why the
  SKR at 0 km is limited to ~120 kbps, not by the chip;
- noise clicks carry a random bit (50% wrong), and double clicks get a random bit.

Result: for every μ and distance, the **gain Q** (clicks per pulse) and **QBER E**.

### Step 4: From QBER to secret key  (the security formulas)

Every key-rate formula has the same structure:

```
 secret key  =  (single-photon detections) × [1 − h(e₁)]   −   (all detections) × f_EC × h(E)
                └──────── what survives privacy amplification ┘   └── cost of error correction ──┘
```

- The first term is the only source of key: single-photon detections, minus what
  privacy amplification removes (set by the phase error e₁).
- The second term is paid on *every* detection, because error correction reveals
  parity bits for the whole sifted key.
- Key exists only while the first term is larger than the second.

The three protocol cases differ in **how well Y₁ and e₁ are known**:

| Case | How Y₁ and e₁ are obtained | Function | Hardware needed |
|---|---|---|---|
| **No decoy** (GLLP) | Worst case: assume every multi-photon pulse that *could* be detected was, and Eve knows it all. μ must be tiny (~η), so the rate scales like η². | `no_decoy_rate` | Current chip |
| **Ideal decoy** (asymptotic bound) | Assume Y₁ and e₁ are known exactly (infinitely many decoys, infinitely long key). An upper bound. | `infinite_decoy_rate` | — |
| **Decoy, finite key** (Lim et al. 2014) | 3 intensities; Y₁ and e₁ are *bounded* from the click statistics of a finite block, with confidence ε. Rate scales like η. | `key_length`, `finite_key_rate` | Intensity modulator |

For the finite-key case, the free settings (two intensities, how often each is sent,
and the basis bias p_Z) are optimized at every distance (`maximize_protocol`): a random
search followed by local refinement. Fixed hand-picked values would understate the system.

### Step 5: Polarization drift  (`load_device(misalign_deg)`)

Fiber birefringence slowly rotates the polarization (a rotation on the Poincaré sphere).
A polarization-encoded link needs active compensation, so the model only asks what the
**residual** misalignment costs. The rotation is applied to the field between Alice and
Bob, so it adds to the chip error as an amplitude, not as a probability. Result: 5° of
residual misalignment costs ~4 km (158 → 154 km), because at long distance dark counts,
not polarization errors, limit the range.

---

## 4. How each part was verified

| Part | Checked against | Result |
|---|---|---|
| Chip model (`bob_model.py`) | INTERCONNECT sweeps (Alice: 4 states; Bob: 8 states, both bases) | Bright outputs within 0.1 dB, ~40 dB nulls within ~1 dB |
| Optics + photon statistics + sifting (`qkd_model.py`) | **Perceval** (Quandela), built from its own components: `Unitary`, `LC` loss, `PBS`, `HWP`, `PR`; Fock-state inputs n = 1…8 | SKR within 0.005%, QBER within 0.04%; reproduces textbook 50% (mismatched bases) and 25% (intercept-resend attack) |
| Security formulas | Published curves, digitized from the papers' own figure files (`benchmark_data/`); see section 5 | see section 5 |
| Detector noise model | Standard model (same as Ma 2005 / Lim 2014) | Not independently simulated: Perceval has no dark-count model |

**What Perceval does and does not check.** Perceval receives only the list of optical
components and computes photon-number probabilities with its own simulator. None of
the formulas in Steps 2–3 is given to it. It cannot check the security formulas
(Step 4); those are pure maths from the papers, which is why they are benchmarked
against published results instead.

---

## 5. Benchmarks against published work

The security formulas cannot be checked by a physics simulator, so they were run with
the parameters of the original papers and compared with the papers' own results.
The curves were recovered directly from the vector figure files in the arXiv sources
(`benchmark_data/digitize_eps.py`); the digitization reproduces every endpoint stated in
the papers' text (e.g. 141.9 vs 142.05 km, 134.6 vs "135 km"). Script:
`benchmark_published.py`, figure: `figures/fig_benchmark_published.png`.

**[A] Ma, Qi, Zhao, Lo 2005 (GYS experiment, asymptotic decoy).** Fully specified, and
the cleanest check of the decoy formulas and the detector model:

| Quantity | This model | Paper |
|---|---|---|
| Max distance, ideal decoy | 142.02 km | 142.05 km |
| Max distance, vacuum + weak decoy | 140.62 km | 140.55 km |
| Rate along the curves | within 0.1–3% | — |
| Optimal μ (f = 1.22 / f = 1) | 0.482 / 0.547 | 0.48 / 0.54 |

**[B] Lo, Ma, Chen 2005 (no decoy, GLLP).** The formula is the same as the paper's
equation, term by term. The paper does not state μ for its no-decoy curve (~32 km). With μ
optimized at every distance (as in this model) the same formula gives 43 km; a fixed
μ ≈ 0.012 with standard sifting reproduces the published curve (30.5 km). So this one is a
consistency check, not a strict benchmark, and it is left out of the figure.

**[C] Lim et al. 2014 (finite key, the formula used for the talk curves).**

| Key block n_X | Reach (rate falls to 10⁻⁷) this model / paper | Rate vs paper |
|---|---|---|
| 10⁹ | 179.8 / 179.6 km | +4 to +7% |
| 10⁸ | ≥175 / 177.1 km (5 km grid) | +3 to +5% |
| 10⁷ | 172.0 / 173.2 km | −2 to +3% |
| 10⁶ | ≥165 / 167.5 km (5 km grid) | −25 to −4% |
| 10⁵ | 150.6 / 156.1 km | −64 to −17% |
| 10⁴ | 118.6 / 134.6 km | −53% at 0 km, worse near the end |
| asymptotic (n_X → ∞) | 185.1 / 191.9 km | −2% at short distance, −60% near the end |

Their text states that at 100 km the rate for n_X = 10⁹ is "about 1.75×" that for 10⁷; the
model gives 1.79×.

For blocks of 10⁷ and above, the model reproduces the paper within a few percent. For
small blocks it gives *less* key than the published figure. This is not the optimizer
(three independent searches agree to all digits) and not a simple variant of the
statistical term (tested); the exact cause is unknown, possibly a detail of how the
published figure was generated. **What matters for the talk:** the main curve (17-minute
block) uses key blocks of 10⁸–10⁷ out to ~120 km, where the benchmark agrees within a
few percent; near its end (150–158 km) the block is ~10⁶, where this model is 4–25%
below the published curve. Every disagreement is on the conservative side: the true range
can only be equal or slightly longer, never shorter.

---

## 6. Results for this chip (talk numbers)

InGaAs SPAD (η = 20%, dark 10⁻⁶/gate, afterpulse 1%, dead time 10 µs), 100 MHz,
0.2 dB/km, f_EC = 1.16, ε_sec = 10⁻¹⁰, Bob chip loss 4.48 dB.

| Case | Max distance | SKR at 0 km | SKR at 100 km |
|---|---|---|---|
| No decoy (current chip) | 58 km | 58 kbps | 0 |
| Decoy, finite key, 10 s block (N = 10⁹) | 122 km | 86 kbps | 1.6 kbps |
| Decoy, finite key, 17 min block (N = 10¹¹) | **158 km** | 118 kbps | 7.2 kbps |
| Ideal decoy (asymptotic bound) | 172 km | 142 kbps | 16.5 kbps |
| Decoy, 17 min, SNSPD instead of SPAD | 226 km | — | — |

Other findings:

- The chip's polarization error contributes **0.012%** QBER (key basis), one to two orders
  of magnitude below detector noise. The detectors, not the chip, limit the range.
- Equivalent-ER sensitivity: even a 20 dB circuit would lose only ~6 km of range. High ER
  mainly helps at short and metro distances (lower QBER, cheaper error correction).
- An optimized PSGC (~70% coupling instead of 38%) would roughly double the SKR at 100 km
  and add ~12 km of range. Only Bob's coupling matters for the key rate.
- QBER never reaches the 11% "BB84 threshold". 11% is the limit for ideal single photons
  with perfect error correction; with laser pulses the key ends earlier, when dark counts
  swamp the single-photon signal (here at 5–7% QBER).

---

## 7. Likely questions in the talk

**"Is this a real quantum simulation?"**
The optics are simulated at the level of individual photons (Poisson photon-number
statistics, threshold detectors, multi-photon events), and this part was cross-checked
with Perceval. The key rate then follows from published security proofs (GLLP, decoy
state, Lim et al. finite-key), benchmarked against the original papers.

**"Why does your paper say 137 km?"**
The paper's estimate used the ideal-decoy formula with a fixed μ = 0.2, counted Alice's
loss as channel loss, and ignored detector dead time. The full model gives 158 km for a
practical decoy protocol with a 17-minute key block, and 58 km without decoys.

**"Why do you need decoy states?"**
Laser pulses sometimes contain two or more photons; Eve can split one off without being
noticed (PNS attack). Decoy states let Alice and Bob measure how many detections came
from true single photons. Hardware: an intensity modulator before the MZM (it can use
the same segmented-shifter idea to make three levels from NRZ drives).

**"What about polarization drift in the fiber?"**
Like every polarization-encoded system, the link needs active polarization compensation.
With residual misalignment up to 5°, the range drops by only ~4 km.

**"Why is the chip error so small but the range still limited?"**
Range is set by loss and detector noise: at long distance, dark counts become comparable
to real photon clicks. The chip contributes 0.012% QBER; the detectors contribute
several percent near the maximum distance.

---

## 8. FAQ: concepts behind the model

### The eavesdropper

**How does Eve actually eavesdrop? Does she have to cut the fiber?**
She needs physical access, but that is easy: a fiber can be cut and spliced in minutes
(it looks like a maintenance outage), and classical fibers can even be tapped by bending
them. Security proofs do not argue about access at all: they assume Eve **controls the
entire channel** and has any technology allowed by physics. Main attack types:

- *Intercept-resend*: measure each photon in a random basis, send a new one. Wrong basis
  half the time gives **25% QBER**, so it is easily caught.
- *Beam splitting*: tap off part of the light. Bob only sees extra loss. A single photon
  cannot be split, but a 2-photon pulse can, giving Eve a copy with no errors.
- *Photon-number splitting (PNS)*: Eve measures the photon number without disturbing the
  polarization, keeps one photon from every multi-photon pulse and blocks single photons
  so the total loss looks normal. Not possible with today's technology, but allowed in
  the proofs.
- *Device attacks ("quantum hacking")*: e.g. blinding Bob's detectors with bright light.
  These attack the hardware, not the protocol, and are the most practical real threat.

**Why is all fiber loss attributed to Eve?**
Bob cannot distinguish "lost in the fiber" from "taken by Eve", so the proof assumes the
worst. This is also why Bob's chip loss counts against the key rate.

**Why can't Eve just copy the photons?**
The *no-cloning theorem*: an unknown quantum state cannot be copied. Any measurement in the
wrong basis disturbs the state and creates errors, so the QBER bounds what Eve can know.

**Why can't the signal be amplified to go further?**
An amplifier would be a copying machine, which no-cloning forbids (it adds noise that
destroys the encoding). So the range is set by loss and detector noise. Going much further
needs trusted relay nodes, satellites, or future quantum repeaters.

### Photons

**Do quantum systems always use exactly one photon? How do they know the number?**
Most practical QKD, including this design, uses **attenuated laser pulses**. Nobody knows
the photon number of a given pulse; only the **average** μ is known, set from power,
photon energy and repetition rate:

```
photon energy at 1550 nm:  E = hc/λ = 1.28×10⁻¹⁹ J
μ = P / (E × f_rep)   →   μ = 0.5 at 100 MHz is P = 6.4 pW (−82 dBm) at Alice's output
```

In practice: measure strong light with a power meter, then add a calibrated attenuator.
True single-photon sources exist (quantum dots; "heralded" sources where detecting one
photon of a pair announces the other) and are used in photonic quantum computing, but they
are harder to build. Bob's detectors cannot count photons either; they only click.

**What is phase randomization, and does my laser need it?**
Each laser pulse must have a random optical phase. Then a pulse is mathematically a
*mixture* of definite photon numbers with Poisson weights, which is what every proof here
assumes. A gain-switched (directly pulsed) laser randomizes the phase naturally; a CW laser
carved into pulses would need an extra phase modulator driven by random numbers.

### Decoy states

**What exactly is uncertain without decoys, and what does "asymptotic" mean?**
The uncertain quantities are the **single-photon yield Y₁** (chance that a 1-photon pulse
is detected) and the **single-photon error e₁**. They are never measured directly; they
are inferred. Two separate things limit the inference:

1. **Too few intensities.** Each intensity gives one equation, but there are infinitely many
   unknowns (Y₀, Y₁, Y₂, ...). With 3 intensities you only get *bounds*. With **infinitely
   many decoys** you could solve exactly.
2. **Finite data.** Measured click rates fluctuate around the true values, so confidence
   intervals are needed. With **infinitely many pulses** ("asymptotic", or infinite key)
   the fluctuations vanish.

The "ideal decoy (asymptotic bound)" curve removes both. Ma et al. show that with infinite
data but only 3 intensities, the range drops only from 142.05 to 140.55 km, so most of the
gap between the ideal curve (172 km) and the practical one (158 km) is the finite data.

**Which "rates" do decoys compare?**
The **gain**: probability that a pulse of a given intensity produces a click (clicks per
pulse), measured separately for signal, decoy and near-vacuum pulses. The QBER per
intensity is also used.

**How do decoys increase the key rate? Do they create more single photons?**
No. Decoys do not change how many single photons are sent, and Bob still cannot tell single
from multi-photon pulses one by one. What they give is a **statistical measurement** of how
many detections came from single photons:

- Eve sees only the photon number of a pulse, not which intensity Alice chose, so a
  1-photon pulse gets through with the same probability Y₁ whether it was "signal" or "decoy".
- Each intensity mixes photon numbers in a different, known (Poisson) proportion; weak
  decoys contain relatively more single photons.
- Comparing the gains at the three intensities pins down Y₁ and e₁. A PNS attack (blocking
  single photons, passing multi-photon ones) would make the gains inconsistent and be detected.

Without decoys you must assume Eve got every multi-photon pulse. To keep that small, μ must
be tiny (about the channel transmission), so few single photons arrive and the rate scales
like loss². With decoys, μ ≈ 0.5 is allowed and the rate scales like loss. That is the
difference between 58 km and 158 km.

### Key-rate bookkeeping

**What does "cost" mean in the formulas?**
Bits of key that have to be given up:

- *Error correction*: `f_EC × h(E)` per sifted bit, the parity bits exchanged publicly to fix
  errors (Eve sees them).
- *Privacy amplification*: `h(e₁)` per single-photon bit, how much the key is compressed
  by hashing to erase what Eve might know.
- *Finite-key*: statistical uncertainty is paid in bits too. The proof uses the worst values
  in the confidence interval (fewer guaranteed single photons, higher phase-error bound),
  plus a fixed overhead of about 280 bits per block for ε_sec and ε_cor.

**Is the key block n_X the number of pulses per block?**
No. n_X is the number of **detections in the key basis** that survive sifting (the raw key
length processed at once). The number of pulses N is much larger:
`N ≈ n_X / (P(both choose the key basis) × P(Bob clicks))`. For the talk model, N = 10¹¹ pulses
(17 minutes at 100 MHz) gives n_Z ≈ 1.5×10⁸ at 0 km and ≈ 6×10⁵ at 156 km. Lim et al. fix
n_X (post-processing memory); this model fixes N (acquisition time). Same physics,
different convention.

**Why is Alice's chip loss irrelevant but Bob's is not?**
μ is defined at Alice's output, and Alice's lab is trusted: she just attenuates less.
Everything after her output (fiber and Bob's chip) is treated as Eve-accessible loss.

**What does "secure" mean in numbers (ε)?**
ε_sec = 10⁻¹⁰ means: the probability that the final key is distinguishable from a perfectly
secret one (i.e. Eve learned anything) is at most 10⁻¹⁰. ε_cor = 10⁻¹⁵ bounds the chance that
Alice's and Bob's keys differ. Smaller values cost a little more key.

**Does QKD work on its own?**
No. Alice and Bob also need a **classical channel** for sifting, error correction and
privacy amplification, and that channel must be **authenticated** (so Eve cannot pose as
Bob). Authentication needs a small pre-shared secret key, so QKD is strictly a "key
growing" scheme. The final key is then used for encryption (e.g. a one-time pad or AES).

### Your design

**Why must the link have polarization compensation?**
Fiber birefringence rotates polarization randomly and slowly (temperature, stress). The PSGC
measures in fixed axes, so without compensation the error would average around 25% or more.
Every polarization-encoded system uses a feedback polarization controller (driven by
reference pulses or by the measured QBER). The model asks only what the *residual*
misalignment costs: 5° → about 4 km of range.

**Why polarization encoding instead of time-bin encoding?**
Time-bin encoding is immune to fiber polarization drift but needs unbalanced
interferometers (long delay lines) on chip, which are large and thermally sensitive.
Polarization encoding with a 2D PSGC is compact, at the price of needing drift compensation.
This trade-off is the core argument of your architecture.

**What do "MDI-QKD" and "twin-field QKD" mean (if someone asks about future work)?**
*Measurement-device-independent* QKD moves the detectors to an untrusted middle node, which
removes all detector attacks. *Twin-field* QKD also uses a middle node and scales like
√loss, reaching 500–1000 km in experiments. Both need interference between Alice's and
Bob's pulses, so they are more demanding than BB84.

---

## 9. Assumptions and limitations

- **Single wavelength** (1550 nm). The laser linewidth (~0.01 nm for ns pulses) is far
  narrower than the PSGC bandwidth.
- **Detector values are typical, not a specific product.** Results scale strongly with
  them (see the SNSPD case).
- **The PSGC S-matrix data is not reciprocal** (WG→fiber 30.6%, fiber→WG 38.4%). It is
  used as-is, as INTERCONNECT does. Physically both directions should be equal; this is
  most of the reason Alice's IL (5.55 dB) and Bob's (4.48 dB) differ. Worth checking in
  the Tidy3D setup.
- **Decoy curves assume an intensity modulator** that the current chip does not have.
- The **no-decoy** curve is asymptotic (infinite key), so it is an optimistic version of
  what the current chip achieves.
- Detector dead time is modelled as an average availability; detectors are assumed
  identical (no efficiency mismatch beyond the chip's ~0.1 dB basis-dependent loss).

---

## 10. Files

| File | What it does |
|---|---|
| `../lumerical/interconnect/bob_test/bob_model.py` | Circuit model (S-parameter network), heater calibration, INTERCONNECT validation |
| `../lumerical/interconnect/bob_test/export_reference.py` | Reads INTERCONNECT sweep results (read-only, needs Lumerical) |
| `qkd_model.py` | Photon/detector model, security formulas, optimizer |
| `make_figures.py` | Talk figures `figures/fig_qkd_main.png`, `figures/fig_qkd_range.png` (`--rerun` after model changes) |
| `perceval_check.py` | Perceval cross-check, `figures/fig_perceval_check.png` |
| `benchmark_published.py` | Benchmarks vs Ma 2005, Lo-Ma-Chen 2005, Lim 2014, `figures/fig_benchmark_published.png` |
| `benchmark_data/` | Digitized published curves + the digitizer |
| `references/` | PDFs of the papers (see its README) |
