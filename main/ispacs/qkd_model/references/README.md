# References used by the QKD model

arXiv preprints of the papers whose formulas are implemented in `../qkd_model.py`
and checked in `../benchmark_published.py`.

| File | Paper | Used for |
|---|---|---|
| `GLLP2004_security_imperfect_devices_arXiv0212066.pdf` | D. Gottesman, H.-K. Lo, N. Lütkenhaus, J. Preskill, "Security of quantum key distribution with imperfect devices", *Quantum Inf. Comput.* 4, 325 (2004) | Key rate with multi-photon pulses (`no_decoy_rate`: every multi-photon detection is treated as known to Eve) |
| `LoMaChen2005_decoy_state_QKD_arXiv0411004.pdf` | H.-K. Lo, X. Ma, K. Chen, "Decoy state quantum key distribution", *Phys. Rev. Lett.* 94, 230504 (2005) | The decoy-state idea and the asymptotic key rate with known single-photon yield (`infinite_decoy_rate`) |
| `Ma2005_practical_decoy_state_arXiv0503005.pdf` | X. Ma, B. Qi, Y. Zhao, H.-K. Lo, "Practical decoy state for quantum key distribution", *Phys. Rev. A* 72, 012326 (2005) | Fiber/detector model (Sec. III); vacuum + weak decoy bounds; **benchmark**: GYS parameters (Table I), 142.05 km and 140.55 km maximum distances (Fig. 3) |
| `Lim2014_finite_key_decoy_BB84_arXiv1311.7129.pdf` | C. C. W. Lim, M. Curty, N. Walenta, F. Xu, H. Zbinden, "Concise security bounds for practical decoy-state quantum key distribution", *Phys. Rev. A* 89, 022307 (2014) | Finite-key, 3-intensity decoy BB84, Eqs. (1)-(5) (`key_length`); **benchmark**: Fig. 1 (key rate vs distance for block sizes 10^4-10^9) |
