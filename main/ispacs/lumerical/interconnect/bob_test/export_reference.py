"""Export INTERCONNECT reference data used to validate bob_model.py.

Opens the projects hidden through lumapi and reads the *saved* sweep results with
getsweepresult. Nothing is simulated and the projects are never saved.

    Alice : alice_test/alice_test_freq.icp :: sweep_basis  (fiber X/Y gain + angle)
    Bob   : bob_test/bob_test_freq.icp     :: sweep_base   (out gain1/gain2 vs rotation)

Output: interconnect_reference.npz next to this file.
Requires Lumerical v241 (Python API) and a license seat while it runs (~1 min).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.append(r"C:\Program Files\Lumerical\v241\api\python")
import lumapi  # noqa: E402

HERE = Path(__file__).resolve().parent
INTC = HERE.parent

SOURCES = {
    "alice": (INTC / "alice_test" / "alice_test_freq.icp", "sweep_basis"),
    "bob": (INTC / "bob_test" / "bob_test_freq.icp", "sweep_base"),
}

data = {}
for tag, (project, sweep) in SOURCES.items():
    with lumapi.INTERCONNECT(filename=str(project), hide=True) as ic:
        for result in ic.getsweepresult(sweep).split("\n"):
            result = result.strip()
            for key, value in ic.getsweepresult(sweep, result).items():
                if isinstance(value, np.ndarray) and key not in ("null", "index_1"):
                    data[f"{tag}/{result}/{key}"] = value

out = HERE / "interconnect_reference.npz"
np.savez(out, **data)
print(f"saved {len(data)} arrays to {out}")
