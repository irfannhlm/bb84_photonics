"""Recover curve data from MATLAB-generated EPS figures (arXiv source files).

MATLAB writes each line as relative 'rlineto' steps followed by "x0 y0 n MP"; the line
style (SO solid, DA dashed, DO dotted, DD dash-dot) is set before it. Axis calibration
comes from the tick marks in the same file (see CALIB below).

Usage (sources from https://arxiv.org/e-print/1311.7129 and .../quant-ph/0503005):
    python digitize_eps.py Figure1_D1.eps lim2014
    python digitize_eps.py ComWang.eps ma2005
    python digitize_eps.py Decoy.eps lmc2005      (https://arxiv.org/e-print/quant-ph/0411004)
Writes <name>_digitized.csv next to this script: style, curve index, distance_km, rate.
"""
import re
import sys
from pathlib import Path

import numpy as np

# (x pixel at x0, x0 km, x pixel at x1, x1 km, y pixel of 1e-7, pixels per decade)
CALIB = {  # decade spacing from the major ticks 1e-7 ... 1e-2
    "lim2014": (1306, 0, 13266, 200, 7874, (7874 - 493) / 5),  # Lim et al. 2014, Fig. 1
    "ma2005": (899, 0, 6255, 160, 4000, (4000 - 389) / 5),  # Ma et al. 2005, Fig. 3
    "lmc2005": (899, 0, 6255, 160, 4133, (4133 - 389) / 5),  # Lo, Ma, Chen 2005, Fig. 1
}


def polylines(eps_text):
    """Yield (style, points) for every MP polyline inside the plot clip region."""
    # data lines follow the plot clip ("gs x y w h rc" or "gs x y w h MR c np") set after the tick labels
    clip = re.compile(r"gs\s+\d+\s+\d+\s+\d+\s+\d+\s+(rc|MR c np)")
    body = eps_text[clip.search(eps_text, eps_text.index("(0) s")).end():]
    stack, style = [], "SO"
    for t in body.split():
        if re.fullmatch(r"-?\d+(\.\d+)?", t):
            stack.append(float(t))
            continue
        if t in ("SO", "DA", "DO", "DD"):
            style = t
        if t == "MP" and len(stack) >= 3:
            n = int(stack[-1])
            x0, y0 = stack[-3], stack[-2]
            steps = np.array(stack[-3 - 2 * (n - 1):-3]).reshape(-1, 2)[::-1]
            yield style, np.vstack([[x0, y0], [x0, y0] + np.cumsum(steps, axis=0)])
        stack = []


def digitize(path, name):
    xp0, xv0, xp1, xv1, y7, dec = CALIB[name]
    rows, seen = [], []
    for style, p in polylines(Path(path).read_text(encoding="latin-1")):
        km = xv0 + (p[:, 0] - xp0) / (xp1 - xp0) * (xv1 - xv0)
        rate = 10 ** (-7 + (y7 - p[:, 1]) / dec)
        key = (style, round(km[0], 3), round(rate[0], 12), len(km))
        if key in seen:  # MATLAB sometimes draws the same line twice
            continue
        seen.append(key)
        rows += [(style, len(seen) - 1, k, r) for k, r in zip(km, rate)]
    out = Path(__file__).with_name(f"{name}_digitized.csv")
    with open(out, "w") as f:
        f.write(f"# digitized from {Path(path).name} ({name}); style SO/DA/DO/DD = MATLAB line style\n")
        f.write("style,segment,distance_km,rate_per_pulse\n")
        f.writelines(f"{s},{i},{k:.4f},{r:.6e}\n" for s, i, k, r in rows)
    return out


if __name__ == "__main__":
    print("wrote", digitize(sys.argv[1], sys.argv[2]))
