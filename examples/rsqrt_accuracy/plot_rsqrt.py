#!/usr/bin/env python3
"""Plot rel error of vrsqrt14 + 2 Newton iters as function of r."""
import os
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "out.txt")
PNG  = os.path.join(HERE, "rsqrt_error.png")
PDF  = os.path.join(HERE, "rsqrt_error.pdf")

a, e0, e1, e2 = np.loadtxt(DATA, comments="#", unpack=True)

EPS = np.finfo(np.float64).eps  # 2.22e-16
floor = 0.5 * EPS               # plot zeros there so the log doesn't break

e0_plot = np.where(e0 > 0, e0, floor)
e1_plot = np.where(e1 > 0, e1, floor)
e2_plot = np.where(e2 > 0, e2, floor)

fig, ax = plt.subplots(figsize=(9, 5.5))

ax.loglog(a, e0_plot, lw=1.0, color="#1f77b4",
          label=r"$y_0 = \mathrm{vrsqrt14pd}(r)$")
ax.loglog(a, e1_plot, lw=1.0, color="#2ca02c",
          label=r"$y_1$ after 1 Newton iter")
ax.loglog(a, e2_plot, lw=1.0, color="#d62728",
          label=r"$y_2$ after 2 Newton iters")

ax.axhline(2.0**-14, color="#1f77b4", ls=":", alpha=0.5,
           label=r"$2^{-14}$")
ax.axhline(EPS, color="black", ls=":", alpha=0.5,
           label=r"$\varepsilon_{\mathrm{double}} \approx 2.22\!\times\!10^{-16}$")

ax.set_xlabel(r"$r$")
ax.set_ylabel(r"$|\,y/(1/\sqrt{r}) - 1\,|$")

ax.set_xlim(a.min(), a.max())
ax.set_ylim(0.1 * EPS, 1e-3)
ax.grid(which="both", ls=":", alpha=0.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper center", frameon=False, ncol=1, fontsize=9,
          bbox_to_anchor=(0.5, -0.13))

fig.tight_layout()
fig.savefig(PNG, dpi=140, bbox_inches="tight")
fig.savefig(PDF, bbox_inches="tight")
print(f"wrote {PNG}")
print(f"wrote {PDF}")

for label, arr in (("y0", e0), ("y1", e1), ("y2", e2)):
    nz = arr[arr > 0]
    print(f"{label}: max rel err = {arr.max():.3e},  max nonzero = {nz.max() if nz.size else 0:.3e}")
