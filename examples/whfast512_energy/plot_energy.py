
import sys, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base_dir = sys.argv[1] if len(sys.argv) > 1 else "run_baseline"
fix_dir  = sys.argv[2] if len(sys.argv) > 2 else "run_fixed"
outpng   = sys.argv[3] if len(sys.argv) > 3 else "energy_vs_time.png"

def load(d):
    runs = []
    for f in sorted(glob.glob(os.path.join(d, "seed*.txt"))):
        if os.path.basename(f) == "seed0.txt":      
            continue
        try:
            a = np.loadtxt(f, comments="#")
        except Exception:
            continue
        if a.ndim == 2 and a.shape[0] >= 3:
            runs.append(a)
    return runs

def median_curve(runs):
    ref = max(runs, key=lambda a: a.shape[0])
    tg = ref[:, 1]
    M = np.vstack([np.interp(tg, r[:, 1], r[:, 2]) for r in runs])
    return tg, np.median(M, axis=0)

base, fix = load(base_dir), load(fix_dir)
print(f"loaded baseline={len(base)} fixed={len(fix)} runs")
tb, mb = median_curve(base)
tf, mf = median_curve(fix)


ta = 1e7
ia = int(np.argmin(np.abs(tb - ta)))
t0, e0 = tb[ia], mb[ia]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True, sharey=True)
panels = [
    (axes[0], base, tb, mb, (tf, mf), "baseline (no fix)", "C3"),
    (axes[1], fix,  tf, mf, (tb, mb), "with fix",          "C2"),
]
for ax, runs, tg, med, (to, mo), title, col in panels:
    for r in runs:
        ax.loglog(r[:, 1], r[:, 2], color=col, alpha=0.10, lw=0.6)
    ax.loglog(to, mo, color="0.5", lw=1.6, ls="--", label="other median")
    ax.loglog(tg, med, color=col, lw=2.6, label="median")
    tt = np.array([t0, tg[-1]])
    ax.loglog(tt, e0*np.sqrt(tt/t0), "k:",  lw=1.2, label=r"$\propto\sqrt{t}$")
    ax.loglog(tt, e0*(tt/t0),        "k--", lw=1.0, label=r"$\propto t$")
    ax.set_title(title); ax.set_xlabel("time [yrs]")
    ax.grid(True, which="both", alpha=0.2)
    ax.set_ylim(3e-14, 2e-9)
    ax.legend(loc="upper left", fontsize=8, ncol=2)
axes[0].set_ylabel("relative energy error")
fig.tight_layout()
fig.savefig(outpng, dpi=140)
print("wrote", outpng)
