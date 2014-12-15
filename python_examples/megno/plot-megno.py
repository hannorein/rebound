#!/usr/bin/env python
import math
import os
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as pl
import matplotlib.colors as mcolors

heatmap = np.load("megno.npy")
shape    = heatmap.shape
print shape

y = heatmap[0,:,0]
x = heatmap[:,0,1]
extent = [x.min(), x.max(), y.min(), y.max()]
np.seterr(divide='ignore')
zMegno = (heatmap[:,:,2])

pl.imshow(zMegno,vmin=0., vmax=4., origin='lower', aspect='auto', interpolation='nearest', cmap="Blues", extent=extent)
cb1 = pl.colorbar()
cb1.solids.set_rasterized(True)

cb1.set_label("MEGNO stability indicator $\\langle Y \\rangle$")
pl.xlim(x.min(),x.max())
pl.ylim(y.min(),y.max())
pl.xlabel("$P$")
pl.ylabel("$e$")

pl.savefig("megno.pdf")

os.system("open megno.pdf")


