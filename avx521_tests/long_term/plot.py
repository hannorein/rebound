import numpy as np
import matplotlib.pyplot as plt

logs = []
minl = 10000000
for i in range(240):
    l = np.loadtxt("/scratch/rein/whfast512_tests/log_%05d.txt"%i,delimiter="\t")
    logs.append(l)
    if len(l) < minl:
        minl = len(l)

logs2 = []
for i in range(240):
        logs2.append(logs[i][:minl])

logs = np.array(logs2)

fig, ax = plt.subplots(2,1,figsize=(5,7),sharex=True)
ax[0].set_ylabel("relative energy error")
ax[1].set_ylabel("max e")
ax[1].set_xlabel("time [Gyr]")
ax[0].set_xscale("log")
ax[0].set_yscale("log")
ax[1].set_xscale("log")
time = logs[0,:,0]/2/np.pi/1e9
percentiles = np.percentile(logs[:,:,1], [10,50,90], axis=0)
median = np.median(logs[:,:,1], axis=0)

ax[0].plot(time, median, label="6 day", color="green");
ax[0].fill_between(time, percentiles[0],percentiles[2], color="green", alpha = 0.2);

ax[1].plot(time, logs[:,:,3].T, color="green", alpha=0.2);
fig.tight_layout()
fig.savefig("out.pdf")
