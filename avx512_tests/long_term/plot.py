import numpy as np
import matplotlib.pyplot as plt

logs = []
minl = 10000000
for i in range(240):
    l = np.loadtxt("6day/log_%05d.txt"%i,delimiter="\t")
    logs.append(l)
    if len(l) < minl:
        minl = len(l)

logs2 = []
for i in range(240):
        logs2.append(logs[i][:minl])

logs = np.array(logs2)

logs1 = []
minl1 = 10000000
for i in range(240):
    l = np.loadtxt("1day/log_1day_%05d.txt"%i,delimiter="\t")
    logs1.append(l)
    if len(l) < minl1:
        minl1 = len(l)

logs2 = []
for i in range(240):
        logs2.append(logs1[i][:minl1])

logs1 = np.array(logs2)

print("unstable", np.sum(np.max(logs[:,:,3], axis=1)>0.6), "/ 240")
print("unstable", np.sum(np.max(logs1[:,:,3], axis=1)>0.7), "/ 240")

fig, ax = plt.subplots(2,1,figsize=(5,7),sharex=True)
ax[0].set_ylabel("relative energy error")
ax[1].set_ylabel("max e")
ax[1].set_xlabel("time [Gyr]")
ax[0].set_xscale("log")
ax[0].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_ylim([0,1])
ax[0].set_xlim([1e-6,5])
time = logs[0,:,0]/2/np.pi/1e9
percentiles = np.nanpercentile(logs[:,:,1], [10,50,90], axis=0)
median = np.nanmedian(logs[:,:,1], axis=0)
time1 = logs1[0,:,0]/2/np.pi/1e9
percentiles1 = np.nanpercentile(logs1[:,:,1], [10,50,90], axis=0)
median1 = np.nanmedian(logs1[:,:,1], axis=0)

ax[0].plot(time1, median1, label="1 day", color="blue");
ax[0].fill_between(time1, percentiles1[0],percentiles1[2], color="blue", alpha = 0.2);
ax[0].plot(time, median, label="6 day", color="green");
ax[0].fill_between(time, percentiles[0],percentiles[2], color="green", alpha = 0.2);

ax[1].plot(time1, logs1[:,:,3].T, color="blue", alpha=0.2);
ax[1].plot(time, logs[:,:,3].T, color="green", alpha=0.2);

ax[0].plot(time, 1e-16*np.sqrt(time*(365e9/6)), color="red",ls=":", label="6day error");
ax[0].plot(time, 1e-16*np.sqrt(time*365e9), color="red", label="1day error");
fig.tight_layout()
fig.savefig("out.pdf")
