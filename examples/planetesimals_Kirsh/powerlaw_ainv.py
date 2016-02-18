#Quick macro to double check that I'm drawing from the correct CDF to get a a^-1 powerlaw.

import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import random

Nsamples = 100000
data = np.zeros(Nsamples)
min = 14.5
max = 35.5
random.seed(11)

for i in xrange(0,Nsamples):
    data[i] = np.exp(random.random()*np.log(max/min) + np.log(max/min))

n, bins, patches = plt.hist(data, 50, facecolor='green', alpha=0.75)
plt.plot(bins, (n[0]/bins[0]**(-1))*bins**(-1), color='black')
plt.show()