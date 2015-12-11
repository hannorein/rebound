__author__ = 'Alysa Obertas'
__email__ = 'obertas@astro.utoronto.ca'

# plot.py
# plots orbit of planet and test particle in intertial and
# co-rotating frame, along with contours of the Jacobi constant
# and the zero-velocity curves
# Written by Alysa Obertas December 2015

import numpy as np
import matplotlib.pyplot as plt

###########################################################
## load file

infile = "three_horseshoe.npz"
data = np.load(infile)

mu = 0.1

times = data['t']

a = data['a']
e = data['e']
x = data['x']
y = data['y']
vxtest = data['vxtest']
vytest = data['vytest']

r = np.sqrt(np.square(x) + np.square(y))

###########################################################
## convert interial coordinates to co-rotating coordinates

theta = np.arctan2(y[0],x[0]) # phase angle of planet

xr = x * np.cos(theta) + y * np.sin(theta) # rotates x by -theta
yr = -x * np.sin(theta) + y * np.cos(theta)

###########################################################
## trajectory of planet and test particle in interial frame

plt.figure()
plt.plot(-mu,0,'*')
plt.plot(x[0],y[0])
plt.plot(x[1],y[1])
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.axis('equal')
plt.show()

###########################################################
## contours of the co-rotating potential (zero velocity curves)

def jacobi(x,y,v=0):
	# calculates the jacobi constant given test particle
	# positions and speed in the co-rotating frame
	r1 = np.sqrt(np.square(-mu-x)   + np.square(y)) # distance between test particle and star
	r2 = np.sqrt(np.square(1.-mu-x) + np.square(y)) # distance between test particle and planet

	return -(np.square(x) + np.square(y))/2. - (1.-mu)/r1 - mu/r2 + np.square(v)/2.

xx = np.arange(-4,4,0.005)
yy = np.arange(-4,4,0.005)
XX,YY = np.meshgrid(xx,yy)

JACOBI_GRID = jacobi(XX,YY)

plt.figure()
plt.plot(-mu,0,'*')
plt.plot(xr[0][0],yr[0][0],'o')
TT = plt.contour(XX,YY,JACOBI_GRID,[-7.3,-3.5,-2,-1.7,-1.5,-1.6])
plt.clabel(TT, fontsize=9, inline=1)
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.axis('equal')
plt.show()

###########################################################
## trajectory of test particle in co-rotating frame

def v_iner_to_rot(x,y,vx,vy):
	# calculates the velocity in the co-rotating frame given
	# the positions and velocities in the inertial frame
	# assumes mean motion of 1
	vxrot = vx + y
	vyrot = vy - x

	return np.sqrt(np.square(vxrot) + np.square(vyrot))

vrot = v_iner_to_rot(x[1][0],y[1][0],vxtest,vytest)

C = jacobi(x[1][0],y[1][0],vrot)

plt.figure()
plt.plot(-mu,0,'*')
plt.plot(xr[0][0],yr[0][0],'o')
plt.plot(xr[1],yr[1])
CS = plt.contour(XX,YY,JACOBI_GRID,[C])
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.axis('equal')
plt.show()
