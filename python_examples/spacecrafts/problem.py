#import rebound module
import rebound
import numpy as np
import matplotlib.pyplot as plt


#include objects from solar system in simulation
solar_system_objects = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "-125544", "-98","-24",  "-31", "-32", "-23"]
objects_name = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "ISS", "New Horizons", "Pioneer 11", "Voyager 1", "Voyager 2", "Pioneer 10"] #only used for labelling. Therefore, it is important that they are in the same order as the objects above. 
sim = rebound.Simulation()
sim.add(solar_system_objects)

#integration
sim.integrator = "whfast"
sim.set_dt = 0.01

#creation x and y arrays.
Nout = 1000
times = np.linspace(0,70.*np.pi,Nout) # integration for 35 years 
x = np.zeros((sim.N,Nout))
y = np.zeros((sim.N,Nout))
ps = sim.particles

#create test array for correct output of the time when spacecraft is leaving solar system
size = np.size(solar_system_objects) 
test = np.zeros(size)

    
for ti,t in enumerate(times):
    sim.integrate(t)
    for i, p in enumerate(ps):
        x[i][ti] = p.x
        y[i][ti] = p.y
        r=(p.x**2+p.y**2+p.z**2)**0.5
#check if spacecraft leaves solarsystem in this timestep
        if r>123 and test[i]<1:
            print '%s leaves the Solar System in %.2f years.' % (objects_name[i],t/(2*np.pi))
            test[i]=10

#plot spacecraft trajectories
fig = plt.figure(figsize=(6,5))
def plot(zoom):
    ax.set_xlim([-zoom,zoom])
    ax.set_ylim([-zoom,zoom])
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    for i in xrange(0,sim.N):
        plt.plot(x[i],y[i])
        if x[i][-1]*x[i][-1]+y[i][-1]*y[i][-1]>0.01*zoom*zoom or i==0:
            ax.annotate(objects_name[i], xy=(x[i][-1], y[i][-1]),horizontalalignment="center")

            
circle=plt.Circle((0,0),123, color='r',fill=False)
#possible second, zoomed plot for better/different resolution if needed
#ax = plt.subplot(121)
#plot(zoom=1.6)
ax = plt.subplot(111)
plot(zoom=220.)
plt.gcf().gca().add_artist(circle)
plt.savefig("trajectories.pdf")

