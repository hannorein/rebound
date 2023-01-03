import rebound
import unittest
import math
import warnings

class TestRotations(unittest.TestCase):
    def test_init(self):
        r = rebound.Rotation(angle=math.pi/2, axis=[0,0,1.])
        res = r*[1,0,0] 
        self.assertAlmostEqual(res[0], 0, delta=1e-15)
        self.assertAlmostEqual(res[1], 1, delta=1e-15)
        self.assertAlmostEqual(res[2], 0, delta=1e-15)
    
    def test_inverse(self):
        r = rebound.Rotation(ix=0.1, iy=0.2, iz=0.3, r=0.5)
        res = r*r.inverse()*[1,1,1]
        self.assertAlmostEqual(res[0], 1, delta=1e-15)
        self.assertAlmostEqual(res[1], 1, delta=1e-15)
        self.assertAlmostEqual(res[2], 1, delta=1e-15)
    
    def test_from_to(self):
        fromv = [0.5,0.5,-0.5]
        tov = [0.5,0.5,0.5]
        r = rebound.Rotation.from_to(fromv=fromv, tov=tov)
        res = r*fromv
        self.assertAlmostEqual(res[0], tov[0], delta=1e-15)
        self.assertAlmostEqual(res[1], tov[1], delta=1e-15)
        self.assertAlmostEqual(res[2], tov[2], delta=1e-15)

    def test_from_to_edge(self):
        q = rebound.Rotation.from_to([1,0,0], [-1,0,0])
        res = q*[1,0,0] 
        self.assertAlmostEqual(res[0], -1, delta=1e-15)
        self.assertAlmostEqual(res[1], 0, delta=1e-15)
        self.assertAlmostEqual(res[2], 0, delta=1e-15)

    def test_to_orbit(self):
        sim = rebound.Simulation()
        a, e, inc, Omega, omega = 1, 0.1, 0.2, 0.3, 0.4
        sim.add(m=1)
        sim.add(a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=0)
        sim.add(a=a,e=e)
        # when we rotate our particle's xyz orbit (x=toward peri, z=orb normal) we should get (a(1-e), 0, 0)
        r = rebound.Rotation.orbit(Omega=Omega, inc=inc, omega=omega)
        res = r*sim.particles[2].xyz
        self.assertAlmostEqual(res[0], sim.particles[1].x, delta=1e-15)
        self.assertAlmostEqual(res[1], sim.particles[1].y, delta=1e-15)
        self.assertAlmostEqual(res[2], sim.particles[1].z, delta=1e-15)
    
    def test_to_new_axes(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(a=1,inc=math.pi/2,Omega=math.pi/2)
        orbnorm = [1,0,0] # along x axis for above orbit
        ascnode = [0,1,0] # along y axis for above orbit
        # since omega=f=0, particle is at node. If we rotate into a frame with newx along the node, and z along orbnorm, should get (1,0,0)
        r = rebound.Rotation.to_new_axes(newz=orbnorm, newx=ascnode)
        res = r*sim.particles[1].xyz
        self.assertAlmostEqual(res[0], 1, delta=1e-15)
        self.assertAlmostEqual(res[1], 0, delta=1e-15)
        self.assertAlmostEqual(res[2], 0, delta=1e-15)
    
    def test_to_new_axes_not_perp(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(a=1,inc=math.pi/2,Omega=math.pi/2)
        orbnorm = [1,0,0] # along x axis for above orbit
        # since omega=f=0, particle is at node. If we rotate into a frame with newx along the node, and z along orbnorm, should get (1,0,0)
        delta = 0.01 # deviation from orthogonal newx ([0,1,0])
        r = rebound.Rotation.to_new_axes(newz=orbnorm, newx=[0,1,delta])
        res = r*sim.particles[1].xyz
        self.assertAlmostEqual(res[0], 1, delta=delta)
        self.assertAlmostEqual(res[1], 0, delta=delta)
        self.assertAlmostEqual(res[2], 0, delta=delta)
    
    
    def test_to_new_axes_no_x_unnormalized(self):
        newz = [-1.,-1.,-1.]
        mag = (newz[0]**2 + newz[1]**2 + newz[2]**2)**0.5
        r = rebound.Rotation.to_new_axes(newz=newz)
        res = r*newz
        self.assertAlmostEqual(res[0], 0, delta=1e-15)
        self.assertAlmostEqual(res[1], 0, delta=1e-15)
        self.assertAlmostEqual(res[2], mag, delta=1e-15)
    
    def test_mul(self):
        r = rebound.Rotation(angle=math.pi/2, axis=[0,0,1.])
        rhalf = rebound.Rotation(angle=math.pi/4, axis=[0,0,1.])
        res1 = r*[1,0,0] 
        res2 = (rhalf*rhalf)*[1,0,0]
        self.assertAlmostEqual(res1[0]-res2[0], 0, delta=1e-15)
        self.assertAlmostEqual(res1[1]-res2[1], 0, delta=1e-15)
        self.assertAlmostEqual(res1[2]-res2[2], 0, delta=1e-15)
    
    def test_to_from_spherical(self):
        mag, theta, phi = 3, math.pi/3, -math.pi/4
        vec = rebound.spherical_to_xyz(mag, theta, phi)
        mag2, theta2, phi2 = rebound.xyz_to_spherical(vec)
        self.assertAlmostEqual(mag, mag2, delta=1e-15)
        self.assertAlmostEqual(theta, theta2, delta=1e-15)
        self.assertAlmostEqual(phi, phi2, delta=1e-15)
        
if __name__ == "__main__":
    unittest.main()
