import rebound
import numpy as np

def compute_satellite_orbits(self, outname='satellite_orbits.dat'):
   '''Compute the semimajor axis, eccentricity, and inclination of 
   each satellite in the simulation. The resulting orbital elements 
   are saved to a unique file corresponding to the timestep at which 
   this function is called.

   outname = string representing the filename that the satellite orbits 
   are saved to.
   '''
   # Define Jupiter as a special particle
   Jup = self.particles[1]
   
   # Compute each satellite's semimajor axis, eccentricity, and inclination
   Ntest = len(self.particles)-2
   orbitalparams = np.zeros((Ntest,4))
   for i in range(2,Ntest+2):
      
      # Move to planetocentric coordinates
      self.particles[i].x  -= Jup.x
      self.particles[i].y  -= Jup.y
      self.particles[i].z  -= Jup.z
      self.particles[i].vx -= Jup.vx
      self.particles[i].vy -= Jup.vy
      self.particles[i].vz -= Jup.vz

      # Compute planetocentric orbital parameters 
      mu = self.G * Jup.m
      hx = (self.particles[i].y*self.particles[i].vz - 
            self.particles[i].z*self.particles[i].vy)
      hy = (self.particles[i].z*self.particles[i].vx - 
            self.particles[i].x*self.particles[i].vz)
      hz = (self.particles[i].x*self.particles[i].vy - 
            self.particles[i].y*self.particles[i].vx)
      h  = np.sqrt(hx*hx + hy*hy + hz*hz)
      v  = np.sqrt(self.particles[i].vx**2 +
                   self.particles[i].vy**2 +
                   self.particles[i].vz**2 )
      r  = np.sqrt(self.particles[i].x**2 +
                   self.particles[i].y**2 +
                   self.particles[i].z**2 )
      vr = (self.particles[i].x * self.particles[i].vx +
            self.particles[i].y * self.particles[i].vy + 
            self.particles[i].z * self.particles[i].vz )/r
      ex = 1./mu * ((v*v-mu/r)*self.particles[i].x - r*vr*self.particles[i].vx)
      ey = 1./mu * ((v*v-mu/r)*self.particles[i].y - r*vr*self.particles[i].vy)
      ez = 1./mu * ((v*v-mu/r)*self.particles[i].z - r*vr*self.particles[i].vz)
      
      # Compute the desired orbital elements 
      a   = -mu/(v*v-2.*mu/r)
      ecc = np.sqrt(ex*ex + ey*ey + ex*ez)
      inc = np.arccos(hz/h)
      if (a > 0 and ecc <= 1):
         orbitalparams[i-2,1] = a
         orbitalparams[i-2,2] = ecc
         orbitalparams[i-2,3] = inc
      else:
         orbitalparams[i-2,1] = np.nan
         orbitalparams[i-2,2] = np.nan
         orbitalparams[i-2,3] = np.nan     
         
      # Revert to heliocentric coordinates to continue the integration
      self.particles[i].x  += Jup.x
      self.particles[i].y  += Jup.y
      self.particles[i].z  += Jup.z
      self.particles[i].vx += Jup.vx
      self.particles[i].vy += Jup.vy
      self.particles[i].vz += Jup.vz

   # Save satellite orbits
   orbitalparams[:,0] = np.repeat(self.t, Ntest)
   np.savetxt(outname, orbitalparams, delimiter='\t', fmt='%.6e',
	      header='Only Satellite orbits:\ntime (days)\nsemimajor axis (AU)\neccentricity\ninclination (radians)')
