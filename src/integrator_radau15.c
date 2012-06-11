/**
 * @file 	integrator.c
 * @brief 	RADAU15 integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the radau15 integration scheme.  
 * This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * See Everhart, 1985, ASSL Vol. 115, IAU Colloq. 83, Dynamics of 
 * Comets, Their Origin and Evolution, 185.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Dave Spiegel.
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "problem.h"
#ifdef TREE
#error RADAU15 integrator not working with TREE module.
#endif
#ifdef MPI
#error RADAU15 integrator not working with MPI.
#endif


void integrator_part1(){
	// Do nothing here. This is for the first drift part in a leapfrog-like DKD integrator.
}

// This function updates the acceleration on all particles. 
// It uses the current position and velocity data in the (struct particle*) particles structure.
// Note: this does currently not work with MPI or any TREE module.
void integrator_update_acceleration(){
  gravity_calculate_acceleration();
  if (problem_additional_forces) problem_additional_forces();
}

void integrator_part2(){
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Available variables at the beginning of this function:
  // t: 			current time
  // dt: 			suggested times step (doesn't have to be fixed)
  // particles:		x,y,z,vx,vy,vz,ax,ay,az of all particles at the beginning of the timestep.
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Input variables in the fortran program: x,v,TF,xl,ll,nv,NCLASS,order(NOR)
  double tval,pw;
  double x[1],v[1];
  double f1[18],fj[18];
  double c[21],d[21],r[21],y[18],z[18];
  double b[7,18],g[7][18],E[7][18],bd[7][18];
  double h[8],w[7],u[7];
  double nw[8] = {0,0,1,3,6,10,15,21};
  bool npq,nsf,nper,ncl,nes;
  double zero=0.,half=0.5,one=1.,sr=1.4;
  int i,j,k,l,m,n,la,lb,lc;
  double w,ww,pw,xl,dir;

  // These H values are the Gauss-Radau spacings,
  // scaled to the range 0 to 1, for integrating to order 15.  H[0] = 0. always.
  h = { 0., 0.05626256053692215, 0.18024069173689236, \
	0.35262471711316964, 0.54715362633055538, 0.73421017721541053, \
	0.88532094683909577, 0.97752061356128750 };
  // The sum of these values should be 3.7333333333333333 (checked, it is)
  nper=false;
  nsf=false;
  ncl=(nclass==1);
  npq=(nclass<2);
  // y'  = F[y,t],    NCL=true.  y'' = F[y,t],    NCL=false. y'' = F[y',y,t], NCL=false.
  // NCLASS=1,        NPQ=true.  NCLASS=-2,       NPQ=true.  NCLASS=2,        NPQ=false.
  // NSF is false on starting sequence, otherwise true
  // NPER is true only on last sequence of the integration
  // NES is true only if ll is negative.  Then the sequence size is xl.
  dir=one;
  if (tf < zero)
    dir=-one;
  new = (ll<0);
  xl = abs(xl)*dir;
  pw = 1./9.;
  // Evaluate the constants in the W-,U-,C-,D-, and R-vectors
  for (n=1; n<=7; ++n) {
    ww = n + n*n;
    if (ncl)
      ww=n;
    w[n-1] = one/ww;
    ww = n;
    u[n-1] = one/ww;
  }
  for (k=0; k<=(nv-1); ++k) {
    if (ncl)
      v[k] = zero;
    for (l=0; l<=6; ++l) {
      bd[l,k] = zero;
      b[l,k]  = zero;
    }
  }
  w1 = half;
  if (ncl)
    w1 = one;
  c[0] = -h[1];
  d[0] =  h[1];
  r[0] = one / (h[2] - h[1]);
  la = 0; // should these be zero or one?
  lc = 0;
  for (k = 2; k<=6; ++k) {
    lb = la;
    la = lc + 1;
    lc = nw[k+1];
    c[la] = -h[k]*c[lb];
    c[lc] =  c[la-1]-h[k];
    d[la] =  h[1]*d[lb];
    d[lc] = -c[lc];
    r[la] = one / (h[k+1] - h[1]);
    r[lc] = one / (h[k+1] - h[k]);
    if (k == 3)
      break;
    for (l = 3; l<=(k-1); ++l) {
      ld    = la + l-3;
      le    = lb + l-4;
      c[ld] = c[le] - h[k]*c[le+1];
      d[ld] = d(le) + h[l-1]*d[le+1];
      r[ld] = one / (h[k+1] - h[l-1]);
    }
  }
  ss = pow(10.,(-ll));
  //  The statements above are used only once in an integration to set up the
  //  constants. They use less than a second of execution time.  Next set in
  //  a reasonable estimate to TP based on experience. Same sign as DIR.
  //  An initial first sequence size can be set with XL even with LL positive.
  tp = 0.1d0*dir;
  if (xl.ne.zero) tp=xl;
  if (nes) tp=xl;
  if (tp/tf > half) tp=half*tf;
  ncount=0;
  //     WRITE (*,3)
  //*  3  FORMAT(/' No. of calls, Every 10th seq.X(1),T,TM,TF')
  //  An * is the symbol for writing on the monitor. The printer is unit 4.
  //  Line 4000 is the starting place of the first sequence.
  //
  // How do I handle "go to" statements??
  // 4000
  ns = 0;
  nf = 0;
  ni = 6;
  tm = zero;
  // DSS 8 June 2012: Replace the following
  //  call force (x, v, zero, f1)
  nf = nf+1;
  // Line 722 is begins every sequence after the first. First find new beta-
  //  values from the predicted B-values, following Eq. (2.7) in text.
  // 722
  for (k=0; k<=(nv-1); ++k) {
    g[0,k] = b[0,k] + d[0]*b[1,k] + d[ 1]*b[2,k] + d[3]*b[3,k]  + d[6]*b[4,k] + d[10]*b[5,k] + d[15]*b[6,k];
    g[1,k] =               b(1,k) + d[ 2]*b[2,k] + d[4]*b[3,k]  + d[7]*b[4,k] + d[11]*b[5,k] + d[16]*b[6,k];
    g[2,k] = b[2,k] + d[5]*b[3,k] + d[ 8]*b[4,k] + d[12]*b[5,k] + d[17]*b[6,k];
    g[3,k] =               b[3,k] + d[ 9]*b[4,k] + d[13]*b[5,k] + d[18]*b[6,k];
    g[4,k] =                              b[4,k] + d[14]*b[5,k] + d[19]*b[6,k];
    g[5,k] =                                             b[5,k] + d[20]*b[6,k];
    g[6,k]=                                                             b[6,k];
    }
  t = tp;
  t2 = t*t;
  if (ncl) t2=t;
  tval = abs(t);
  // DSS 8 June 2012: Done up to here 
  //     IF(NS/10*10.EQ.NS) WRITE(*,7) NF,NS,X(1),X(2),T,TM,TF
  //  7  FORMAT(1X,2I6,3F12.5,1P,2E10.2)
  //  Loop 175 is 6 iterations on first sequence and two iterations therafter.
      do 175 m=1,ni
	   //  Loop 174 is for each substep within a sequence.
      do 174 j=2,8
      jd=j-1
      jdm=j-2
      s=h(j)
      q=s
      if(ncl) q=one
*  Use Eqs. (2.9) and (2.10) of text to predict positions at each aubstep.
*  These collapsed series are broken into two parts because an otherwise
*  excellent  compiler could not handle the complicated expression.
      do 130 k=1,nv
      a=w(3)*b(3,k)+s*(w(4)*b(4,k)+s*(w(5)*b(5,k)+s*(w(6)*b(6,k)+
     v   s*w(7)*b(7,k))))
      y(k)=x(k)+q*(t*v(k)+t2*s*(f1(k)*w1+s*(w(1)*b(1,k)+s*(w(2)*b(2,k)
     x  +s*a))))
      if(npq) go to 130
*  Next are calculated the velocity predictors need for general class II.
      a=u(3)*b(3,k)+s*(u(4)*b(4,k)+s*(u(5)*b(5,k)+s*(u(6)*b(6,k)+
     t    s*u(7)*b(7,k))))
      z(k)=v(k)+s*t*(f1(k)+s*(u(1)*b(1,k)+s*(u(2)*b(2,k)+s*a)))
 130  continue
*  Find forces at each substep.
      call force(y,z,tm+s*t,fj)
      nf=nf+1
      do 171 k=1,nv
*  Find G-value for the force FJ found at the current substep. This
*  section, including the many-branched GOTO, uses Eq. (2.4) of text.
      temp=g(jd,k)
      gk=(fj(k)-f1(k))/s
      go to (102,102,103,104,105,106,107,108),j
 102  g(1,k)=gk
      go to 160
 103  g(2,k)=(gk-g(1,k))*r(1)
      go to 160
 104  g(3,k)=((gk-g(1,k))*r(2)-g(2,k))*r(3)
      go to 160
 105  g(4,k)=(((gk-g(1,k))*r(4)-g(2,k))*r(5)-g(3,k))*r(6)
      go to 160
 106  g(5,k)=((((gk-g(1,k))*r(7)-g(2,k))*r(8)-g(3,k))*r(9)-
     x     g(4,k))*r(10)
      go to 160
 107  g(6,k)=(((((gk-g(1,k))*r(11)-g(2,k))*r(12)-g(3,k))*r(13)-
     x     g(4,k))*r(14)-g(5,k))*r(15)
      go to 160
 108  g(7,k)=((((((gk-g(1,k))*r(16)-g(2,k))*r(17)-g(3,k))*r(18)-
     x     g(4,k))*r(19)-g(5,k))*r(20)-g(6,k))*r(21)
*  Upgrade all B-values.
 160  temp=g(jd,k)-temp
      b(jd,k)=b(jd,k)+temp
*  TEMP is now the improvement on G(JD,K) over its former value.
*  Now we upgrade the B-value using this dfference in the one term.
*  This section is based on Eq. (2.5).
      go to (171,171,203,204,205,206,207,208),j
 203  b(1,k)=b(1,k)+c(1)*temp
      go to 171
 204  b(1,k)=b(1,k)+c(2)*temp
      b(2,k)=b(2,k)+c(3)*temp
      go to 171
 205  b(1,k)=b(1,k)+c(4)*temp
      b(2,k)=b(2,k)+c(5)*temp
      b(3,k)=b(3,k)+c(6)*temp
      go to 171
 206  b(1,k)=b(1,k)+c(7)*temp
      b(2,k)=b(2,k)+c(8)*temp
      b(3,k)=b(3,k)+c(9)*temp
      b(4,k)=b(4,k)+c(10)*temp
      go to 171
 207  b(1,k)=b(1,k)+c(11)*temp
      b(2,k)=b(2,k)+c(12)*temp
      b(3,k)=b(3,k)+c(13)*temp
      b(4,k)=b(4,k)+c(14)*temp
      b(5,k)=b(5,k)+c(15)*temp
      go to 171
 208  b(1,k)=b(1,k)+c(16)*temp
      b(2,k)=b(2,k)+c(17)*temp
      b(3,k)=b(3,k)+c(18)*temp
      b(4,k)=b(4,k)+c(19)*temp
      b(5,k)=b(5,k)+c(20)*temp
      b(6,k)=b(6,k)+c(21)*temp
 171  continue
 174  continue
      if(nes.or.m.lt.ni) go to 175
*  Integration of sequence is over. Next is sequence size control.
      hv=zero
      do 635 k=1,nv
 635  hv=dmax1(hv,dabs(b(7,k)))
      hv=hv*w(7)/tval**7
 175  continue
      if (nsf) go to 180
      if(.not.nes) tp=(ss/hv)**pw*dir
      if(nes) tp=xl
      if(nes) go to 170
      if(tp/t.gt.one) go to 170
   8  format (2x,2i2,2d18.10)
      tp=.8d0*tp
      ncount=ncount+1
      if(ncount.gt.10) return
*     IF(NCOUNT.GT.1) WRITE (4,8) NOR,NCOUNT,T,TP
*  Restart with 0.8x sequence size if new size called for is smaller than
*  originally chosen starting sequence size on first sequence.
      go to 4000
 170  nsf=.true.
* Loop 35 finds new X and V values at end of sequence using Eqs. (2.11),(2.12)
 180  do 35 k=1,nv
      x(k)=x(k)+v(k)*t+t2*(f1(k)*w1+b(1,k)*w(1)+b(2,k)*w(2)+b(3,k)*w(3)
     x    +b(4,k)*w(4)+b(5,k)*w(5)+b(6,k)*w(6)+b(7,k)*w(7))
      if(ncl) go to 35
      v(k)=v(k)+t*(f1(k)+b(1,k)*u(1)+b(2,k)*u(2)+b(3,k)*u(3)
     v    +b(4,k)*u(4)+b(5,k)*u(5)+b(6,k)*u(6)+b(7,k)*u(7))
  35  continue
      tm=tm+t
      ns=ns+1
*  Return if done.
      if(.not.nper) go to 78
*     WRITE(*,7) NF,NS,X(1),X(2),T,TM,TF
*     WRITE(4,7) NF,NS
      return
*  Control on size of next sequence and adjust last sequence to exactly
*  cover the integration span. NPER=.TRUE. set on last sequence.
78    call force (x,v,tm,f1)
      nf=nf+1
      if(nes) go to 341
      tp=dir*(ss/hv)**pw
      if(tp/t.gt.sr) tp=t*sr
 341  if(nes) tp=xl
      if(dir*(tm+tp).lt.dir*tf-1.d-8) go to 77
      tp=tf-tm
      nper=.true.
*  Now predict B-values for next step. The predicted values from the preceding
*  sequence were saved in the E-matrix. Te correction BD between the actual
*  B-values found and these predicted values is applied in advance to the
*  next sequence. The gain in accuracy is significant. Using Eqs. (2.13):
  77  q=tp/t
      do 39 k=1,nv
      if(ns.eq.1) go to 31
      do 20 j=1,7
  20  bd(j,k)=b(j,k)-e(j,k)
  31  e(1,k)=      q*(b(1,k)+ 2.d0*b(2,k)+ 3.d0*b(3,k)+
     x           4.d0*b(4,k)+ 5.d0*b(5,k)+ 6.d0*b(6,k)+ 7.d0*b(7,k))
      e(2,k)=                q**2*(b(2,k)+ 3.d0*b(3,k)+
     y           6.d0*b(4,k)+10.d0*b(5,k)+15.d0*b(6,k)+21.d0*b(7,k))
      e(3,k)=                             q**3*(b(3,k)+
     z           4.d0*b(4,k)+10.d0*b(5,k)+20.d0*b(6,k)+35.d0*b(7,k))
      e(4,k)=   q**4*(b(4,k)+ 5.d0*b(5,k)+15.d0*b(6,k)+35.d0*b(7,k))
      e(5,k)=                q**5*(b(5,k)+ 6.d0*b(6,k)+21.d0*b(7,k))
      e(6,k)=                             q**6*(b(6,k)+ 7.d0*b(7,k))
      e(7,k)=                                           q**7*b(7,k)
      do 39 l=1,7
  39  b(l,k)=e(l,k)+bd(l,k)
*  Two iterations for every sequence after the first.
      ni=2
      go to 722
      end
			      // End of the Everhart bit.  
    
  // Do the real work here.
	
  // Note that when calling integrator_update_acceleration() within this function, the 
  // correct position and velocities must be stored in the particles[] array. 
  // The updated accelerations will then also be stored in the particles[] array.
  // It's like if the function F(y,y',t) can only take particles[] as an argument, nothing else.
  // You may, however, change the pointer to the particles array temporarily to help you with that. 
  // For example do something like this:
  struct particle* old_particles = particles;		// Save old particle array pointer.
  particles = malloc(sizeof(struct particle)*N);	// Create space for temporary array.
  for (int i=0;i<N;i++){
    particles[i] = old_particles[i];		// Copy old particle data.
    // Let's do an implicit backward Euler method as an example.
    // Note: this is a bad integrator which need a really small timestep.
    // It's not time-reversible and particles fall into the star quickly.
    particles[i].x  += dt*particles[i].vx;		// First guess.	
    particles[i].y  += dt*particles[i].vy;		
    particles[i].z  += dt*particles[i].vz;		
    particles[i].vx += dt*particles[i].ax;		
    particles[i].vy += dt*particles[i].ay;		
    particles[i].vz += dt*particles[i].az;		
  }

  // Copying comments from Everhart preamble
  //
  // Integrator by E Everhart, Physics Department, University of Denver
  // This 15th-order version is called RA15.  Order NOR is 15.
  // y'  = F[y,t]    is NCLASS = 1,    y'' = F[y,t] is NCLASS = -2
  // y'' = F[y',y,t] is NCLASS = 2
  // TF is t[final] - t[initial]. (Negative when integrating backward.)
  // NV = the number of simultaneous differential equations.
  // Change dimensioning if NV is greater than 18. (DSS comment: why??)
  // LL controls accuracy.  Thus SS=10.^(-LL) controls the size of the last term in a series.
  // Try LL=8 and work up or down from there.
  // However, if LL<0, then XL is the constant sequence size used.
  // A non-zero XL sets the size of the first sequence regardless of LL's sign.
  // X and V enter as the starting position-velocity vector (values of y and y' at t=0)
  // and they output as the final position-velocity vector.
  // Integration is in double precision.  A 64-bit double word is assumed.

  for (int iterations=0;iterations<1000;iterations++){
    integrator_update_acceleration(); // This updates ax,ay,az in the particles[] array.
    for (int i=0;i<N;i++){
      particles[i].x  = old_particles[i].x  + dt*particles[i].vx;	// Each iter	
      particles[i].y  = old_particles[i].y  + dt*particles[i].vy;		
      particles[i].z  = old_particles[i].z  + dt*particles[i].vz;		
      particles[i].vx = old_particles[i].vx + dt*particles[i].ax;		
      particles[i].vy = old_particles[i].vy + dt*particles[i].ay;		
      particles[i].vz = old_particles[i].vz + dt*particles[i].az;		
    }
  }
  for (int i=0;i<N;i++){
    old_particles[i] = particles[i];		// Copy particles back to the old particle data structure.
  }
  free(particles);					// Release memory of temporary array.
  particles = old_particles;				// Change pointer back to old particle array pointer.


  // Advance the time. Because the function is using an adaptive timestep, use the dt that has been computed self-consistely.
  t+=dt; 	



  /////////////////////////////////////////////////////////////////////////////////////////////////
  // These vailables should be set at the end of the timestep:
  // dt: 			actual timestep done
  // particles:		x,y,z,vx,vy,vz of all particles at the end of the timestep.
  /////////////////////////////////////////////////////////////////////////////////////////////////
}
	

