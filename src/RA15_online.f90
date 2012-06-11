!  Integrator RADAU by E. Everhart, Physics Department, University of De
!  This 15th-order version, called RA15, is written out for faster execu
!  y'=F(y,t) is  NCLASS=1,  y"=F(y,t) is NCLASS= -2,  y"=F(y',y,t) is NC
!  TF is t(final) - t(initial). It may be negative for backward integrat
!  NV = the number of simultaneous differential equations.              
!  The dimensioning below assumes NV will not be larger than NVX.       
!  LL controls sequence size. Thus SS=10**(-LL) controls the size of a t
!  A typical LL-value is in the range 6 to 12 for this order 11 program.
!  However, if LL.LT.0 then XL is the constant sequence size used.      
!  X and V enter as the starting position-velocity vector, and are outpu
!  the final position-velocity vector.                                  
!  Integration is in double precision. A 64-bit double-word is assumed. 
!                                                                       
      subroutine ra15b(xx,vv,tstart,dt,xl,ll,nv,nclass) 
      implicit double precision (a-h,o-z) 
                                                                        
      include 'dim.h90' 
                                                                        
      parameter (nvx=3*nbtx) 
                                                                        
      dimension xx(nv),vv(nv),f1(nvx),fj(nvx),c(21),d(21),r(21),        &
     &     y(nvx),z(nvx),x(nvx),v(nvx),                                 &
     &     b(7,nvx),g(7,nvx),e(7,nvx),bd(7,nvx),h(8),w(7),u(7),nw(8)    
      logical npq,nsf,nper,ncl,nes 
      data nw/0,0,1,3,6,10,15,21/ 
      data zero, half, one,sr/0.0d0, 0.5d0, 1.0d0,1.4d0/ 
!  These H values are the Gauss-Radau spacings, scaled to the range 0 to
!  for integrating to order 15.                                         
      data h/         0.d0, .05626256053692215d0, .18024069173689236d0, &
     &.35262471711316964d0, .54715362633055538d0, .73421017721541053d0, &
     &.88532094683909577d0, .97752061356128750d0/                       
      if(nv.gt.nvx)stop' **** RA15: NV > NVX ****' 
      if(abs(dt).lt.1.d-8) then 
          do 79 i=1,nv 
          xx(i)=xx(i)+vv(i)*dt 
   79     continue 
          return 
      end if 
      DO k=1,nv 
        x(k)=xx(k) 
        v(k)=vv(k) 
      ENDDO 
      tf=dt 
!  The sum of the H-values should be 3.73333333333333333                
      nper=.false. 
      nsf=.false. 
      ncl=nclass.eq.1 
      npq=nclass.lt.2 
!  y'=F(y,t)  NCL=.TRUE.    y"=F(y,t)  NCL=.FALSE.   y"=F(y',y,t) NCL=.F
!  NCLASS=1   NPQ=.TRUE.    NCLASS= -2 NPQ=.TRUE.    NCLASS= 2    NPQ=.F
!  NSF is .FALSE. on starting sequence, otherwise .TRUE.                
!  NPER is .TRUE. only on last sequence of the integration.             
!  NES is .TRUE. only if LL is negative. Then the sequence size is XL.  
      dir=one 
      if(tf.lt.zero) dir=-one 
      nes=ll.lt.0 
      xl=dir*dabs(xl) 
      pw=1.d0/9.d0 
!  Evaluate the constants in the W-, U-, C-, D-, and R-vectors          
      do 14 n=2,8 
      ww=n+n*n 
      if(ncl) ww=n 
      w(n-1)=one/ww 
      ww=n 
   14 u(n-1)=one/ww 
      do 22 k=1,nv 
      if(ncl) v(k)=zero 
      do 22 l=1,7 
      bd(l,k)=zero 
   22 b(l,k)=zero 
      w1=half 
      if(ncl) w1=one 
      c(1)=-h(2) 
      d(1)=h(2) 
      r(1)=one/(h(3)-h(2)) 
      la=1 
      lc=1 
      do 73 k=3,7 
      lb=la 
      la=lc+1 
      lc=nw(k+1) 
      c(la)=-h(k)*c(lb) 
      c(lc)=c(la-1)-h(k) 
      d(la)=h(2)*d(lb) 
      d(lc)=-c(lc) 
      r(la)=one/(h(k+1)-h(2)) 
      r(lc)=one/(h(k+1)-h(k)) 
      if(k.eq.3) go to 73 
      do 72 l=4,k 
      ld=la+l-3 
      le=lb+l-4 
      c(ld)=c(le)-h(k)*c(le+1) 
      d(ld)=d(le)+h(l-1)*d(le+1) 
   72 r(ld)=one/(h(k+1)-h(l-1)) 
   73 continue 
      ss=10.**(-ll) 
!  The statements above are used only once in an integration to set up t
!  constants. They use less than a second of execution time.  Next set i
!  a reasonable estimate to TP based on experience. Same sign as DIR.   
!  An initial first sequence size can be set with XL even with LL positi
      tp=0.1d0*dir 
      if(xl.ne.zero) tp=xl 
      if(nes) tp=xl 
      if(tp/tf.gt.half) tp=half*tf 
      ncount=0 
!     WRITE (*,3)                                                       
!  3  FORMAT(/' No. of calls, Every 10th seq.X(1),T,TM,TF')             
!  An * is the symbol for writing on the monitor. The printer is unit 4.
!  Line 4000 is the starting place of the first sequence.               
 4000 ns=0 
      nf=0 
      ni=6 
      tm=zero 
      call forceb (x, v, tstart+zero, f1) 
      nf=nf+1 
! Line 722 is begins every sequence after the first. First find new beta
!  values from the predicted B-values, following Eq. (2.7) in text.     
  722 do 58 k=1,nv 
      g(1,k)=b(1,k)+d(1)*b(2,k)+                                        &
     &  d(2)*b(3,k)+d(4)*b(4,k)+d( 7)*b(5,k)+d(11)*b(6,k)+d(16)*b(7,k)  
      g(2,k)=            b(2,k)+                                        &
     &  d(3)*b(3,k)+d(5)*b(4,k)+d( 8)*b(5,k)+d(12)*b(6,k)+d(17)*b(7,k)  
      g(3,k)=b(3,k)+d(6)*b(4,k)+d( 9)*b(5,k)+d(13)*b(6,k)+d(18)*b(7,k) 
      g(4,k)=            b(4,k)+d(10)*b(5,k)+d(14)*b(6,k)+d(19)*b(7,k) 
      g(5,k)=                         b(5,k)+d(15)*b(6,k)+d(20)*b(7,k) 
      g(6,k)=                                      b(6,k)+d(21)*b(7,k) 
   58 g(7,k)=                                                   b(7,k) 
      t=tp 
      t2=t*t 
      if(ncl) t2=t 
      tval=dabs(t) 
!     IF(NS/10*10.EQ.NS) WRITE(*,7) NF,NS,X(1),X(2),T,TM,TF             
!  7  FORMAT(1X,2I6,3F12.5,1P,2E10.2)                                   
!  Loop 175 is 6 iterations on first sequence and two iterations theraft
      do 175 m=1,ni 
!  Loop 174 is for each substep within a sequence.                      
      do 174 j=2,8 
      jd=j-1 
      jdm=j-2 
      s=h(j) 
      q=s 
      if(ncl) q=one 
!  Use Eqs. (2.9) and (2.10) of text to predict positions at each aubste
!  These collapsed series are broken into two parts because an otherwise
!  excellent  compiler could not handle the complicated expression.     
      do 130 k=1,nv 
      a=w(3)*b(3,k)+s*(w(4)*b(4,k)+s*(w(5)*b(5,k)+s*(w(6)*b(6,k)+       &
     &   s*w(7)*b(7,k))))                                               
      y(k)=x(k)+q*(t*v(k)+t2*s*(f1(k)*w1+s*(w(1)*b(1,k)+s*(w(2)*b(2,k)  &
     &  +s*a))))                                                        
      if(npq) go to 130 
!  Next are calculated the velocity predictors need for general class II
      a=u(3)*b(3,k)+s*(u(4)*b(4,k)+s*(u(5)*b(5,k)+s*(u(6)*b(6,k)+       &
     &    s*u(7)*b(7,k))))                                              
      z(k)=v(k)+s*t*(f1(k)+s*(u(1)*b(1,k)+s*(u(2)*b(2,k)+s*a))) 
  130 continue 
!  Find forces at each substep.                                         
      call forceb(y,z,tstart+tm+s*t,fj) 
      nf=nf+1 
      do 171 k=1,nv 
!  Find G-value for the force FJ found at the current substep. This     
!  section, including the many-branched GOTO, uses Eq. (2.4) of text.   
      temp=g(jd,k) 
      gk=(fj(k)-f1(k))/s 
      go to (102,102,103,104,105,106,107,108),j 
  102 g(1,k)=gk 
      go to 160 
  103 g(2,k)=(gk-g(1,k))*r(1) 
      go to 160 
  104 g(3,k)=((gk-g(1,k))*r(2)-g(2,k))*r(3) 
      go to 160 
  105 g(4,k)=(((gk-g(1,k))*r(4)-g(2,k))*r(5)-g(3,k))*r(6) 
      go to 160 
  106 g(5,k)=((((gk-g(1,k))*r(7)-g(2,k))*r(8)-g(3,k))*r(9)-             &
     &     g(4,k))*r(10)                                                
      go to 160 
  107 g(6,k)=(((((gk-g(1,k))*r(11)-g(2,k))*r(12)-g(3,k))*r(13)-         &
     &     g(4,k))*r(14)-g(5,k))*r(15)                                  
      go to 160 
  108 g(7,k)=((((((gk-g(1,k))*r(16)-g(2,k))*r(17)-g(3,k))*r(18)-        &
     &     g(4,k))*r(19)-g(5,k))*r(20)-g(6,k))*r(21)                    
!  Upgrade all B-values.                                                
  160 temp=g(jd,k)-temp 
      b(jd,k)=b(jd,k)+temp 
!  TEMP is now the improvement on G(JD,K) over its former value.        
!  Now we upgrade the B-value using this dfference in the one term.     
!  This section is based on Eq. (2.5).                                  
      go to (171,171,203,204,205,206,207,208),j 
  203 b(1,k)=b(1,k)+c(1)*temp 
      go to 171 
  204 b(1,k)=b(1,k)+c(2)*temp 
      b(2,k)=b(2,k)+c(3)*temp 
      go to 171 
  205 b(1,k)=b(1,k)+c(4)*temp 
      b(2,k)=b(2,k)+c(5)*temp 
      b(3,k)=b(3,k)+c(6)*temp 
      go to 171 
  206 b(1,k)=b(1,k)+c(7)*temp 
      b(2,k)=b(2,k)+c(8)*temp 
      b(3,k)=b(3,k)+c(9)*temp 
      b(4,k)=b(4,k)+c(10)*temp 
      go to 171 
  207 b(1,k)=b(1,k)+c(11)*temp 
      b(2,k)=b(2,k)+c(12)*temp 
      b(3,k)=b(3,k)+c(13)*temp 
      b(4,k)=b(4,k)+c(14)*temp 
      b(5,k)=b(5,k)+c(15)*temp 
      go to 171 
  208 b(1,k)=b(1,k)+c(16)*temp 
      b(2,k)=b(2,k)+c(17)*temp 
      b(3,k)=b(3,k)+c(18)*temp 
      b(4,k)=b(4,k)+c(19)*temp 
      b(5,k)=b(5,k)+c(20)*temp 
      b(6,k)=b(6,k)+c(21)*temp 
  171 continue 
  174 continue 
      if(nes.or.m.lt.ni) go to 175 
!  Integration of sequence is over. Next is sequence size control.      
      hv=zero 
      do 635 k=1,nv 
  635 hv=dmax1(hv,dabs(b(7,k))) 
      hv=hv*w(7)/tval**7 
  175 continue 
      if (nsf) go to 180 
      if(.not.nes) tp=(ss/hv)**pw*dir 
      if(nes) tp=xl 
      if(nes) go to 170 
      if(tp/t.gt.one) go to 170 
    8 format (2x,2i2,2d18.10) 
      tp=.8d0*tp 
      ncount=ncount+1 
      IF(ncount.gt.10) THEN 
          WRITE(*,*)' ra15b: integration failed' 
          stop 
      ENDIF 
!     IF(NCOUNT.GT.1) WRITE (4,8) NOR,NCOUNT,T,TP                       
!  Restart with 0.8x sequence size if new size called for is smaller tha
!  originally chosen starting sequence size on first sequence.          
      go to 4000 
  170 nsf=.true. 
! Loop 35 finds new X and V values at end of sequence using Eqs. (2.11),
  180 do 35 k=1,nv 
      x(k)=x(k)+v(k)*t+t2*(f1(k)*w1+b(1,k)*w(1)+b(2,k)*w(2)+b(3,k)*w(3) &
     &    +b(4,k)*w(4)+b(5,k)*w(5)+b(6,k)*w(6)+b(7,k)*w(7))             
      if(ncl) go to 35 
      v(k)=v(k)+t*(f1(k)+b(1,k)*u(1)+b(2,k)*u(2)+b(3,k)*u(3)            &
     &    +b(4,k)*u(4)+b(5,k)*u(5)+b(6,k)*u(6)+b(7,k)*u(7))             
   35 continue 
      tm=tm+t 
      ns=ns+1 
!  Return if done.                                                      
      if(.not.nper) go to 78 
!     WRITE(*,7) NF,NS,X(1),X(2),T,TM,TF                                
!     WRITE(4,7) NF,NS                                                  
      DO k=1,nv 
        xx(k)=x(k) 
        vv(k)=v(k) 
      ENDDO 
      return 
!  Control on size of next sequence and adjust last sequence to exactly 
!  cover the integration span. NPER=.TRUE. set on last sequence.        
   78 call forceb (x,v,tstart+tm,f1) 
      nf=nf+1 
      if(nes) go to 341 
      tp=dir*(ss/hv)**pw 
      if(tp/t.gt.sr) tp=t*sr 
  341 if(nes) tp=xl 
      if(dir*(tm+tp).lt.dir*tf-1.d-8) go to 77 
      tp=tf-tm 
      nper=.true. 
!  Now predict B-values for next step. The predicted values from the pre
!  sequence were saved in the E-matrix. Te correction BD between the act
!  B-values found and these predicted values is applied in advance to th
!  next sequence. The gain in accuracy is significant. Using Eqs. (2.13)
   77 q=tp/t 
      do 39 k=1,nv 
      if(ns.eq.1) go to 31 
      do 20 j=1,7 
   20 bd(j,k)=b(j,k)-e(j,k) 
   31 e(1,k)=      q*(b(1,k)+ 2.d0*b(2,k)+ 3.d0*b(3,k)+                 &
     &           4.d0*b(4,k)+ 5.d0*b(5,k)+ 6.d0*b(6,k)+ 7.d0*b(7,k))    
      e(2,k)=                q**2*(b(2,k)+ 3.d0*b(3,k)+                 &
     &           6.d0*b(4,k)+10.d0*b(5,k)+15.d0*b(6,k)+21.d0*b(7,k))    
      e(3,k)=                             q**3*(b(3,k)+                 &
     &           4.d0*b(4,k)+10.d0*b(5,k)+20.d0*b(6,k)+35.d0*b(7,k))    
      e(4,k)=   q**4*(b(4,k)+ 5.d0*b(5,k)+15.d0*b(6,k)+35.d0*b(7,k)) 
      e(5,k)=                q**5*(b(5,k)+ 6.d0*b(6,k)+21.d0*b(7,k)) 
      e(6,k)=                             q**6*(b(6,k)+ 7.d0*b(7,k)) 
      e(7,k)=                                           q**7*b(7,k) 
      do 39 l=1,7 
   39 b(l,k)=e(l,k)+bd(l,k) 
!  Two iterations for every sequence after the first.                   
      ni=2 
      go to 722 
      END subroutine ra15b
