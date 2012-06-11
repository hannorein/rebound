      program main
      stop
      end

      subroutine RA15(x,v,tf,xl,ll,nv,nclass,nor)
      implicit real*8 (a-h,o-z)
      real*4 tval,pw
      dimension x(1),v(1),f1(18),fj(18),c(21),d(21),r(21),y(18),z(18),
     a     b(7,18),g(7,18),e(7,17),bd(7,18),h(8),w(7),u(7),nw(8)
      logical npq,nsf,nper,ncl,nes
      data nw/0,0,1,3,6,10,15,21/
      data zero,half,one,sr/0.0d0,0.5d0,1.0d0,1.4d0/
      data h/0.d0,0.05626256053692215d0,0.18024069173689236d0,
     a 0.35262471711316964d0,0.54715362633055538d0,0.73421017721541053,
     b 0.88532094683909577d0,0.97752061356128750d0/
      nper=.false.
      nsf=.false.
      ncl=nclass.eq.1
      npq=nclass.lt.2
      dir=one
      if(tf.lt.zero) dir=-one
      nes=ll.lt.0
      end
