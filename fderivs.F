#include "fintrf.h"
C======================================================================
#if 0
C     
C     fderivs.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     fderivs.f
C
C     Finds df_{ab}(p,q)/dt.   
C======================================================================
C     The gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      implicit none

C     in/out variables:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
C     Pointers to input/output mxArrays:
      mwPointer F_ptr, p_ptr, q_ptr, a_ptr, b_ptr, l_ptr, c_ptr, df_ptr
C     Input arguments for computational routine:
      integer p, q, alpha, beta
      real*8  F(501,501,3,3), lambda, cnu, df
C     Input arguments for common block:
      mwPointer z_ptr, n_ptr, t_ptr, e_ptr, lm_ptr
      integer Z, N
      real*8 taue,extdot,lambdam
C     Fixed values for common block:
      real*8 retshift, ds, dssqinv, ds2inv
      real*8 epsilon, eps2inv, pi, lambdam2, const
C     F size variables:
      mwPointer mrows, ncols
      mwSize size

C     Function declarations:
      mwPointer mxGetPr, mxGetM, mxGetN
      mwPointer mxCreateDoubleMatrix

C     Declare common library:
      common/retractionshift/RetShift
      common/segment/ds,dssqinv,ds2inv
      common/epsilon/epsilon,eps2inv
      common/taue/taue       
      common/Pi/pi 
      common/extension/extdot
      common/entanglements/Z
      common/points/N
      common/lambdamaxsq/lambdam2
      common/constant/const

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 12) then
         call mexErrMsgIdAndTxt ('MATLAB:fderivsT:nInput',
     +                           'Twelve inputs required.')
      elseif(nlhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fderivsT:nOutput',
     +                           'One output required.')
      endif

C     Get the size of the input array.
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      size = mrows*ncols

C     Create Fortran variables from the input arguments.
      F_ptr = mxGetPr(prhs(1))
      p_ptr = mxGetPr(prhs(2))
      q_ptr = mxGetPr(prhs(3))
      a_ptr = mxGetPr(prhs(4))
      b_ptr = mxGetPr(prhs(5))
      l_ptr = mxGetPr(prhs(6))
      c_ptr = mxGetPr(prhs(7))
      call mxCopyPtrToReal8(F_ptr,F,size)
      call mxCopyPtrToInteger8(p_ptr,p,1)
      call mxCopyPtrToInteger8(q_ptr,q,1)
      call mxCopyPtrToInteger8(a_ptr,alpha,1)
      call mxCopyPtrToInteger8(b_ptr,beta,1)
      call mxCopyPtrToReal8(l_ptr,lambda,1)
      call mxCopyPtrToReal8(c_ptr,cnu,1)

C     Save inputs in common library.
      t_ptr = mxGetPr(prhs(8))
      e_ptr = mxGetPr(prhs(9))
      z_ptr = mxGetPr(prhs(10))
      n_ptr = mxGetPr(prhs(11))
      lm_ptr = mxGetPr(prhs(12))
      call mxCopyPtrToReal8(t_ptr,taue,1)
      call mxCopyPtrToReal8(e_ptr,extdot,1)
      call mxCopyPtrToInteger8(z_ptr,Z,1)
      call mxCopyPtrToInteger8(n_ptr,N,1)
      call mxCopyPtrToReal8(lm_ptr,lambdam,1)

C     Set up common library.
      pi = 3.14159265359
      RetShift = 2.0
      ds=(1.0*Z)/(1.0*N)
      dssqinv = 1.0/(ds**2)
      ds2inv = 0.5/ds
      epsilon = 0.001
      eps2inv = 0.5/epsilon
      lambdam2 = lambdam*lambdam
      const = (lambdam2 - 1.0)/(lambdam2 - 1.0/3.0)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
      df_ptr = mxGetPr(plhs(1))

C     Call the computational subroutine.
      call derivs(p,q,F,alpha,beta,lambda,cnu,df)

C     Load the data into y_pr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(df, df_ptr, 1)

      return
      end
C-----------------------------------------------------------------------
C     Computational subroutine
      subroutine derivs(p,q,F,alpha,beta,lambda,cnu,df)
        implicit none

        integer p,q,alpha,beta,mini
        integer N
        integer nmax

        parameter(nmax=500)

        double precision F(0:nmax,0:nmax,3,3)
        double precision cnu,df,retshift,ds,nu,lambda
        double precision sum,sumret,sumcon,sumccr,sumrep
        double precision pi,taue, extdot,dssqinv,ds2inv
        double precision pp,pm,qp,qm,epsilon,eps2inv
        double precision fpq,fpp1q,fpm1q,fpqp1,fpqm1
        double precision fpp1qp1,fpm1qm1
        double precision feqpq,feqpp1q,feqpm1q,feqpqp1,feqpqm1
        double precision trpp,trqq
        double precision trppp1,trppm1,trqqp1,trqqm1
        double precision trppsr,trqqsr
        double precision trppp1sr,trppm1sr,trqqp1sr,trqqm1sr
        double precision flampp,flamqq
        double precision flamppp1,flamppm1,flamqqp1,flamqqm1
        double precision flambylpp,flambylqq
        double precision flambylppp1,flambylppm1
        double precision flambylqqp1,flambylqqm1
        double precision trmsr,trmp1sr,trmm1sr,dclfpq

        double precision TrF,Feq,Dclf,flambda,flambdabyl
	    integer findMini
        external TrF,Feq,Dclf,flambda,flambdabyl, findMini

        common/retractionshift/RetShift
        common/segment/ds,dssqinv,ds2inv
        common/epsilon/epsilon,eps2inv
        common/taue/taue       
        common/Pi/pi 
        common/extension/extdot
        common/points/N

        nu=Cnu*lambda
        sum=0.0
        sumret = 0.0
        sumcon = 0.0
        sumccr = 0.0
        sumrep = 0.0

        fpq = f(p,q,alpha,beta)
        fpp1q = f(p+1,q,alpha,beta)
        fpm1q = f(p-1,q,alpha,beta)
        fpqp1 = f(p,q+1,alpha,beta)
        fpqm1 = f(p,q-1,alpha,beta)

        feqpq = feq(p,q,alpha,beta)
        feqpp1q = feq(p+1,q,alpha,beta)
        feqpm1q = feq(p-1,q,alpha,beta)
        feqpqp1 = feq(p,q+1,alpha,beta)
        feqpqm1 = feq(p,q-1,alpha,beta)

        trpp = TrF(p,p,F)
        trqq = TrF(q,q,F)
        trppp1 = TrF(p+1,p+1,F)
        trppm1 = TrF(p-1,p-1,F)
        trqqp1 = TrF(q+1,q+1,F)
        trqqm1 = Trf(q-1,q-1,F)
        
        trppsr = dsqrt(TrF(p,p,F))
        trqqsr = dsqrt(TrF(q,q,F))
        trppp1sr = dsqrt(TrF(p+1,p+1,F))
        trppm1sr = dsqrt(TrF(p-1,p-1,F))
        trqqp1sr = dsqrt(TrF(q+1,q+1,F))
        trqqm1sr = dsqrt(Trf(q-1,q-1,F))

        flampp = flambda(p,p,F)
        flamqq = flambda(q,q,F)
        flamppp1 = flambda(p+1,p+1,F)
        flamppm1 = flambda(p-1,p-1,F)
        flamqqp1 = flambda(q+1,q+1,F)
        flamqqm1 = flambda(q-1,q-1,F)

        flambylpp = flambdabyl(p,p,F)
        flambylqq = flambdabyl(q,q,F)
        flambylppp1 = flambdabyl(p+1,p+1,F)
        flambylppm1 = flambdabyl(p-1,p-1,F)
        flambylqqp1 = flambdabyl(q+1,q+1,F)
        flambylqqm1 = flambdabyl(q-1,q-1,F)

C     =======Retraction term======================
   
C       ###  (df/dp)*(1/lambda_p)*(dflambda/dp) term ###

        sumret=sumret+
     &     ds2inv*(fpp1q-fpm1q)
     &    *1.0/trppsr
     &    *ds2inv*(flamppp1-flamppm1) 
									    
        sumret=sumret+
     &     ds2inv*(fpqp1-fpqm1)
     &    *1.0/trqqsr
     &    *ds2inv*(flamqqp1-flamqqm1)  

									
C       ###  f*d(1/lambda_p)/dp*(dflambda/dp) term ###

        sumret=sumret+fpq
     &    *ds2inv*(1.0/trppp1sr-1.0/trppm1sr)
     &    *ds2inv*(flamppp1-flamppm1)

        sumret=sumret+fpq
     &   *ds2inv*(1.0/trqqp1sr-1.0/trqqm1sr)
     &   *ds2inv*(flamqqp1-flamqqm1)  


C       ###  f*(1/lambda_p)*(d2/dp^2)flambda term  ###

        sumret=sumret+fpq
     &   *1.0/trppsr
     &   *dssqinv*(flamppp1+flamppm1-2.0*flampp)     

        sumret=sumret+fpq
     &   *1.0/trqqsr
     &   *dssqinv*(flamqqp1+flamqqm1-2.0*flamqq)

        
        sumret=sumret/(Pi**2*taue)*RetShift

C     =======Convection term======================

C=====For extension=====
C       ###  Fxx  ###
        if (alpha.EQ.1) then
	       if (beta.EQ.1) then
              sumcon=sumcon+2.0*extdot*F(p,q,1,1)
	       endif
	    endif

C       ###  Fyy  ###
        if (alpha.EQ.2) then
           if (beta.EQ.2) then
              sumcon=sumcon-extdot*F(p,q,2,2)
           endif
        endif

C=====For shear=====
C       ###  Fxx  ###
        !if (alpha.EQ.1) then
        !   if (beta.EQ.1) then
        !      sumcon=sumcon+2.0*extdot*F(p,q,1,2)
        !   endif
        !endif

C       ###  Fxy,yx term ###
        !if (alpha.EQ.1) then
        !   if (beta.EQ.2) then
        !      sumcon=sumcon+extdot*F(p,q,2,2)
        !   endif
        !endif

C     =======CCr term======================
 
C       ###   f(s,s') terms   ###

C       ###   (Flambda/(lambdap2))*(d^2/ds^2)F term   ###

        sumccr=sumccr+
     &     dssqinv*(fpp1q+fpm1q-2.0*fpq)
     &     *flampp/Trpp
        
        sumccr=sumccr+
     &     dssqinv*(fpqp1+fpqm1-2.0*fpq)
     &     *flamqq/Trqq


C       ###   flambda*(d(1/lambdap2)/ds)*d/ds(F)   ###

        sumccr=sumccr+flampp
     &     *ds2inv*(1.0/Trppp1-1.0/Trppm1) 
     &     *ds2inv*(fpp1q-fpm1q)
    
        sumccr=sumccr+flamqq
     &     *ds2inv*(1.0/Trqqp1-1.0/Trqqm1) 
     &     *ds2inv*(fpqp1-fpqm1)


C       ###   d(flambda)/ds*(1/lambdap2)*d/ds(F)   ###

        sumccr=sumccr+
     &      ds2inv*(flamppp1-flamppm1)
     &     *1.0/Trpp 
     &     *ds2inv*(fpp1q-fpm1q)

        sumccr=sumccr+
     &      ds2inv*(flamqqp1-flamqqm1)
     &     *1.0/Trqq
     &     *ds2inv*(fpqp1-fpqm1)


C       ###   Extra terms  ###

C       ###   f/lambdap*d2(flambdabyl)/dp2

        sumccr=sumccr+fpq/trppsr
     &     *dssqinv*(flambylppp1+flambylppm1-2.0*flambylpp)

        sumccr=sumccr+fpq/trqqsr
     &     *dssqinv*(flambylqqp1+flambylqqm1-2.0*flambylqq)
                

C       ###  (df/ds)*(1/lambdap)*(dflambdabyl/dp)

        sumccr=sumccr+1.0/trppsr
     &     *ds2inv*(fpp1q-fpm1q)
     &     *ds2inv*(flambylppp1-flambylppm1)

        sumccr=sumccr+1.0/trqqsr
     &     *ds2inv*(fpqp1-fpqm1)
     &     *ds2inv*(flambylqqp1-flambylqqm1)


C       ###  f*(d(1/lambdap)/dp)*(dflambdabyl/dp)

        sumccr=sumccr+fpq
     &   *ds2inv*(1.0/trppp1sr-1.0/trppm1sr)
     &   *ds2inv*(flambylppp1-flambylppm1)

        sumccr=sumccr+fpq
     &   *ds2inv*(1.0/trqqp1sr-1.0/trqqm1sr)
     &   *ds2inv*(flambylqqp1-flambylqqm1)


C       ###   feq terms   ###

C       ###   1/sqrt(TrFpp)(d^2/ds^2)Feq term   ###

        sumccr=sumccr-
     &     dssqinv*(Feqpp1q+Feqpm1q-2.0*Feqpq)
     &     *1.0/trppsr
     
        sumccr=sumccr-
     &     dssqinv*(Feqpqp1+Feqpqm1-2.0*Feqpq)
     &     *1.0/trqqsr


C       ###   d/ds(1/sqrt(TrFpp)*d/ds(Feq)   ###

        sumccr=sumccr-
     &     ds2inv*(1.0/trppp1sr-1.0/trppm1sr)
     &     *ds2inv*(Feqpp1q-Feqpm1q)
     
        sumccr=sumccr-
     &     ds2inv*(1.0/trqqp1sr-1.0/trqqm1sr) 
     &     *ds2inv*(Feqpqp1-Feqpqm1)


        sumccr=sumccr*1.5*nu

C     =======Reptation + CLF term======================
 
C       find point closest to chain end
        mini=findMini(p,q)

        pp=ds*p+epsilon
        pm=ds*p-epsilon 
        qp=ds*q+epsilon
        qm=ds*q-epsilon
 
        fpp1qp1 = f(p+1,q+1,alpha,beta)
        fpm1qm1 = f(p-1,q-1,alpha,beta)

        trmsr = dsqrt(trF(mini,mini,F))
        trmp1sr = dsqrt(TrF(mini+1,mini+1,F))
        trmm1sr = dsqrt(TrF(mini-1,mini-1,F))

        dclfpq = Dclf(ds*p,ds*q)

C       ###   (d/dp+d/dq)Dclf*(1/sqrt(Trfmin))*(d/dp+d/dq)F   ###

        sumrep=sumrep+
     &     eps2inv*(Dclf(pp,qp)-Dclf(pm,qm))
     &     *1.0/trmsr
     &     *ds2inv*(fpp1qp1-fpm1qm1)


C       ###  Dclf*(d/dp+d/dq)(1/sqrt(Trfmin))*(d/dp+d/dq)F

        sumrep=sumrep+
     &     Dclfpq
     &     *ds2inv*(1.0/trmp1sr-1.0/trmm1sr) 
     &     *ds2inv*(fpp1qp1-fpm1qm1)


C       ###  Dclf*(1/sqrt(Trfmin))*(d/dq+d/dq)**2F  ###

        sumrep=sumrep+
     &     Dclfpq/trmsr
     &     *dssqinv*(fpp1qp1+fpm1qm1-2.0*fpq)

        sumrep=sumrep*1.0/(3.0*Pi**2*taue)*(1.0/trmsr)

C=====Output df=====

      df = sumret + sumcon + sumccr + sumrep
      return
      end

C=====TrF function=====

      Double Precision function TrF(p,q,F)
        implicit none

        integer p,q
        integer nmax

        parameter(nmax=500)

        double precision  F(0:nmax,0:nmax,3,3)

      TrF=F(p,q,1,1)+2.0*F(p,q,2,2)
      end function

C=====flambda function=====
        
      Double precision function flambda(p,q,F)
        implicit none
  
        integer p,q
        integer nmax
        double precision term,lambdam2,const,trace

        double precision TrF
        external TrF

        common/lambdamaxsq/lambdam2
        common/constant/const
        
        parameter(nmax=500)

        double precision F(0:nmax,0:nmax,3,3)
 
      trace = TrF(p,q,F)
      if (trace .gt. lambdam2) then
         write(*,*) "Lambda exceeding lambda_max !!!!!!"
         write(*,*) lambdam2,trace
         stop
      endif
      term=(lambdam2-trace/3.0)/(lambdam2-trace)
      term=dsqrt(trace)*term
      flambda=term*const
      end function

C=====flambdabyl function=====

      Double precision function flambdabyl(p,q,F)
        implicit none
  
        integer p,q
        integer nmax
        double precision term,lambdam2,const,trace

        double precision TrF
        external TrF

        common/lambdamaxsq/lambdam2
        common/constant/const
        
        parameter(nmax=500)

        double precision F(0:nmax,0:nmax,3,3)
 
      trace = TrF(p,q,F)
      term=(lambdam2-trace/3.0)/(lambdam2-trace)
      flambdabyl=term*const
      end function        

C=====Feq function=====

      Double Precision function Feq(p,q,alpha,beta)
        implicit none

        integer p,q,alpha,beta
        double precision dum,ds,preal,qreal,dssqinv,ds2inv

        common/segment/ds,dssqinv,ds2inv

      dum=0.0
      preal = 1.0*p
      qreal = 1.0*q
      if(dabs(preal-qreal).lt.ds2inv) then
         if(alpha.EQ.beta) then
            dum=1.0/3.0
         endif
      endif
      Feq=dum
      end function

C=====Dclf function====(p,q inputs here are double precision!)
       
      Double Precision Function Dclf(p,q)
        implicit none

        Double Precision  alphad,s,p,q,zreal
        integer z
        
        common/entanglements/Z

      zreal = 1.0*z
      alphad=1.15
      s=Min(1.0*p,1.0*q,1.0*(Z-p),1.0*(Z-q))
 
      if (s.lt.alphaD*dsqrt(zreal)) then
         Dclf=(alphaD/s)**2
      else
         if(s.gt.zreal-alphaD*dsqrt(zreal)) then
            Dclf=(alphaD/(s-zreal))**2
         else
            Dclf=1.0/zreal
         endif
      endif

      end function

C=====findMini function=====

      Integer function findMini(p,q)

	    implicit none
	    integer :: p,q, answ,N

  	    common/points/N
	
	  if(p==q) then
	     answ=p
	  endif
	
	  if( (1.0*p< 1.0*N/2.0) .AND. (1.0*q<1.0*N/2.0)) then
	     answ = min( p,q)
	  endif
	
	  if( (1.0*p > 1.0*N/2.0) .AND. (1.0*q > 1.0*N/2.0)) then
	     answ = max(p,q)
	  endif
	
	  if( (1.0*p < 1.0*N/2.0) .AND. (1.0*q > 1.0*N/2.0)) then
	     if( p< N-q) then
	        answ = p
	     else
	        answ = q
	     endif
	  endif
	
	  if( (1.0*p > 1.0*N/2.0) .AND. (1.0*q < 1.0*N/2.0)) then
	     if( q< N-p) then
	        answ = q
	     else
	        answ = p
	     endif
	  endif
	
	  if(p==q) then
	     answ = p
	  endif
	
	  if(p==N/2) then
	     answ = q
	  endif
	
	  if(q==N/2) then
	     answ = p
	  endif
	
	  !print*,"mini",p,q,N, N/2, answ
	
	  findMini = answ
	
	  end function findMini