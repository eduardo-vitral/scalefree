      PROGRAM FITVP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C Read a Velocity Profile from file. Then fiot a Gaussian and determine
C the Gauss-Hermite moments.
C
C This code was simplified/extracted from scalefree.f,
C cresated by RvdM around 1995-1996.      
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION harr(0:100), velar(1001),vpval(1001), weight(1001)
C
CCCCCCCCCCCCCCCCCCCC
C      
      COMMON /vpcur/ vvelar(1001),vvpval(1001),wweight(1001),nvvel
      CHARACTER(500) :: fileplace
C
C Common block in which the current velocity profile is stored, for use in
C REAL FUNCTION CHI2H.
C
CCCCCCCCCCCCCCCCCCCC
C
C Asks the user to specify a path for the input file
C     
      WRITE (*,*) 'Please specify a path where to find the input file.'
      WRITE (*,*) 'Expected file shape:'
      WRITE (*,*) '--> First line (header): "	v  VP(v)"'
      WRITE (*,*) '--> 1st column: x / 2nd column: f(x)'
      READ (*,'(A)') fileplace
      WRITE (*,*) ' '
C
C
CCCCCCCCCCCCCCCCCCCC      
C
C Read input file      
C      
      OPEN (UNIT=11,FILE=fileplace,STATUS='OLD')
      READ (11,*)
C
      i    = 0
C      
 51   READ (UNIT=11,FMT=*,ERR=52,END=52) vcur, vpcur
C
      i = i+1
      vvelar(i) = vcur 
      vvpval(i) = vpcur
C
      GOTO 51
C
C When we get here we are done reading the input file
C      
 52   CLOSE (UNIT=11)
      nvvel = i
C
      WRITE (*,*) 'Number of gridpoints read:'
      WRITE (*,*) nvvel
C
CCCCCCCCCCCCCCCCCCCC
C
C Set the weights for the points and calculate the mean and dispersion      
C
      wweight(1) = (vvelar(2)-vvelar(1))
      DO i=2,nvvel-1      
        wweight(i) = 0.5D0*(vvelar(i+1)-vvelar(i-1))
      END DO   
      wweight(nvvel) = (vvelar(nvvel)-vvelar(nvvel-1))
C
      sum0 = 0.0D0
      sum1 = 0.0D0
      sum2 = 0.0D0
C
      DO i=1,nvvel
        sum0 = sum0 + (vvpval(i)*wweight(i))
        sum1 = sum1 + (vvelar(i)*vvpval(i)*wweight(i))
        sum2 = sum2 + ((vvelar(i)**2)*vvpval(i)*wweight(i))
      END DO
C
      vmean = (sum1/sum0)
      vdisp = SQRT((sum2/sum0)-(vmean**2))
C
      WRITE (*,*) 'True normalization, mean, dispersion:'
      WRITE (*,'(3F12.4)') sum0, vmean, vdisp
C
CCCCCCCCCCCCCCCCCCCC
C
C Make a copy of the array
C
      nvel = nvvel
      DO i=1,nvel
        weight(i) = wweight(i)
        velar(i)  = vvelar(i)
        vpval(i)  = vvpval(i)
      END DO
C
CCCCCCCCCCCCCCCCCCCC
C
C Find the best fitting Gaussian (by searching for those values that
C come as closely as possible to generating h0=1, h1=h2=0).
C Set initial guesses before running.
C      
      gam  = sum0
      Vgau = vmean
      sig  = vdisp
C      
      CALL FITGAUSS (gam,Vgau,sig)
C
      WRITE (*,*) 'Best Gauss-fit normalization, mean, dispersion:'
      WRITE (*,'(3F12.4)') gam,Vgau,sig
C
CCCCCCCCCCCCCCCCCCCC
C      
C Now get all the Gauss-Hermite coefficients up to order nord
C
      nord = 10
      CALL GAUHERM (velar,vpval,weight,nvel,gam,Vgau,sig,harr,nord)
C
      WRITE (*,*) 'GH moments order 0-10'
      WRITE (*,'(11F7.3)') (harr(k),k=0,10)
C      
CCCCCCCCCCCCCCCCCCCC
C      
      END


      SUBROUTINE FITGAUSS (gam,Vgau,sig)
C
C Find the best fitting Gaussian parameters (gam,Vgau,sig) for
C the VP in the common block /vpcur/. The parameters must be preset at
C initial guesses. The function CHI2H below is minimized using the Numerical
C Recipes routine amoeba.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (eps = 1.0D-16)
C
      DIMENSION P(4,3),Y(4),help(3)
C
C Starting simpleces for routine AMOEBA, and the function values
C in the starting simpleces.
C
      EXTERNAL CHI2H
C
C Initialize the fits
C
      epsmal = 0.6D0
      epsmal = MIN(0.9D0,ABS(epsmal))
      pl1    = 1.0D0 + epsmal
      xmn1   = 1.0D0 - epsmal
C
C Determines the size of the starting simplex
C
      P(1,1) = gam
      P(1,2) = Vgau
      P(1,3) = sig * pl1
      P(2,1) = gam * xmn1
      P(2,2) = Vgau - (sig*epsmal)
      P(2,3) = sig * xmn1
      P(3,1) = gam * pl1
      P(3,2) = Vgau - (sig*epsmal)
      P(3,3) = sig * xmn1
      P(4,1) = gam
      P(4,2) = Vgau + (sig*epsmal)
      P(4,3) = sig * xmn1
C
C Initialize
C
      DO i=1,4
        DO j=1,3
          help(j) = P(i,j)
        END DO
        Y(i) = CHI2H(help)
      END DO
C
      CALL AMOEBA(P,Y,4,3,3,eps,CHI2H,iter)
C
      gam  = (ABS(P(1,1))+ABS(P(2,1))+ABS(P(3,1))+
     &        ABS(P(4,1)))/4.0D0
      Vgau = (P(1,2)+P(2,2)+P(3,2)+P(4,2))/4.0D0
      sig  = (ABS(P(1,3))+ABS(P(2,3))+ABS(P(3,3))+
     &        ABS(P(4,3)))/4.0D0
C
      END


      REAL*8 FUNCTION CHI2H (y)
C
C Calculates the chih^2 = (h0-1)^2 + (h1^2) + (h2^2) for a Gaussian
C with parameters gam = |y(1)|, Vgau = y(2), sig = |y(3)|, 
C for the VP in the common block /vpcur/ 
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION y(3),harr(0:100)
C
      COMMON /vpcur/ velar(1001),vpval(1001),weight(1001),nvel
C
C Avoid values of gamma and sigma too close to zero. 
C
      gam  = MAX(1.0D-3,ABS(y(1)))
      Vgau = y(2)
      sig  = MAX(1.0D-3,ABS(y(3)))
C
      CALL GAUHERM (velar,vpval,weight,nvel,gam,Vgau,sig,harr,2)
C
      CHI2H = 1.0D0 + ((harr(0)-1.0D0)**2.0D0) + 
     &        (harr(1)**2.0D0) + (harr(2)**2.0D0)
C
      END


      SUBROUTINE GAUHERM (velar,vpval,weight,nvel,
     &                    gam,Vgau,sig,harr,nhord)
C
C Given the VP as an array vpval, corrsponding to nvel velocities velar,
C calculate the Gauss-Hermite coefficients harr up to order nhord that belong
C to the given reference values gam,Vgau,sig.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (pi=3.14159265358979D0)
C
      DIMENSION velar(1001),vpval(1001),weight(1001),harr(0:100)
C
      DO l=0,nhord
        harr(l) = 0.0D0
      END DO
C
      DO i=1,nvel
        w     = (velar(i)-Vgau)/sig
        DO l=0,nhord
          harr(l) = harr(l) + 
     &      (vpval(i)*weight(i)*SDGAUSS(w)*H_POL(l,w))
        END DO
      END DO
C
      DO l=0,nhord
        harr(l) = harr(l) * 2.0D0 * SQRT(pi) / gam 
      END DO
C
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The function SDGAUSS, HE_POL and H_POL calculate the standard
C Gaussian and Hermite polynomials.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION SDGAUSS (x)
C
C Returns the standard Gaussian as function of x
C
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (pi=3.14159265358979D0)
      SDGAUSS = (1.0D0/SQRT(2.0D0*pi)) * EXPP(-0.5D0*x*x)
      END


      REAL*8 FUNCTION HE_POL (l,x)
C
C Returns the value of the Hermite polynomial He_l(x) as defined
C in Appendix A of van der Marel & Franx.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      HE_POL = 0.0D0
      dl = DBLE(l)
      DO j=l,0,-1
        IF (2*((j+l)/2).EQ.j+l) THEN
          dj = DBLE(j)
          pfacln = (0.5D0*gammaln(dl+1.0D0)) - 
     &             gammaln(dj+1.0D0) - gammaln(dl-dj+1.0D0) +
     &             gammaln(0.5D0*(dl-dj+1.0D0)) - 
     &             gammaln(0.5D0) + (0.5D0*(dl-dj)*LOG(2.0D0))
          pfac = ((-1.0D0)**((l-j)/2)) * EXPP(pfacln)
        ELSE
          pfac = 0.0D0
        END IF
        HE_POL = pfac + (x*HE_POL)
      END DO
C
      END


      REAL*8 FUNCTION H_POL (l,x)
C
C Returns the value of the Hermite polynomial H_l(x) as defined
C in Appendix A of van der Marel & Franx.
C
      IMPLICIT REAL*8 (a-h,o-z)
      H_POL = HE_POL(l,x*SQRT(2.0D0))
      END

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The function EXPP is used throughout rather than the normal function
C EXP, to avoid underflow.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION EXPP(x)
C
C Modified exponential function that will not underflow
C
      IMPLICIT REAL*8 (a-h,o-z)
      IF (x.GE.-500D0) THEN      
        EXPP = EXP(x)
      ELSE
        EXPP = 0.0D0
      END IF
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The functions gammaln, betaln, pochln and binomln calculate (the ln of) 
C mathematical expressions that are often needed. They are all based 
C indirectly on the function gammln in Numerical recipes which is accurate
C approximately to 1.0D-10. 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      REAL*8 FUNCTION gammaln(x)
C
C The logarithm of the gamma function
C
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (pi=3.14159265358979D0)
      IF (x.GE.1.0D0) THEN
        gammaln = gammln(x)
      ELSE IF (x.GE.0.0D0) THEN
        z = 1.0D0 - x
        gammaln = (LOG((pi*z)/SIN(pi*z))) - gammln(2.0D0-x)
      ELSE
        STOP 'x < 0 in gammaln'
      END IF
      END

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Required routines from Numerical recipes. All were modified to
C be double precision. The parameter eps in QROMO determines the
C speed of the program in calculating the VP parameters. Best results
C are obtained with eps=1.0D-8, but this leads to a slow program.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      REAL*8 FUNCTION gammln(x)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER j
      DIMENSION cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END


      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL*8 ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=5000)
      EXTERNAL funk
CU    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL*8 rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1     do 12 n=1,ndim
        sum=0.0D0
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
C
C Possibility to write intermediate results
C
C      WRITE (*,'(I5,4F15.8)') iter,(p(ilo,j),j=1,3),y(ilo)
C
      if (iter.ge.ITMAX) then
        write (*,*) 'ITMAX exceeded in amoeba'
        return
      endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0D0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0D0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5D0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5D0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END


      REAL*8 FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL*8 fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
CU    USES funk
      INTEGER j
      REAL*8 fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.0D0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END


