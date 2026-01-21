      PROGRAM FITVP_STDOUT
C
C Read a Velocity Profile from file, fit a Gaussian and compute
C Gauss-Hermite moments. Emits results to STDOUT in a structured way.
C
C INPUT (stdin):
C   One line: path to the input file containing 2 columns (v, VP(v)).
C   A header line is allowed; it will be ignored if non-numeric.
C
C OUTPUT (stdout):
C   # gauss_moments norm <...> mean <...> dispersion <...>
C   # gauss_fit     norm <...> mean <...> dispersion <...>
C   # gh_moments 0-10 <h0> <h1> ... <h10>
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER ios, i, nvel, nvvel, nord
      CHARACTER*500 fileplace
      CHARACTER*500 line
      REAL*8 vcur, vpcur
      REAL*8 sum0, sum1, sum2
      REAL*8 vmean, vdisp, gam, Vgau, sig
      REAL*8 harr(0:100), velar(1001), vpval(1001), weight(1001)
      REAL*8 vvelar(1001), vvpval(1001), wweight(1001)
C
      COMMON /vpcur/ vvelar, vvpval, wweight, nvvel
C
C --- Read filepath from stdin (single line)
      READ (*,'(A)') fileplace
C
C --- Open and read file; allow a header line
      OPEN (UNIT=11, FILE=fileplace, STATUS='OLD', IOSTAT=ios)
      IF (ios .NE. 0) THEN
        WRITE (*,*) '# error could_not_open_file'
        STOP 1
      END IF
C
      nvvel = 0
C
C Read first line as text and try to parse numbers from it
      READ (11,'(A)',END=900) line
      READ (line,*,IOSTAT=ios) vcur, vpcur
      IF (ios .EQ. 0) THEN
        nvvel = nvvel + 1
        vvelar(nvvel) = vcur
        vvpval(nvvel) = vpcur
      END IF
C
C Read remaining lines as numeric
  50  CONTINUE
      READ (11,*,IOSTAT=ios,END=900) vcur, vpcur
      IF (ios .NE. 0) GOTO 900
      nvvel = nvvel + 1
      IF (nvvel .GT. 1001) GOTO 900
      vvelar(nvvel) = vcur
      vvpval(nvvel) = vpcur
      GOTO 50
C
 900  CONTINUE
      CLOSE (UNIT=11)
C
      IF (nvvel .LT. 3) THEN
        WRITE (*,*) '# error too_few_points'
        STOP 2
      END IF
C
C --- Weights and raw moments (normalization/mean/dispersion)
      wweight(1) = (vvelar(2)-vvelar(1))
      DO i=2,nvvel-1
        wweight(i) = 0.5D0*(vvelar(i+1)-vvelar(i-1))
      END DO
      wweight(nvvel) = (vvelar(nvvel)-vvelar(nvvel-1))
C
      sum0 = 0.0D0
      sum1 = 0.0D0
      sum2 = 0.0D0
      DO i=1,nvvel
        sum0 = sum0 + (vvpval(i)*wweight(i))
        sum1 = sum1 + (vvelar(i)*vvpval(i)*wweight(i))
        sum2 = sum2 + ((vvelar(i)**2)*vvpval(i)*wweight(i))
      END DO
      vmean = (sum1/sum0)
      vdisp = SQRT((sum2/sum0)-(vmean**2))
C
C Copy arrays into local working arrays
      nvel = nvvel
      DO i=1,nvel
        weight(i) = wweight(i)
        velar(i)  = vvelar(i)
        vpval(i)  = vvpval(i)
      END DO
C
C --- Best-fitting Gaussian (amoeba search on h0,h1,h2)
      gam  = sum0
      Vgau = vmean
      sig  = vdisp
      CALL FITGAUSS (gam, Vgau, sig)
C
C --- Compute GH moments (0..10) for the fitted Gaussian reference
      nord = 10
      CALL GAUHERM (velar, vpval, weight, nvel, gam, Vgau, sig, harr,
     &              nord)
C
C --- Structured stdout (machine-parseable; single-line records)
      WRITE(*,1001) sum0, vmean, vdisp
      WRITE(*,1002) gam,  Vgau,  sig
      WRITE(*,1003) (harr(i), i=0,10)

 1001 FORMAT('# gauss_moments norm ',1P,E25.16,' mean ',1P,E25.16,
     &       ' dispersion ',1P,E25.16)
 1002 FORMAT('# gauss_fit norm ',1P,E25.16,' mean ',1P,E25.16,
     &       ' dispersion ',1P,E25.16)
 1003 FORMAT('# gh_moments 0-10',11(1X,1P,E25.16))
C
      END


      SUBROUTINE FITGAUSS (gam,Vgau,sig)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (eps = 1.0D-16)
      DIMENSION P(4,3),Y(4),help(3)
      EXTERNAL CHI2H
      epsmal = 0.6D0
      epsmal = MIN(0.9D0,ABS(epsmal))
      pl1    = 1.0D0 + epsmal
      xmn1   = 1.0D0 - epsmal
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
      DO i=1,4
        DO j=1,3
          help(j) = P(i,j)
        END DO
        Y(i) = CHI2H(help)
      END DO
      CALL AMOEBA(P,Y,4,3,3,eps,CHI2H,iter)
      gam  = (ABS(P(1,1))+ABS(P(2,1))+ABS(P(3,1))+
     &        ABS(P(4,1)))/4.0D0
      Vgau = (P(1,2)+P(2,2)+P(3,2)+P(4,2))/4.0D0
      sig  = (ABS(P(1,3))+ABS(P(2,3))+ABS(P(3,3))+
     &        ABS(P(4,3)))/4.0D0
      END


      REAL*8 FUNCTION CHI2H (y)
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION y(3),harr(0:100)
      COMMON /vpcur/ velar(1001),vpval(1001),weight(1001),nvel
      gam  = MAX(1.0D-3,ABS(y(1)))
      Vgau = y(2)
      sig  = MAX(1.0D-3,ABS(y(3)))
      CALL GAUHERM (velar,vpval,weight,nvel,gam,Vgau,sig,harr,2)
      CHI2H = 1.0D0 + ((harr(0)-1.0D0)**2.0D0) +
     &        (harr(1)**2.0D0) + (harr(2)**2.0D0)
      END


      SUBROUTINE GAUHERM (velar,vpval,weight,nvel,
     &                    gam,Vgau,sig,harr,nhord)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (pi=3.14159265358979D0)
      DIMENSION velar(1001),vpval(1001),weight(1001),harr(0:100)
      DO l=0,nhord
        harr(l) = 0.0D0
      END DO
      DO i=1,nvel
        w     = (velar(i)-Vgau)/sig
        DO l=0,nhord
          harr(l) = harr(l) +
     &      (vpval(i)*weight(i)*SDGAUSS(w)*H_POL(l,w))
        END DO
      END DO
      DO l=0,nhord
        harr(l) = harr(l) * 2.0D0 * SQRT(pi) / gam
      END DO
      END


      REAL*8 FUNCTION SDGAUSS (x)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (pi=3.14159265358979D0)
      SDGAUSS = (1.0D0/SQRT(2.0D0*pi)) * EXPP(-0.5D0*x*x)
      END


      REAL*8 FUNCTION HE_POL (l,x)
      IMPLICIT REAL*8 (a-h,o-z)
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
      END


      REAL*8 FUNCTION H_POL (l,x)
      IMPLICIT REAL*8 (a-h,o-z)
      H_POL = HE_POL(l,x*SQRT(2.0D0))
      END


      REAL*8 FUNCTION EXPP(x)
      IMPLICIT REAL*8 (a-h,o-z)
      IF (x.GE.-500D0) THEN
        EXPP = EXP(x)
      ELSE
        EXPP = 0.0D0
      END IF
      END


      REAL*8 FUNCTION gammaln(x)
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


      REAL*8 FUNCTION gammln(x)
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
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL*8 ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=5000)
      EXTERNAL funk
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
      if (iter.ge.ITMAX) return
      iter=iter+1
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0D0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0D0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5D0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo) then
              do 15 j=1,ndim
                psum(j)=0.5D0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          goto 1
        endif
      endif
      goto 2
      END


      REAL*8 FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL*8 fac,p(mp,np),psum(ndim),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
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
