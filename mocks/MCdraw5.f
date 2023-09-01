      PROGRAM MCDRAW
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Create a Monte-Carlo realization of the DFs in de Bruijne 
C et al. (1996), and write 6D cartesian particle phase space coordinates to file.
C      
C All routines are written in double precision. Note that this requires
C the statement IMPLICIT REAL*8 (a-h,o-z) at the beginning of each
C program part. Don't mix single and real precision ! This generally leads
C to disaster.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
CCCCCCCCCCCCCCCCCCCC
C
      PARAMETER (pi=3.14159265358979D0)
C
CCCCCCCCCCCCCCCCCCCC
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /DFcase/ icase
C
C icase=1 for the case I DFs, icase=2 for the case II DFs.
C
      COMMON /param/ gamma, beta, q, alpha, eta
      CHARACTER(500) :: fileplace
C
CCCCCCCCCCCCCCCCCCCC
C
C Set model parameters
C
CCCCCCCCCCCCCCCCCCCC      
C
      WRITE (*,*) 'Please answer the following questions about'
      WRITE (*,*) 'the parameters of the model:'
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Kepler (1) or Logarithmic (2) Potential ?'
      READ (*,*) ipot
      WRITE (*,*) ' '
C
      delta = 2.0D0 - DBLE(ipot)
C      
      WRITE (*,*) 'Power-law slope gamma of the mass density'
      READ (*,*) gamma
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Intrinsic axial ratio q of the mass density'
      READ (*,*) q
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Case I (1) or Case II (2) DF ?'
      READ (*,*) icase
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Anisotropy parameter beta of the DF'
      READ (*,*) beta
      WRITE (*,*) ' '
C
C The parameters s and t of the odd part of the DF are called
C eta and alpha in the coding convention of this program.
C
      WRITE (*,*) 'Odd part parameters s and t for the DF'
      READ (*,*) eta, alpha
      WRITE (*,*) ' '
C
CCCCCCCCCCCCCCCCCCCC
C      
C Fill the arrays with the coefficients of the power series that yield
C the intrinsic velocity moments.
C
      WRITE (*,*) 'Calculating all series coefficients'
      maxord = 20
      CALL FILLARRAYS(maxord)
      WRITE (*,*) ' '
C
CCCCCCCCCCCCCCCCCCCC
C
C Set parameters for Monte-Carlo drawing      
C
CCCCCCCCCCCCCCCCCCCC
C
      WRITE (*,*) 'Give # Monte-Carlo tracers per dataset'
      READ (*,*) Nsamp
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Give an integer to initialize the random sequence'
      READ (*,*) IDUM
      WRITE (*,*) 
C
C Make sure that the integer is negative as required for initialization
C      
      IDUM = (-1*ABS(IDUM)) - 1
C
      WRITE (*,*) 'Give inner and outer radius for Monte-Carlo drawings'
      WRITE (*,*) 'in dimensionless units'
      READ (*,*) ri, ro
      WRITE (*,*) ' '
C
C Calculate derived factors.
C      
      gi = ri**(3.0D0-gamma)
      go = ro**(3.0D0-gamma)
C
CCCCCCCCCCCCCCCCCCCC
C
C Asks the user to specify a path for the input file
C     
      WRITE (*,*) 'Please specify a path where to save the output file.'
      READ (*,'(A)') fileplace
      WRITE (*,*) ' '
C
C
CCCCCCCCCCCCCCCCCCCC 
C
C Open output file
C      
      OPEN (UNIT=11,FILE=fileplace,STATUS='UNKNOWN')
C 
CCCCCCCCCCCCCCCCCCCC
C
C Loop over Monte-Carlo points
C
CCCCCCCCCCCCCCCCCCCC      
C
C Initialize second velocity moments
C      
      vr2sum = 0.0D0
      vt2sum = 0.0D0
      v2sum  = 0.0D0
C      
      DO i=1,Nsamp
C
CCCCCCCCCCCCCCCCCCCC
C           
C Draw random radii.
C
C The number density scales as r^-gamma. That means that dM/dr scales as       
C r^{2-gamma}. The integral of that scales as r^{3-gamma}. We evaluate the 
C integral between ri and ro. So the intregand 
C   xi = (r^{3-gamma}-ri^{3-gamma}) / (ro^{3-gamma}-ri^{3-gamma})
C ranges from 0 to 1 . So we pick a random xi, and then calculate the 
C corresponding radius r by inverting the equation. 
C
C Note gamma=3 provides an exception. Then the integrand scales as log(r).
C So the intregand                                                  
C   xi = (log(r/ri)) / (log(ro/ri))
C ranges from 0 to 1 . So we pick a random xi, and then calculate the
C corresponding radius r by inverting the equation.             
C        
CCCCCCCCCCCCCCCCCCCC
C      
        xi = RAN1(IDUM)
        IF (gamma.EQ.3.0D0) THEN
          rad = ri * ((ro/ri)**xi)
        ELSE    
          rad = (gi + (xi*(go-gi)))**(1.0D0/(3.0D0-gamma))
        END IF         
C
CCCCCCCCCCCCCCCCCCCC
C
C Draw theta as described in draw_particles.txt.
C When done, we add a random sign.      
C
CCCCCCCCCCCCCCCCCCCC
C
 201    st    = RAN1(IDUM)
        th    = ASIN(st)
        ct    = COS(th)
        fmass = (((ct**2)+((st/q)**2))**(-0.5D0*gamma))
        xi    = RAN1(IDUM) 
        IF (xi.GT.fmass) GOTO 201
        xxi   = RAN1(IDUM)
        IF (xxi.LE.0.5D0) THEN
          th = -1.0D0*th
        END IF   
C
CCCCCCCCCCCCCCCCCCCC
C        
C Set up 3D positions for all particles. Calculate the value thcur that is
C (as in the paper) calculated from the symmerty axis, instead of from        
C the symmetry plane (as is th).
C
CCCCCCCCCCCCCCCCCCCC
C        
        rcur    = rad
        thcur   = (pi/2.0D0) - th
        phicur  = 2.0D0*pi*RAN1(IDUM)
C     
CCCCCCCCCCCCCCCCCCCC
C
C Calculate the velocity moments using the numerical machinery of this program.
C
C Results of RHOVELMOM are at radius of 1 in dimensionless units.
C They need to be scaled to the current radius using the scale-free
C nature of the models.
C 
        rho     = RHOVELMOM(thcur,0,0,0)
C        
        vradsec = (1.0D0/(rcur**delta))*RHOVELMOM(thcur,2,0,0)/rho
        vthsec  = (1.0D0/(rcur**delta))*RHOVELMOM(thcur,0,2,0)/rho
        vphsec  = (1.0D0/(rcur**delta))*RHOVELMOM(thcur,0,0,2)/rho
C
CCCCCCCCCC
C        
C Draw Gaussian deviates from the intrinsic velocity dispersions
C        
        vrcur  = GASDEV(IDUM) * SQRT(vradsec)
        vthcur = GASDEV(IDUM) * SQRT(vthsec)
        vphcur = GASDEV(IDUM) * SQRT(vphsec)
C
C Write phase-space coordinates of Monte-Carlo points
C
        cth   = COS(thcur)
        sth   = SIN(thcur)
        cph   = COS(phicur)
        sph   = SIN(phicur)
C       
        Rcyl  = rcur*sth
        zz    = rcur*cth
C        
        vRcyl = (vrcur*sth) + (vthcur*cth)        
        vz    = (vrcur*cth) - (vthcur*sth)
C
        xx    = Rcyl*cph
        yy    = Rcyl*sph
C        
        vx    = (vRcyl*cph) - (vphcur*sph)
        vy    = (vRcyl*sph) + (vphcur*cph)
C        
        WRITE (11,'(6E15.5)') xx, yy, zz, vx, vy, vz
C
C Add results to runnings sums needed to calculated the Binney beta
C
        vr2sum  = vr2sum + (vrcur**2)
        vt2sum  = vt2sum + (vthcur**2) + (vphcur**2)
        v2sum   = v2sum  + ((vx**2)+(vy**2)+(vz**2))
C
C End the i-loop over data points
C        
      END DO
C
C Calculate the Binney beta (as opposed to the DF parameter beta)
C
      betabin = 1.0D0 - (vt2sum/(2.0D0*vr2sum))
C      
C Write statistics of Monte-Carlo results
C
      WRITE (*,*) 'Binney beta,        v_RMS(calculated 2 ways):'
      WRITE (*,'(F10.6,A10,2F10.6)') betabin, ' ',
     &   SQRT(((vr2sum+vt2sum)/DBLE(Nsamp))), 
     &   SQRT((v2sum/DBLE(Nsamp)))
C      
CCCCCCCCCCCCCCCCCCCC
C
C Close output file
C
      CLOSE (UNIT=11)
C      
      END


      SUBROUTINE FILLARRAYS (maxxord)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Fill the arrays with the coefficients of the power series that yield
C the intrinsic velocity moments. It is assumed that moments up to order 
C maxxord are required. After this subroutine has been called, 
C the subroutines COFARR and KMAXARR can be used to recover the coefficients, 
C without having to calculate them again.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (maxlength = 2000000 ,
     &           mord      = 200     ,
     &           mev       = mord/2  )
C
      COMMON /all/ allcofs(1:maxlength),
     &             icstart(0:mev,0:mev,0:mord),
     &             kmaxall(0:mev,0:mev,0:mord), maxord
C
      DIMENSION cofs(0:1000)
C
CCCCCCCCCCCCCCCCCCCC
C
      maxord = maxxord
C
      IF (2*(maxord/2).NE.maxord) STOP 'maxord should be even'
      IF (maxord.GT.mord) STOP 'maxord is too large for array'
C
C Fill the array which holds all the coefficients of the power series
C in e^2 SIN^2(theta), required to calculate the intrinsic velocity
C moments up to order maxord.
C
C      WRITE (*,*) 
C     & 'Calculating series coefficients for the velocity moments ...'
C
      ic = 1
      DO ir=0,maxord,2
        DO ith=0,maxord,2
          DO iph=0,maxord
            icstart(ir/2,ith/2,iph) = ic
            CALL CALCCOEFFLN(cofs,kmax,ir,ith,iph)
            kmaxall(ir/2,ith/2,iph) = kmax
            DO k=0,kmax
              allcofs(ic) = cofs(k)
              ic = ic+1
            END DO
            IF (ic.GE.maxlength) STOP 'ic too large for array'
          END DO
        END DO
C        WRITE (*,'(2I5,I10)') ir,maxord,ic-1
      END DO
C      WRITE (*,*) ' '
C
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The function RHOVELMOM calculates rho times intrinsic velocity momemnt
C of a given order. This is a prefactor times a power series in
C e^2 SIN^2(theta). The (ln of the) prefactor is calculated by PREFACLN.
C The subroutines COFARR and KMAXARR allow the power series 
C to be reconstructed, by reading values from the arrays allcofs and 
C kmaxall, in the common block /all/. These arrays must have been filled
C previously in the main program, using the subroutine CALCCOEFFLN.
C This subroutine calls the function COEFFLN for each term in the series.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION RHOVELMOM(theta,lr,lth,lph)
C
C Calculates the rho times velocity moment of order (lr,lth,lph).
C The fractional accuracy is eps. The polar angle theta must be given 
C in radians. The dimensionless radius is assumed to be unity.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (eps = 1.0D-10)
C
      COMMON /param/ gamma, beta, q, alpha, eta
C
      DIMENSION cofs(0:1000)
C
      IF ((2*(lr/2).NE.lr).OR.(2*(lth/2).NE.lth)) THEN
        RHOVELMOM = 0.0D0
        RETURN
      END IF
C
      kmax = KMAXARR(lr,lth,lph)
      DO k=0,kmax
        cofs(k) = COFARR(k,lr,lth,lph)
      END DO
C
      prefln = PREFACLN(lr,lth,lph)
C
      e2 = 1.0D0 - (q*q)
      s2 = (SIN(theta))**2.0D0
C
      k       = 0
      xnewadd = EXPP(cofs(0))
      sum     = xnewadd
      epsab   = eps*EXPP(cofs(0))
C
      IF ((e2*s2).GE.1.0D-12) THEN
C
35      k         = k+1
        oldadd    = xnewadd
        xnewadd   = EXPP(cofs(k)+(DBLE(k)*LOG(e2*s2)))
        sum       = sum + xnewadd
        IF ( (.NOT.((ABS(xnewadd).LE.epsab).AND.
     &              (ABS(oldadd).LE.epsab))) .AND. 
     &       (k.LE.kmax) ) THEN
          GOTO 35
        END IF
C
      END IF
C
      IF (2*(lph/2).EQ.lph) THEN
        RHOVELMOM = EXPP(prefln) * sum
      ELSE
        IF (s2.LE.1.0D-12) THEN
          RHOVELMOM = 0.0D0
        ELSE
          RHOVELMOM = EXPP(prefln) * sum * ((2.0D0*eta)-1.0D0) *
     &         (s2**alpha)
        END IF
      END IF
C
      END


      REAL*8 FUNCTION PREFACLN(lr,lth,lph)
C
C Calculates the logarithm of the pre-factor by which the relevant 
C series must be multiplied to get the velocity moment
C of order (lr,lth,lph). This function does not take into account the factor
C (2 eta - 1) * (sin^2(theta))^{alpha} that must be added for the odd
C moments.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /param/ gamma, beta, q, alpha, eta
C
      IF ((2*(lr/2).NE.lr).OR.(2*(lth/2).NE.lth)) THEN
        PREFACLN = -100.0D0
        RETURN
      END IF
C
      IF (2*(lph/2).EQ.lph) THEN
        alp = 0.0D0
      ELSE
        alp = alpha
      END IF
C
      dlr  = DBLE(lr)
      dlth = DBLE(lth)
      dlph = DBLE(lph)
C
      PREFACLN = (gamma*LOG(q)) + 
     &  betaln(1.0D0+alp-beta+(0.5D0*(dlth+dlph)),
     &                (0.5D0*(dlr+1.0D0))) +
     &  betaln(alp+(0.5D0*(dlph+1.0D0)),(0.5D0*(dlth+1.0D0))) -
     &  betaln(1.0D0-beta,0.5D0) - betaln(0.5D0,0.5D0)
C
      IF (ipot.EQ.1) THEN
        PREFACLN = PREFACLN + (alp*LOG(4.0D0)) +
     &    betaln(1.5D0-beta+alp+(0.5D0*(dlr+dlth+dlph)),
     &           gamma-beta-0.5D0+alp) -
     &    betaln(1.5D0-beta,gamma-beta-0.5D0)
      ELSE IF (ipot.EQ.2) THEN
        PREFACLN = PREFACLN + (alp*(1.0D0+LOG(2.0D0))) +
     &    gammaln(1.5D0-beta+alp+(0.5D0*(dlr+dlth+dlph))) -
     &    gammaln(1.5D0-beta)
      ELSE
        STOP 'ipot wrong value'
      END IF
C
      END


      REAL*8 FUNCTION COFARR(k,lr,lth,lph)
C
C Get a value from the array allcofs
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (maxlength = 2000000 ,
     &           mord      = 200     ,
     &           mev       = mord/2  )
C
      COMMON /all/ allcofs(1:maxlength),
     &             icstart(0:mev,0:mev,0:mord),
     &             kmaxall(0:mev,0:mev,0:mord),maxord
C
      IF ((lr.GT.maxord).OR.(lth.GT.maxord).OR.
     &    (lph.GT.maxord)) THEN
        STOP 'order to large in cofarr'
      END IF
C
      IF ((2*(lr/2).NE.lr).OR.(2*(lth/2).NE.lth)) THEN
        COFARR = -100.0D0
      ELSE
        COFARR = allcofs(k+icstart(lr/2,lth/2,lph))      
      END IF
C
      END


      INTEGER FUNCTION KMAXARR(lr,lth,lph)
C
C Get a value from the array kmaxall
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (maxlength = 2000000 ,
     &           mord      = 200     ,
     &           mev       = mord/2  )
C
      COMMON /all/ allcofs(1:maxlength),
     &             icstart(0:mev,0:mev,0:mord),
     &             kmaxall(0:mev,0:mev,0:mord),maxord
C
      IF ((lr.GT.maxord).OR.(lth.GT.maxord).OR.
     &    (lph.GT.maxord)) THEN
        STOP 'order too large in kmaxarr'
      END IF
C
      IF ((2*(lr/2).NE.lr).OR.(2*(lth/2).NE.lth)) THEN
        KMAXARR = 0
      ELSE
        KMAXARR = kmaxall(lr/2,lth/2,lph)
      END IF
C
      END


      SUBROUTINE CALCCOEFFLN(cofs,kmax,lr,lth,lph)
C
C Calculates the logarithm of coefficients in the hypergeometric series
C that occurs in the expression for the velocity moment
C of order (lr,lth,lph). The array cofs is filled up to order kmax,
C which is determined by this subroutine so as to allow the
C power series to be evaluated to a fractional accuracy eps.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (eps = 1.0D-10)
C
      DIMENSION cofs(0:1000)
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /order/ llr,llth,llph
C
C Fill the common block order
C
      llr  = lr
      llth = lth
      llph = lph
C
      e2 = 1.0D0 - (q*q)
C
      k       = 0
      cofs(0) = COEFFLN(0)
      xnewadd = EXPP(cofs(0))
      sum     = xnewadd
      epsab   = eps*EXPP(cofs(0))
C
      IF (e2.GE.1.0D-20) THEN
15      k         = k+1
        oldadd    = xnewadd
        cofs(k)   = COEFFLN(k)
        xnewadd   = EXPP(cofs(k)+(DBLE(k)*LOG(e2)))
        sum       = sum + xnewadd
        IF ( (.NOT.((ABS(xnewadd).LE.epsab).AND.
     &            (ABS(oldadd).LE.epsab))) .AND. 
     &     (k.LE.1000) ) THEN
          GOTO 15
        ELSE
          kmax = k
        END IF
      ELSE
        kmax = 0
      END IF
C
      END


      REAL*8 FUNCTION COEFFLN(k)
C
C The logarithm of the k-th coefficient in the hypergeometric series
C that occurs in the expression for the velocity moment
C of order (lr,lth,lph) (as contained in the common block /order/)
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /DFcase/ icase
C
C icase=1 for the case I DFs, icase=2 for the case II DFs.
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /order/ lr,lth,lph
C
      IF ((2*(lr/2).NE.lr).OR.(2*(lth/2).NE.lth)) THEN
        COEFFLN = 0.0D0
        RETURN
      END IF
C
      IF (2*(lph/2).EQ.lph) THEN
        alp = 0.0D0
      ELSE
        alp = alpha
      END IF
C
      dk   = DBLE(k)
      dlr  = DBLE(lr)
      dlth = DBLE(lth)
      dlph = DBLE(lph)
C
      COEFFLN = pochln(0.5D0*gamma,k) - gammaln(dk+1.0D0) +
     &  pochln(alp+(0.5D0*(dlph+1.0D0)),k) -
     &  pochln(0.5D0,k) +
     &  pochln(1.0D0,k) -
     &  pochln(1.0D0+alp+(0.5D0*(dlth+dlph)),k)
C
      IF (icase.EQ.1) THEN
C
        IF (ipot.EQ.1) THEN
          COEFFLN = COEFFLN +
     &      pochln(1.0D0-beta+alp+(0.5D0*(dlth+dlph)),k) -
     &      pochln(1.0D0-beta,k) +
     &      pochln(gamma-0.5D0-beta+alp,k) -
     &      pochln(gamma-0.5D0-beta,k) +
     &      pochln((0.5D0*(gamma+1.0D0))-beta,k) -
     &      pochln((0.5D0*(gamma+1.0D0))-beta+alp+
     &             (0.25D0*(dlr+dlth+dlph)),k) +
     &      pochln(1.0D0-beta+(0.5D0*gamma),k) -
     &      pochln(1.0D0-beta+(0.5D0*gamma)+alp+
     &             (0.25D0*(dlr+dlth+dlph)),k) 
        ELSE IF (ipot.EQ.2) THEN
          COEFFLN = COEFFLN +
     &      pochln(1.0D0-beta+alp+(0.5D0*(dlth+dlph)),k) -
     &      pochln(1.0D0-beta,k) +
     &      ( (beta-1.5D0-dk-alp-(0.5D0*(dlr+dlth+dlph)))*
     &        LOG(gamma-(2.0D0*beta)+(2.0D0*dk)+(2.0D0*alp)) ) -
     &      ( (beta-1.5D0-dk)*
     &        LOG(gamma-(2.0D0*beta)+(2.0D0*dk)) )
        ELSE
          STOP 'ipot wrong value'
        END IF
C
      ELSE IF (icase.EQ.2) THEN
C
        IF (ipot.EQ.1) THEN
C
C In this case the prefactor is the entire value of coeffln
C
        ELSE IF (ipot.EQ.2) THEN
C
C Note that this factor could also be put in the prefactor PREFACLN because
C it is independent of k (or dk=DOUBLE(k))
C
          COEFFLN = COEFFLN +
     &    ( (beta-1.5D0-alp-(0.5D0*(dlr+dlth+dlph)))*
     &      LOG(gamma-(2.0D0*beta)+(2.0D0*alp)) ) -
     &    ( (beta-1.5D0)*
     &      LOG(gamma-(2.0D0*beta)) )
C
        ELSE
          STOP 'ipot wrong value'
        END IF
C
      ELSE
C
        STOP 'icase wrong value'
C
      END IF
C
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


      REAL*8 FUNCTION betaln(x,y)
C
C The logarithm of the beta function
C
      IMPLICIT REAL*8 (a-h,o-z)
      betaln = gammaln(x) + gammaln(y) - gammaln(x+y)
      END

  
      REAL*8 FUNCTION pochln(x,k)
C
C The logarithm of Pochhammer's symbol
C
      IMPLICIT REAL*8 (a-h,o-z)
      pochln = gammaln(x+DBLE(k)) - gammaln(x)
      END


      REAL*8 FUNCTION binomln(j,k)
C
C The logarithm of the binomial coefficient (j over k)
C
      IMPLICIT REAL*8 (a-h,o-z)
      dj = DBLE(j)
      dk = DBLE(k)
      binomln = gammaln(dj+1.0D0) - gammaln(dk+1.0D0) -
     &          gammaln(dj-dk+1.0D0)
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
C Required routines from Numerical recipes. All were modified to
C be double precision. The parameter eps in QROMO determines the
C speed of the program in calculating the VP parameters. Best results
C are obtained with eps=1.0D-8, but this leads to a slow program.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      REAL*8 FUNCTION ran1(idum)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.0D0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2D-7,RNMX=1.0D0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
      

      REAL*8 FUNCTION gasdev(idum)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER idum
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.0D0*ran1(idum)-1.0D0
        v2=2.0D0*ran1(idum)-1.0D0
        rsq=v1**2+v2**2
        if(rsq.ge.1.0D0.or.rsq.eq.0.0D0)goto 1
        fac=sqrt(-2.0D0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

      
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
