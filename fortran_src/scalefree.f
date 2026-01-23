      PROGRAM FLATPOWER
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Scale-free dynamical modeling software as discussed in the file README.
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
      DIMENSION xmom(0:100),darr(0:100),harr(0:100),
     &          velar(101),vpval(101)
      CHARACTER*256 outpath
      INTEGER iout
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
C
C Model parameters
C 
      COMMON /viewing/ xinc
C
C The inclination angle in radians 
C
      COMMON /howint/ iint
C
C Determines how the line of sight integration is done
C
      COMMON /romoeps/ eps
C
C fractional accuracy for Romberg integration
C 
      COMMON /gleg/ qx(300),qw(300),nGL
C
C Common block with Gauss-Legendre coefficients for integration along the
C line of sight.
C
      COMMON /reg/ xlam
C
C Common block with the regularization parameter used in calculating VPs
C
      COMMON /smoo/ epsmoo
C
C Parameter that determines the requested smoothness, as fraction of
C the VP maximum.
C
      COMMON /verbose/ iverb
C
      COMMON /projcase/ iproj
C
C iproj=1 (LOS), 2 (POSR), 3 (POST)
C
C
C Determines whether verbose output should be generated for VP calculations
C
CCCCCCCCCCCCCCCCCCCC
C
C Get the model parameters
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
      WRITE (*,*) 'Viewing inclination i in degrees (90=edge-on)'
      READ (*,*) xinc
      WRITE (*,*) ' '
C
C Transform the inclination to radians
C
      xinc = (xinc/180.0D0)*pi
C
CCCCCCCCCCCCCCCCCCCC
C
C Get information on numerical details
C
CCCCCCCCCCCCCCCCCCCC
C
      WRITE (*,*) 'Please answer the following questions about'
      WRITE (*,*) 'the numerical details of the calculations:'
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Use Romberg (0) or Gauss-Legendre (1) integration'
      WRITE (*,*) 'for line-of-sight projection and integration over'
      WRITE (*,*) 'the meridional or sky plane? (default=1)' 
      READ (*,*) iint
      WRITE (*,*) ' '
C
      IF (iint.EQ.0) THEN
        WRITE (*,*) 'Give the fractional accuracy epsilon (0.0=default)'
        READ (*,*) eps
        WRITE (*,*) ' '
        IF (eps.EQ.0.0D0) eps = 1.0D-7
      ELSE
        WRITE (*,*) 'Give number of quadrature points (0=default)'
        READ (*,*) nGL
        WRITE (*,*) ' '
        IF (nGL.EQ.0) nGL = 100
        CALL GAULEG (0.0D0,1.0D0,qx,qw,nGL)
      END IF    
C
      WRITE (*,*) 'Algorithm to calculate VPs and Gauss-Hermite moments'
      WRITE (*,*) '   1: Solve VanderMonde matrix directly without'
      WRITE (*,*) '      regularization. Resulting VP will be nonsense'
      WRITE (*,*) '      but the GH moments are generally well'
      WRITE (*,*) '      determined.'
      WRITE (*,*) '   2: Use regularization with a fixed regularization'
      WRITE (*,*) '      parameter.'
      WRITE (*,*) '   3: Use regularization. Increase regularization'
      WRITE (*,*) '      parameter until the VP has no more than 3'
      WRITE (*,*) '      significant local maxima. A local maximum is'
      WRITE (*,*) '      significant if it exceeds the value of its'
      WRITE (*,*) '      neighbors on the grid by eps times the'
      WRITE (*,*) '      absolute VP maximum.'
      WRITE (*,*) 'Choose 1 for default.'
      READ (*,*) iVPfit
      WRITE (*,*) ' '
C
      IF (iVPfit.EQ.1) THEN
C
        xlam = -1.0
C
        WRITE (*,*) 'Give the maximum number of projected moments'
        WRITE (*,*) 'to use (should be an even number)' 
        WRITE (*,*) '(0 yields default)'
        READ (*,*) maxord
        WRITE (*,*) ' '
C
      ELSE IF (iVPfit.EQ.2) THEN
C
        WRITE (*,*) 'Give the regularization parameter (> 0)'
        READ (*,*) xlam
        WRITE (*,*) ' '
C
        WRITE (*,*) 'Give the number of projected moments'
        WRITE (*,*) 'to use (should be an even number)'
        WRITE (*,*) '(0 yields default)'
        READ (*,*) maxord
        WRITE (*,*) ' '
C
      ELSE IF (iVPfit.EQ.3) THEN
C
        xlam = 0.0      
C
        WRITE (*,*) 'Give smoothness factor eps (0.0 yields default)'
        READ (*,*) epsmoo
        WRITE (*,*) ' ' 
C
        IF (epsmoo.EQ.0.0) epsmoo = 1.0D-3

        WRITE (*,*) 'Give the number of projected moments'
        WRITE (*,*) 'to use (should be an even number)'
        WRITE (*,*) '(0 yields default)'
        READ (*,*) maxord
        WRITE (*,*) ' '
C
      ELSE 
C
        STOP 'Wrong answer'
C
      END IF
C
      IF (maxord.EQ.0) maxord = 30
C
CCCCCCCCCCCCCCCCCCCC
C
C Fill the arrays with the coefficients of the power series that yield
C the intrinsic velocity moments.
C
CCCCCCCCCCCCCCCCCCCC
C
      WRITE (*,*) 'Calculating all series coefficients up to'
      WRITE (*,*) 'the selected (maximum) order ........'
      CALL FILLARRAYS(maxord)
      WRITE (*,*) ' '
C
CCCCCCCCCCCCCCCCCCCC
C
C Now calculate anything that the user may be interested in
C
CCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCC
C
C Select where to write numerical results (table output)
C
      iproj = 1
      WRITE (*,*) 'Output file for results (blank = STDOUT)'
      READ (*,'(A)') outpath
      lout = LENSTR(outpath)
      IF (lout.GT.0) THEN
        iout = 20
        OPEN(unit=iout,file=outpath(1:lout),status='unknown')
      ELSE
        iout = 6
      END IF
      WRITE (*,*) ' '
C
51    WRITE (*,*) 'Calculate intrinsic (0) or projected (1)'
      WRITE (*,*) 'kinematical quantities ?'
      WRITE (*,*) 'This gives results for a fixed angle in the'
      WRITE (*,*) 'meridional or projected plane.'
      WRITE (*,*) 'Instead, add 2 to get results mass-weighted' 
      WRITE (*,*) 'over angles between 0 and pi/2'
      READ (*,*) iwhat
      WRITE (*,*) ' '
C
      WRITE (*,*) 'Note: all results are at an'
      WRITE (*,*) '(intrinsic or projected) radius of 1 in'
      WRITE (*,*) 'dimensionless units. Results can be scaled to other'
      WRITE (*,*) 'radii using the scale-free nature of the models.'
      WRITE (*,*) ' '
C
      IF (iwhat.EQ.0) THEN
C
        WRITE (*,*) 'Give angle theta in the meridional plane'
        WRITE (*,*) '(in degrees) (0 = symmetry axis)'
        READ (*,*) theta
        WRITE (*,*) ' '
C
C Transform to radians 
C
        theta = (theta/180.0D0)*pi
        rho = RHOVELMOM(theta,0,0,0)
C
        WRITE (*,*) 'Intrinsic velocity moments:'
        WRITE (*,'(5A12)') 'rho','<v_ph>','<v_r^2>','<v_th^2>',
     &                     '<v_ph^2>'   
        WRITE (*,'(6F12.8)') rho,RHOVELMOM(theta,0,0,1)/rho,
     &                           RHOVELMOM(theta,2,0,0)/rho,
     &                           RHOVELMOM(theta,0,2,0)/rho,
     &                           RHOVELMOM(theta,0,0,2)/rho
        WRITE (*,*) ' '
C
      ELSE IF (iwhat.EQ.2) THEN
C
        rho   = RHOVELMOMTHAV(0,0,0)
        rhop1 = RHOVELMOMTHAV(0,0,1)
        rhor2 = RHOVELMOMTHAV(2,0,0)
        rhot2 = RHOVELMOMTHAV(0,2,0)
        rhop2 = RHOVELMOMTHAV(0,0,2)
        betav = 1.0D0 - ((rhot2+rhop2)/(2.0D0*rhor2))
C
        WRITE (*,*) 'Intrinsic velocity moments:'
        WRITE (*,*) 'Mass-weighted average spherical shell:'
        WRITE (*,'(6A12)') 'rho','<v_ph>','<v_r^2>','<v_th^2>',
     &                     '<v_ph^2>','beta'   
        WRITE (*,'(6F12.8)') rho, rhop1/rho,
     &        rhor2/rho, rhot2/rho, rhop2/rho, betav
        WRITE (*,*) ' '
C
      ELSE IF (iwhat.EQ.1) THEN
C
        WRITE (*,*) 'Give angle on the projected plane'
        WRITE (*,*) '(in degrees) (0 = major axis)'
        READ (*,*) xi
        WRITE (*,*) ' '
C
C Transform to radians
C
        xi = (xi/180.0D0)*pi
C
        rhop = RHOPROJ(xi)
        WRITE (iout,'(A)') '# kind=projected_point'
        WRITE (iout,'(A,F16.8)') '# xi_deg ', (xi*180.0D0/pi)
        WRITE (iout,'(A)') '# columns: iproj rho_p v1 v2 v3 v4'
        DO iproj=1,3
          v1 = PROJMOM(xi,1)
          v2 = PROJMOM(xi,2)
          v3 = PROJMOM(xi,3)
          v4 = PROJMOM(xi,4)
          WRITE (iout,'(I3,1X,5E24.16)') iproj, rhop, v1, v2, v3, v4
        END DO
        WRITE (iout,*) ' '
C
      ELSE IF (iwhat.EQ.3) THEN
C
        rhop0 = RHOPROJMOMAV(0)
        rhop1 = RHOPROJMOMAV(1)
        rhop2 = RHOPROJMOMAV(2)
        rhop3 = RHOPROJMOMAV(3)
        rhop4 = RHOPROJMOMAV(4)
C         
        WRITE (iout,'(A)') '# kind=projected_circle_average'
        WRITE (iout,'(A)') '# columns: iproj rho_p v1 v2 v3 v4'
        DO iproj=1,3
          rhop0 = RHOPROJMOMAV(0)
          rhop1 = RHOPROJMOMAV(1)
          rhop2 = RHOPROJMOMAV(2)
          rhop3 = RHOPROJMOMAV(3)
          rhop4 = RHOPROJMOMAV(4)
          WRITE (iout,'(I3,1X,5E24.16)') iproj, rhop0, rhop1/rhop0,
     &                              rhop2/rhop0, rhop3/rhop0, rhop4/rhop0
        END DO
        WRITE (iout,*) ' '
C
      ELSE
C
        STOP 'Wrong answer'
C
      END IF
C
CCCCCCCCCCCCCCCCCCCC
C
C If iwhat=1 or 3, then calculate and write also VP information      
C
CCCCCCCCCCCCCCCCCCCC
C
      IF ((iwhat.EQ.1).OR.(iwhat.EQ.3)) THEN
C
C Set the number of GH moments to calculate
C
        nord = 6
C
C Verbose ?
C
        WRITE (*,*) 'Give verbose output of intermediate steps'
        WRITE (*,*) 'for VP calculation ? (0/1)'
        READ (*,*) iverb
        WRITE (*,*) ' '
C
C Calculate the VP and VP coefficients for all projected components
C
        WRITE (iout,'(A)') '# kind=vp'
        WRITE (iout,'(A)') '# columns: iproj true_gam true_V true_sig'
        WRITE (iout,'(A)') '#          gauss_gam gauss_V gauss_sig'
        WRITE (iout,'(A)') '#         h0 h1 h2 h3 h4 h5 h6'

        DO iproj=1,3
          IF (xlam.LT.0.0) THEN
            CALL VPANALYSE_CONV (iwhat,xi,nord,xmom,gam0,V0,sig0,darr,
     &                      velar,vpval,nvel,gam,Vgau,sig,harr)
          ELSE
            CALL VPANALYSE_FIX (iwhat,xi,maxord,xmom,gam0,V0,sig0,darr,
     &                      velar,vpval,nvel,gam,Vgau,sig,harr)
          END IF
C
          WRITE (iout,'(I3,1X,6E24.16,1X,7F10.5)') iproj, gam0, V0, sig0,
     &            gam, Vgau, sig, harr(0), harr(1), harr(2), harr(3),
     &            harr(4), harr(5), harr(6)
C
C Write the VP solution as a table block (v, VP)
C
          WRITE (iout,'(A,I3)') '# vp_table iproj ', iproj
          WRITE (iout,'(A)') '# columns: v vp'
          DO j=1,nvel
            WRITE (iout,'(2E24.16)') velar(j), vpval(j)
          END DO
          WRITE (iout,*) ' '
        END DO
C
C
C If iwhat=0 or 2, then calculate and write also intrinsic VP information
C (reconstructed from intrinsic moments of the selected velocity component)
C
CCCCCCCCCCCCCCCCCCCC
C
      ELSE IF ((iwhat.EQ.0).OR.(iwhat.EQ.2)) THEN
C
C Set the number of GH moments to calculate
C
        nord = 6
C
C Verbose ?
C
        WRITE (*,*) 'Give verbose output of intermediate steps'
        WRITE (*,*) 'for VP calculation ? (0/1)'
        READ (*,*) iverb
        WRITE (*,*) ' '
C
C For intrinsic VPs we interpret iproj as the intrinsic component:
C   1 = v_r, 2 = v_theta, 3 = v_phi
C
        IF (iwhat.EQ.0) THEN
C
C Structured intrinsic moments at (r=1, theta)
C
          rho_int  = RHOVELMOM(theta,0,0,0)
          vphi_int = RHOVELMOM(theta,0,0,1) / rho_int
          vr2_int  = RHOVELMOM(theta,2,0,0) / rho_int
          vth2_int = RHOVELMOM(theta,0,2,0) / rho_int
          vphi2_int= RHOVELMOM(theta,0,0,2) / rho_int
C
          WRITE (iout,'(A)') '# kind=intrinsic_point'
          WRITE (iout,'(A)') '# columns: rho vphi vr2 vth2 vphi2'
          WRITE (iout,'(5E24.16)') rho_int, vphi_int, vr2_int,
     &                             vth2_int, vphi2_int
        ELSE
C
C Structured intrinsic moments: mass-weighted spherical shell average
C
          rho_int  = RHOVELMOMTHAV(0,0,0)
          vphi_int = RHOVELMOMTHAV(0,0,1) / rho_int
          vr2_int  = RHOVELMOMTHAV(2,0,0) / rho_int
          vth2_int = RHOVELMOMTHAV(0,2,0) / rho_int
          vphi2_int= RHOVELMOMTHAV(0,0,2) / rho_int
          betav    = 1.0D0 - ((vth2_int+vphi2_int)/(2.0D0*vr2_int))
C
          WRITE (iout,'(A)') '# kind=intrinsic_shell_average'
          WRITE (iout,'(A)') '# columns: rho vphi vr2 vth2 vphi2 beta'
          WRITE (iout,'(6E24.16)') rho_int, vphi_int, vr2_int,
     &                             vth2_int, vphi2_int, betav
        END IF
C
C Calculate the VP and VP coefficients for all intrinsic components
C
        WRITE (iout,'(A)') '# kind=vp_intrinsic'
        WRITE (iout,'(A)') '# columns: icomp true_gam true_V true_sig'
        WRITE (iout,'(A)') '#          gauss_gam gauss_V gauss_sig'
        WRITE (iout,'(A)') '#         h0 h1 h2 h3 h4 h5 h6'
C
        DO iproj=1,3
          IF (xlam.LT.0.0) THEN
            IF (iwhat.EQ.0) THEN
              ang = theta
            ELSE
              ang = 0.0D0
            END IF
            CALL VPANALYSE_CONV (iwhat,ang,nord,xmom,gam0,V0,sig0,darr,
     &                      velar,vpval,nvel,gam,Vgau,sig,harr)
          ELSE
            IF (iwhat.EQ.0) THEN
              ang = theta
            ELSE
              ang = 0.0D0
            END IF
            CALL VPANALYSE_FIX (iwhat,ang,maxord,xmom,gam0,V0,sig0,darr,
     &                      velar,vpval,nvel,gam,Vgau,sig,harr)
          END IF
C
          WRITE (iout,'(I3,1X,6E24.16,1X,7F10.5)') iproj, gam0, V0, sig0,
     &            gam, Vgau, sig, harr(0), harr(1), harr(2), harr(3),
     &            harr(4), harr(5), harr(6)
C
C Write the VP solution as a table block (v, VP)
C
          WRITE (iout,'(A,I3)') '# vp_table icomp ', iproj
          WRITE (iout,'(A)') '# columns: v vp'
          DO j=1,nvel
            WRITE (iout,'(2E24.16)') velar(j), vpval(j)
          END DO
          WRITE (iout,*) ' '
        END DO
C
      END IF
C      
CCCCCCCCCCCCCCCCCCCC
C
C Continue ?
C
CCCCCCCCCCCCCCCCCCCC
C
      WRITE (*,*) 'Calculate something else for this model ? (0/1)'
      READ (*,*) imore
      WRITE (*,*) ' '
C
      IF (imore.EQ.1) GOTO 51
C
CCCCCCCCCCCCCCCCCCCC
C
C End of program
C
CCCCCCCCCCCCCCCCCCCC
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
C The function SPHPROJMOM calculatates the results for spherical non-rotating
C models using the results of van der Marel & Franx, for the purpose of testing
C only.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION SPHPROJMOM(nord)
C
C Calculate the projected moment of order nord at projected radius 1 
C (in dimensionless units), using equation (B6a) or (B8a) of van der Marel 
C & Franx (modified to correspond to the units used here).
C Useful for the purpose of testing.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /param/ gamma, beta, q, alpha, eta
C
      IF (2*(nord/2).NE.nord) STOP 'No vdM & F result for odd order'
C
      no2   = nord/2
      dnord = DBLE(nord)
      dno2  = DBLE(no2)
C
      sum = 0.0D0
      DO k=0,no2
        dk = DBLE(k)
        facln = binomln(nord,2*k) + 
     &    betaln(dno2-dk+0.5D0,0.5D0) - betaln(0.5D0,0.5D0) +
     &    betaln(dk+0.5D0,dno2-dk+1.0D0-beta) -
     &    betaln(0.5D0,1.0D0-beta) +
     &    gammaln(dno2+1.5D0-beta) - gammaln(1.5D0-beta)
        IF (ipot.EQ.1) THEN
          facln = facln + 
     &       betaln(0.5D0*(gamma-1.0D0+dno2+dnord-(2.0D0*dk)),
     &           dk+0.5D0) - betaln(0.5D0*(gamma-1.0D0),0.5D0) +
     &       gammaln(gamma-(2.0D0*beta)+1.0D0) -
     &       gammaln(gamma+dno2-(2.0D0*beta)+1.0D0)
        ELSE IF (ipot.EQ.2) THEN
          facln = facln +
     &       betaln(0.5D0*(gamma-1.0D0+dnord-(2.0D0*dk)),
     &           dk+0.5D0) - betaln(0.5D0*(gamma-1.0D0),0.5D0) -
     &       (dno2*LOG(gamma-(2.0D0*beta)))
        ELSE
          STOP 'ipot wrong value'
        END IF
        sum = sum + EXPP(facln)
      END DO
C
      SPHPROJMOM = sum      
C
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The subroutines VPANALYSE_* calculate properties of the VPs on the sky.
C They call FINDVP, which calculates the VP by solution of a 
C (regularized) VanderMonde matrix, calculates the best fitting Gaussian by 
C means of FITGAUSS and CHI2H, and calculates the Gauss-Hermite coefficients
C by means of GAUHERM. The routine GRAMCHARCOF calculates the
C coefficients of the Gram-Charlier series.
C
C   VPANALYSE_FIX  : uses a fixed number of projected moments.
C                    Can be used if the regularization is sufficient.
C   VPANALYSE_CONV : Adds more or less projected momemnts until some
C                    sort of convergence is achieved. Works reasonable
C                    without regularization. The returned VP is nonsense,
C                    but the GH coefficients are reasonably OK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE VPANALYSE_FIX (iwhat,xi,nord,xmom,gam0,V0,sig0,darr,
     &      velar,vpval,nvel,gam,Vgau,sig,harr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C If iwhat=1, then do calculations at a given an angle xi from the major axis.
C If iwhat=3, then integrated over all angles in the first sky quadrant.
C In either case, calculate at projected radius 1
C (in dimensionless units), given an integer nord:
C   - the first nord projected moments, returned in xmom
C   - the lowest order true moments (gam0,V0,sig0) 
C   - the first nord Gram-Charlier coefficients, returned in darr
C   - the velocity profile, approximated at nord discrete velocities velar,
C        returned in vpval
C   - the parameters (gam,Vgau,sig) of the best fitting Gaussian
C   - the first nord Gauss-Hermite moments, returned in harr.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (maxlength = 2000000 ,
     &           mord      = 200     ,
     &           mev       = mord/2  )
C
      DIMENSION xmom(0:100),darr(0:100),velar(101),vpval(101),
     &          harr(0:100)
C
      COMMON /all/ allcofs(1:maxlength),
     &             icstart(0:mev,0:mev,0:mord),
     &             kmaxall(0:mev,0:mev,0:mord), maxord
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /projcase/ iproj
C
      COMMON /verbose/ iverb
C
C Determines whether verbose output should be generated for VP calculations
C
CCCCCCCCC
C
C Set nord to the nearest even integer
C
      nordin = nord
      IF (2*(nord/2).NE.nord) nord=nord-1
C
C Calculate the projected velocity moments
C
      IF (iverb.EQ.1) THEN
          WRITE (*,*) 'Calculating first', nord,
     &                ' projected velocity moments ...'
          WRITE (*,*) ' '
      END IF
C
      IF (iwhat.EQ.1) THEN
        rhopr0 = RHOPROJMOMAV(0)
        DO i=0,nord
          xmom(i) = PROJMOM(xi,i)
        END DO
      ELSE IF (iwhat.EQ.3) THEN
        rhopr0 = RHOPROJMOMAV(0)
        DO i=0,nord
          xmom(i) = RHOPROJMOMAV(i) / rhopr0
        END DO
      ELSE IF (iwhat.EQ.0) THEN
        rho0 = RHOVELMOM(xi,0,0,0)
        xmom(0) = 1.0D0
        DO i=1,nord
          IF (iproj.EQ.1) THEN
            xmom(i) = RHOVELMOM(xi,i,0,0) / rho0
          ELSE IF (iproj.EQ.2) THEN
            xmom(i) = RHOVELMOM(xi,0,i,0) / rho0
          ELSE
            xmom(i) = RHOVELMOM(xi,0,0,i) / rho0
          END IF
        END DO
      ELSE IF (iwhat.EQ.2) THEN
        rho0 = RHOVELMOMTHAV(0,0,0)
        xmom(0) = 1.0D0
        DO i=1,nord
          IF (iproj.EQ.1) THEN
            xmom(i) = RHOVELMOMTHAV(i,0,0) / rho0
          ELSE IF (iproj.EQ.2) THEN
            xmom(i) = RHOVELMOMTHAV(0,i,0) / rho0
          ELSE
            xmom(i) = RHOVELMOMTHAV(0,0,i) / rho0
          END IF
        END DO
      END IF
C
C Calculate the normalization, mean and dispersion
C
      gam0 = xmom(0)
      V0   = xmom(1)
      sig0 = SQRT(xmom(2)-(xmom(1)**2.0D0))
C
C Set initial values for the best fitting Gaussian
C
      gam  = gam0
      Vgau = V0
      sig  = sig0
C
C Calculate the VP, best fitting Gaussian, and Gauss-Hermite coefficients.
C
      CALL FINDVP (nord,xmom,velar,vpval,nvel,gam,Vgau,sig,harr)
C
C Write the resulting VP to the screen
C
      IF (iverb.EQ.1) THEN
C        WRITE (*,'(A31,I3,A9)')
C     &    'VP reconstructed from the first',nord,'moments:'
C        DO i=1,nvel
C          WRITE (*,'(I5,2F20.8)') i,velar(i),vpval(i)
C        END DO
C        WRITE (*,*) ' '
      END IF
C
C Calculate the Gram-Charlier coefficients
C
      CALL GRAMCHARCOF (xmom,nord,V0,sig0,darr)
C
C Reset nord to its input value
C
      nord = nordin
C
      END


      SUBROUTINE VPANALYSE_CONV (iwhat,xi,nord,xmom,gam0,V0,sig0,darr,
     &      velar,vpval,nvel,gam,Vgau,sig,harr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C If iwhat=1, then do calculations at a given an angle xi from the major axis.
C If iwhat=3, then integrated over all angles in the first sky quadrant.
C In either case, calculate at projected radius 1
C (in dimensionless units), given an integer nord:
C   - the first nord or more projected moments, returned in xmom
C   - the lowest order true moments (gam0,V0,sig0) 
C   - the first nord or more Gram-Charlier coefficients, returned in darr
C   - the velocity profile, approximated at nvel discrete velocities velar,
C        returned in vpval
C   - the parameters (gam,Vgau,sig) of the best fitting Gaussian
C   - the first nord Gauss-Hermite moments, returned in harr.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (maxlength = 2000000 ,
     &           mord      = 200     ,
     &           mev       = mord/2  )
C
      DIMENSION xmom(0:100),darr(0:100),velar(101),vpval(101),
     &          harr(0:100),hold(0:100)
C
      COMMON /all/ allcofs(1:maxlength),
     &             icstart(0:mev,0:mev,0:mord),
     &             kmaxall(0:mev,0:mev,0:mord), maxord
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /projcase/ iproj
C
      COMMON /verbose/ iverb
C
C Determines whether verbose output should be generated for VP calculations
C
CCCCCCCCC
C
C Calculate the lowest order moments
C
      IF (iwhat.EQ.1) THEN
        DO i=0,2
          xmom(i) = PROJMOM(xi,i)
        END DO
      ELSE IF (iwhat.EQ.3) THEN
        rhopr0 = RHOPROJMOMAV(0)
        DO i=0,2
          xmom(i) = RHOPROJMOMAV(i) / rhopr0
        END DO
      ELSE IF (iwhat.EQ.0) THEN
        rho0 = RHOVELMOM(xi,0,0,0)
        xmom(0) = 1.0D0
        DO i=1,2
          IF (iproj.EQ.1) THEN
            xmom(i) = RHOVELMOM(xi,i,0,0) / rho0
          ELSE IF (iproj.EQ.2) THEN
            xmom(i) = RHOVELMOM(xi,0,i,0) / rho0
          ELSE
            xmom(i) = RHOVELMOM(xi,0,0,i) / rho0
          END IF
        END DO
      ELSE IF (iwhat.EQ.2) THEN
        rho0 = RHOVELMOMTHAV(0,0,0)
        xmom(0) = 1.0D0
        DO i=1,2
          IF (iproj.EQ.1) THEN
            xmom(i) = RHOVELMOMTHAV(i,0,0) / rho0
          ELSE IF (iproj.EQ.2) THEN
            xmom(i) = RHOVELMOMTHAV(0,i,0) / rho0
          ELSE
            xmom(i) = RHOVELMOMTHAV(0,0,i) / rho0
          END IF
        END DO
      END IF
C
C Calculate the normalization, mean and dispersion
C
      gam0 = xmom(0)
      V0   = xmom(1)
      sig0 = SQRT(xmom(2)-(xmom(1)**2.0D0))
C
C Set initial values for the best fitting Gaussian
C
      gam  = gam0
      Vgau = V0
      sig  = sig0
C
C Now calculate all moments up to order minord. The choice 
C for minord is based on the criterium that delv (see FINDVP below)
C should be approximately equal to sig. It is never chosen smaller 
C than either nord or neversmaller.
C
      IF (ipot.EQ.1) THEN
        vmin = -1.0D0
        vmax = 1.0D0
      ELSE
        vmin = -3.0D0
        vmax = 3.0D0
      END IF
C
      neversmaller = 10
      minord = INT(1.0D0+((vmax-vmin)/(0.5*sig)))
      ncur   = MIN(maxord,MAX(minord,MAX(neversmaller,nord)))
C
C Use only even values for ncur
C
      IF (2*(ncur/2).NE.ncur) ncur=ncur-1
C
      IF (iverb.EQ.1) THEN
          WRITE (*,*) 'Calculating first', ncur,
     &                ' projected velocity moments ...'
          WRITE (*,*) ' '
      END IF
C
      DO i=3,ncur
        IF (iwhat.EQ.1) THEN   
          xmom(i) = PROJMOM(xi,i)
        ELSE IF (iwhat.EQ.3) THEN    
          xmom(i) = RHOPROJMOMAV(i) / rhopr0
        ELSE IF (iwhat.EQ.0) THEN
          IF (iproj.EQ.1) THEN
            xmom(i) = RHOVELMOM(xi,i,0,0) / rho0
          ELSE IF (iproj.EQ.2) THEN
            xmom(i) = RHOVELMOM(xi,0,i,0) / rho0
          ELSE
            xmom(i) = RHOVELMOM(xi,0,0,i) / rho0
          END IF
        ELSE IF (iwhat.EQ.2) THEN
          IF (iproj.EQ.1) THEN
            xmom(i) = RHOVELMOMTHAV(i,0,0) / rho0
          ELSE IF (iproj.EQ.2) THEN
            xmom(i) = RHOVELMOMTHAV(0,i,0) / rho0
          ELSE
            xmom(i) = RHOVELMOMTHAV(0,0,i) / rho0
          END IF
        END IF
      END DO
C
      gamold  = gam
      Vold    = Vgau
      sigold  = sig
      DO i=0,nord
        hold(i) = 0.0D0
      END DO
      delold  = 1.0D10
C
      istop = -1
C
C Calculate the VP, best fitting Gaussian, and Gauss-Hermite coefficients.
C Iterate until convergence is reached.
C
71    CALL FINDVP (ncur,xmom,velar,vpval,nvel,gam,Vgau,sig,harr)
C
C Write the resulting VP to the screen
C
      IF (iverb.EQ.1) THEN
C
C        WRITE (*,'(A31,I3,A9)')
C     &    'VP reconstructed from the first',ncur,'moments:'
C        DO i=1,nvel
C          WRITE (*,'(I5,2F20.8)') i,velar(i),vpval(i)
C        END DO
C        WRITE (*,*) ' '
C
        WRITE (*,'(A44,I3,A9)')
     &   'VP coefficients reconstructed from the first',ncur,'moments:'
        WRITE (*,'(15F8.4)') gam,Vgau,sig,(harr(k),k=0,nord)
        WRITE (*,*) ' '
C
      END IF
C
C If this was the first try, check for large negative values in the VP, 
C indicative of the fact that too many moments have been used.
C
      IF (istop.EQ.-1) THEN      
        ineg = 0
        DO i=1,nvel
          IF (vpval(i).LE.-0.5D0) ineg=1
        END DO
        IF ((ineg.EQ.1).AND.(ncur.GE.10)) THEN
          ncur = ncur-2
          GOTO 71
        END IF
      END IF
C
      IF (istop.NE.1) THEN
C
C Calculate how much the difference is from the previous estimate
C 
        del = (((gam-gamold)/gamold)**2.0D0) + 
     &        (((Vgau-Vold)/sigold)**2.0D0) +
     &        (((sig-sigold)/sigold)**2.0D0)    
        DO i=0,nord
          del = del + ((harr(i)-hold(i))**2.0D0)
        END DO
        del = SQRT(del/DBLE(nord+4))
C
        IF (iverb.EQ.1) THEN
          WRITE (*,'(A36,F12.8)') 
     &      'Change with respect to previous step',del
          WRITE (*,*) ' '
        END IF
C
        IF (del.GT.delold) THEN
C
C Stop whenever things start diverging (the problem of recovering the VP
C from the moments is numerically unstable to round-off error if to many
C moments are used).
C
          gam   = gamold
          Vgau  = Vold
          sig   = sigold
          IF (istop.EQ.-2) THEN
            ncur  = ncur-4
            delold  = 1.0D10
            istop = -1
          ELSE
            ncur  = ncur-2
            istop = 1
          END IF
          GOTO 71
C
        ELSE IF ((del.GE.1.0D-5).AND.(ncur+2.LE.maxord)) THEN
C
C Keep on going by including two more moments
C
          gamold  = gam
          Vold    = Vgau
          sigold  = sig
          DO i=0,nord
            hold(i) = harr(i)
          END DO
          delold  = del
C
          xmom(ncur+1) = PROJMOM(xi,ncur+1)
          xmom(ncur+2) = PROJMOM(xi,ncur+2)
          ncur = ncur+2
          IF (istop.EQ.-1) THEN
            istop = -2
          ELSE
            istop = 0
          END IF
          GOTO 71
C
        END IF
C
      END IF
C
C Calculate the Gram-Charlier coefficients
C
      CALL GRAMCHARCOF (xmom,ncur,V0,sig0,darr)
C
      END


      SUBROUTINE FINDVP (nord,xmom,velar,vpval,nvel,
     &                   gam,Vgau,sig,harr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Given nord projected moments xmom, calculate the VP as an array vpval
C at the nvel=nord+1 velocities velar, by solving the corresponding 
C vanderMonde matrix. Calculate the parameters 
C (gam,Vgau,sig) of the best fitting Gaussian (which must have been
C preset at appropriate initial guesses). The first nord Gauss-Hermite 
C moments are returned in harr.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION xmom(0:100),velar(101),vpval(101),harr(0:100)
      DIMENSION qmom(101),weight(101)
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /vpcur/ vvelar(101),vvpval(101),wweight(101),nvvel
C
C Common block in which the current velocity profile is stored, for use in
C REAL FUNCTION CHI2H.
C
      COMMON /reg/ xlam
C
      COMMON /verbose/ iverb
C
C Determines whether verbose output should be generated for VP calculations
C
CCCCCCCCCCCCCCCCC
C
C We will approximate integrals with quadrature formulae of the form:
C   INT f(x) dx = SUM_{n=0}^{N} w_n f(x_n) 
C where the {x_n} are a set of abscissa, and the w_n a set of weights.      
C
C For the Kepler potential the terminal velocities are -1 and 1. We
C choose Eulerian integration between these limits.
C For the Logarithmic potential we also use Eulerian integration, now
C with limits -3 and 3. Gaussian quadratures were tried, but didn't work
C satisfactorily.
C
      nvel  = nord+1
C
      IF (ipot.EQ.1) THEN
        vmin = -1.0D0
        vmax = 1.0D0
        delv = (vmax-vmin)/DBLE(nvel)
        DO i=1,nvel
          velar(i)  = vmin + ((DBLE(i)-0.5D0)*delv)
          weight(i) = delv
        END DO
      ELSE
        vmin = -3.0D0
        vmax = 3.0D0
        delv = (vmax-vmin)/DBLE(nvel)
        DO i=1,nvel
          velar(i)  = vmin + ((DBLE(i)-0.5D0)*delv)
          weight(i) = delv
        END DO
      END IF
C
C Define a vector with the moments that are to be reproduced
C
      DO i=1,nvel
        qmom(i) = xmom(i-1)
      END DO
C
C Now solve the corresponding Van der Monde matrix (Numerical 
C Recipes, eq. 2.8.2). If xlam<0, then do not use regularization.
C If xlam>0, then use regularization. If xlam=0, then the
C regularization parameter is determined iteratively, so as to yield 
C an acceptably smooth VP.
C
      IF (xlam.LT.0.0D0) THEN
        CALL VANDER(velar,vpval,qmom,nvel)
      ELSE IF (xlam.GT.0.0D0) THEN
        CALL REG_VANDER(velar,vpval,qmom,nvel,xlam)
      ELSE
        xlamreg = 1.0D-30
81      CALL REG_VANDER(velar,vpval,qmom,nvel,xlamreg) 
        CALL TESTSMOOTH(velar,vpval,nvel,iOK,Nmax)
        IF (iverb.EQ.1) THEN
          WRITE (*,'(A35,E15.7,I5)') 
     &     'Regularization, local maxima : ', xlamreg, Nmax
        END IF
        IF (iOK.EQ.0) THEN
          xlamreg = xlamreg * 10.0D0
          GOTO 81
        END IF
        IF (iverb.EQ.1) WRITE (*,*) ' '
      END IF
C
C And divide by the weights to get the approximations to the VP
C
      DO i=1,nvel
        vpval(i) = vpval(i)/weight(i)
      END DO
C
C Copy the results for use in FITGAUSS
C
      nvvel = nvel
      DO i=1,nvel
        vvpval(i)  = vpval(i)
        vvelar(i)  = velar(i)
        wweight(i) = weight(i)
      END DO
C
C Find the best fitting Gaussian (by searching for those values that
C come as closely as possible to generating h0=1, h1=h2=0).
C
      CALL FITGAUSS (gam,Vgau,sig)
C
C Now get all the Gauss-Hermite coefficients up to order nord
C
      CALL GAUHERM (velar,vpval,weight,nvel,gam,Vgau,sig,harr,nord)
C
      END


      SUBROUTINE TESTSMOOTH(velar,vpval,nvel,iOK,Nmax)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Test whether a recoverd VP is acceptably smooth. The test crtiterion is
C that there are no more than 3 `significant' local maxima in the VP. 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION velar(101),vpval(101)
C
      COMMON /smoo/ epsmoo
C
C Parameter that determines the requested smoothness, as fraction of
C the VP maximum.
C
CCCCCCCCCCCCCCCCCCCC
C
C Find the maximum of ABS(VP)
C 
      vpmax = -1.0D30
      DO i=1,nvel
        IF (ABS(vpval(i)).GT.vpmax) THEN
          vpmax = ABS(vpval(i))
        END IF  
      END DO
C
C Find the number of local maxima
C
      dmax = epsmoo * vpmax
      Nmax = 0
      iOK  = 1
C
      DO i=2,nvel-1
        dvm = vpval(i)-vpval(i-1)
        dvp = vpval(i)-vpval(i+1)
        IF ((dvm.GT.dmax).AND.(dvp.GT.dmax)) Nmax = Nmax + 1
      END DO
C
      IF (Nmax.GE.4) iOK=0
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
      COMMON /vpcur/ velar(101),vpval(101),weight(101),nvel
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
      DIMENSION velar(101),vpval(101),weight(101),harr(0:100)
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


      SUBROUTINE GRAMCHARCOF (xmom,nord,V0,sig0,darr)
C
C Calculate the Gram-Charlier coefficients of the velocity profile, given
C the array xmom, containing the first nord moments.
C The reference values (V0,sig0) of the series are required at 
C input separately. These need not necessarily correspond to the true
C mean and dispersion of the VP. The Gram-Charlier coefficients are
C calculated up to order nord, and are returned in darr.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION xmom(0:100),darr(0:100)
C
      DO l=0,nord
        dl = DBLE(l)
        darr(l) = 0.0D0
        DO j=0,l
          dj = DBLE(j)
          IF (2*((j+l)/2).EQ.(j+l)) THEN
            pfacln = (0.5D0*gammaln(dl+1.0D0)) - gammaln(dj+1.0D0) -
     &         gammaln(dl-dj+1.0D0) + 
     &         (gammaln(0.5D0*(dl-dj+1.0D0))) - gammaln(0.5D0) +
     &         (0.5D0*(dl-dj)*LOG(2.0D0))
            pfac = ((-1.0D0)**((l-j)/2)) * EXPP(pfacln)
            DO i=0,j
              darr(l) = darr(l) +
     &           ( ((-1.0D0)**(j-i)) * (V0**(j-i)) * 
     &             ((1.0/sig0)**j) * xmom(i) * pfac *
     &             EXPP(binomln(j,i)) )
            END DO
          END IF
        END DO
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


      REAL*8 FUNCTION RHOPROJMOMAV (nord)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The function RHOPROJMOMAV calculates the
C projected rho times velocity moment of a given order, integrated along
C a circle in the first quadrant of the sky.
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (pi=3.14159265358979D0)
C
      COMMON /howint/ iint
C
      COMMON /curprojorder/ nordcur
C
C Common block that holds copies of the order with which this 
C function is called.
C      
      EXTERNAL TOINTPROJ,MIDPNTB
C      
CCCCCCCCCCCCCCCCCCCC
C
      nordcur  = nord
C      
      IF (iint.EQ.0) THEN
        CALL QROMOB (TOINTPROJ,0.0D0,0.5D0*pi,SS,MIDPNTB)
      ELSE
        CALL QGAUSLEGB (TOINTPROJ,0.0D0,0.5D0*pi,SS)
      END IF
C
      RHOPROJMOMAV = SS / (0.5D0*pi)
C 
      END


      REAL*8 FUNCTION TOINTPROJ (xi)
C
C The function that must be integrated to get the first-quadrant integral
C over a circle on the sky of rho times projected velocity moment for the
C order given by /curprojorder/.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /curprojorder/ nordcur
C
      TOINTPROJ = RHOPROJ(xi) * PROJMOM(xi,nordcur)
C
      END
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The function RHOPROJ calculates the projected mass density analytically.
C The function PROJMOM calculates the projected velocity momemnts by
C numerical 1D integration over the function TOINT. This uses
C RHOVPROJMOM, rho times the n-th projected velocity moment (selected
C point in the galaxy.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION RHOPROJ (xi)
C
C The projected mass density at scaled projected radius equal to 1,
C and angle xi from the major axis. The expression is known
C analytically, because the 3D mass density is simply a power law on 
C spheroids.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /viewing/ xinc
C
      qp2 = (COS(xinc)**2.0D0) + ((q*q)*(SIN(xinc)**2.0D0))
      qp  = SQRT(qp2)
C
      xacc = COS(xi)
      yacc = SIN(xi)
C
      RHOPROJ = (q/qp) * EXPP(betaln(0.5D0,0.5D0*(gamma-1.0D0))) *
     &  ( ((xacc**2.0D0)+((yacc**2.0D0)/qp2))**
     &    (0.5D0*(1.0D0-gamma)) )
C
      END


      REAL*8 FUNCTION PROJMOM (xi,nord)
C
C Calculate the projected n-th order velocity moment on the sky,
C at scaled projected radius equal to 1, and angle xi from the
C major axis.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (pi=3.14159265358979D0)
C
      COMMON /losmom/ xxi,nnord
C
      COMMON /projcase/ iproj
C
C
      COMMON /howint/ iint
C
      EXTERNAL TOINT,MIDPNT
C
      IF ( ((2*(nord/2)).NE.nord) .AND.
     &     (ABS(COS(xi)).LE.1.0D-12) ) THEN
        PROJMOM = 0.0D0
        RETURN
      END IF
C
      IF (nord.EQ.0) THEN
        PROJMOM = 1.0D0
        RETURN
      END IF
C        
      xxi   = xi
      nnord = nord
C
      IF (iint.EQ.0) THEN
        CALL QROMO (TOINT,-0.5D0*pi,0.5D0*pi,SS,MIDPNT)
      ELSE
        CALL QGAUSLEG (TOINT,-0.5D0*pi,0.5D0*pi,SS)
      END IF
C
      PROJMOM = SS / RHOPROJ(xi)
C
      END


      REAL*8 FUNCTION TOINT (tau)
C
C The function that must be integrated to get the n-th order velocity 
C moment on the sky.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /losmom/ xi,nord
C
      COMMON /projcase/ iproj
C
C
      COMMON /potent/ ipot
C
C ipot=1 for Kepler potential and ipot=2 for logarithmic potential
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /viewing/ xinc
C
      xacc = COS(xi)
      yacc = SIN(xi)
      zacc = TAN(tau)
C
      x = (-1.0D0*COS(xinc)*yacc) + (SIN(xinc)*zacc)
      y = xacc
      z = (SIN(xinc)*yacc) + (COS(xinc)*zacc)
      Rc = SQRT((x**2.0D0)+(y**2.0D0))
      r = SQRT((Rc**2.0D0)+(z**2.0D0))
C
      theta = ACOS(z/r)
      phi   = ACOS(x/Rc)
C
      IF (ipot.EQ.1) THEN
        TOINT = (COS(tau)**(gamma-2.0D0+(0.5D0*DBLE(nord)))) *
     &        RHOVPROJMOM(theta,phi,nord)
      ELSE IF (ipot.EQ.2) THEN
        TOINT = (COS(tau)**(gamma-2.0D0)) *
     &        RHOVPROJMOM(theta,phi,nord)
      ELSE
        STOP 'ipot wrong value'
      END IF
C
      END


      REAL*8 FUNCTION RHOVLOSMOM(theta,phi,nord)
C
C Returns rho times the nord-th line-of-sight velocity moment at the
C point with polar coordinates (theta,phi) in radians, and radius 1
C (in dimensionless units).
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /viewing/ xinc
C
      Afac = (COS(phi)*SIN(theta)*SIN(xinc)) + (COS(theta)*COS(xinc))
      Bfac = (COS(phi)*COS(theta)*SIN(xinc)) - (SIN(theta)*COS(xinc))
      Cfac = (-1.0D0*SIN(phi)*SIN(xinc))
C
      IF (Afac.LT.0.0D0) THEN
        Asign = -1.0D0        
      ELSE
        Asign = 1.0D0        
      END IF
      Afac = MAX(1.0D-20,ABS(Afac))
C
      IF (Bfac.LT.0.0D0) THEN
        Bsign = -1.0D0        
      ELSE
        Bsign = 1.0D0        
      END IF
      Bfac = MAX(1.0D-20,ABS(Bfac))
C
      IF (Cfac.LT.0.0D0) THEN
        Csign = -1.0D0        
      ELSE
        Csign = 1.0D0        
      END IF
      Cfac = MAX(1.0D-20,ABS(Cfac))
C
      RHOVLOSMOM = 0.0D0
      DO k=0,nord
        DO j=0,nord-k
          facln = binomln(nord,k) + binomln(nord-k,j) +
     &            (DBLE(k)*LOG(Afac)) + (DBLE(j)*LOG(Bfac)) +
     &            (DBLE(nord-k-j)*LOG(Cfac))
          RHOVLOSMOM = RHOVLOSMOM + ( EXPP(facln)*
     &       (Asign**k)*(Bsign**j)*(Csign**(nord-k-j))*
     &       RHOVELMOM(theta,k,j,nord-k-j) )
        END DO
      END DO
C
      END


      REAL*8 FUNCTION RHOVPOSRMOM(theta,phi,nord)
C
C Returns rho times the nord-th line-of-sight velocity moment at the
C point with polar coordinates (theta,phi) in radians, and radius 1
C (in dimensionless units).
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /viewing/ xinc
C
      SP = SIN(phi)
      CP = COS(phi)
      ST = SIN(theta)
      CT = COS(theta)
      SX = SIN(xinc)
      CX = COS(xinc)
      CTT = COS(2.0D0*theta)
C
      Afac = SQRT(ST*ST*SP*SP+(CT*SX - ST*CX*CP)**2.0D0)
      Bfac = (ST*CT*CX*CX*CP*CP + ST*CT*(SP*SP - SX*SX) -
     &       CTT*SX*CX*CP)/Afac
      Cfac = SP*(ST*CP + ST*CX*CX*(-1.0D0*CP) + 
     &       CT*SX*CX)/Afac
C
      IF (Afac.LT.0.0D0) THEN
        Asign = -1.0D0        
      ELSE
        Asign = 1.0D0        
      END IF
      Afac = MAX(1.0D-20,ABS(Afac))
C
      IF (Bfac.LT.0.0D0) THEN
        Bsign = -1.0D0        
      ELSE
        Bsign = 1.0D0        
      END IF
      Bfac = MAX(1.0D-20,ABS(Bfac))
C
      IF (Cfac.LT.0.0D0) THEN
        Csign = -1.0D0        
      ELSE
        Csign = 1.0D0        
      END IF
      Cfac = MAX(1.0D-20,ABS(Cfac))
C
      RHOVPOSRMOM = 0.0D0
      DO k=0,nord
        DO j=0,nord-k
          facln = binomln(nord,k) + binomln(nord-k,j) +
     &            (DBLE(k)*LOG(Afac)) + (DBLE(j)*LOG(Bfac)) +
     &            (DBLE(nord-k-j)*LOG(Cfac))
          RHOVPOSRMOM = RHOVPOSRMOM + ( EXPP(facln)*
     &       (Asign**k)*(Bsign**j)*(Csign**(nord-k-j))*
     &       RHOVELMOM(theta,k,j,nord-k-j) )
        END DO
      END DO
C
      END



      REAL*8 FUNCTION RHOVPOSTMOM(theta,phi,nord)
C
C Returns rho times the nord-th line-of-sight velocity moment at the
C point with polar coordinates (theta,phi) in radians, and radius 1
C (in dimensionless units).
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /param/ gamma, beta, q, alpha, eta
      COMMON /viewing/ xinc
C
      SP = SIN(phi)
      CP = COS(phi)
      ST = SIN(theta)
      CT = COS(theta)
      SX = SIN(xinc)
      CX = COS(xinc)
      DEN = SQRT(ST*ST*SP*SP+(CT*SX - ST*CX*CP)**2.0D0)
C
      Bfac = SX*SP/DEN
      Cfac = (CT*SX*CP - ST*CX)/DEN
C
      IF (Bfac.LT.0.0D0) THEN
        Bsign = -1.0D0        
      ELSE
        Bsign = 1.0D0        
      END IF
      Bfac = MAX(1.0D-20,ABS(Bfac))
C
      IF (Cfac.LT.0.0D0) THEN
        Csign = -1.0D0        
      ELSE
        Csign = 1.0D0        
      END IF
      Cfac = MAX(1.0D-20,ABS(Cfac))
C
      RHOVPOSTMOM = 0.0D0
      DO j=0,nord
            facln = binomln(nord,0) + binomln(nord,j) +
     &      (DBLE(j)*LOG(Bfac)) +
     &      (DBLE(nord-j)*LOG(Cfac))
            RHOVPOSTMOM = RHOVPOSTMOM + ( EXPP(facln)*
     &      (Bsign**j)*(Csign**(nord-j))*
     &      RHOVELMOM(theta,0,j,nord-j) )
      END DO
C
      END


      REAL*8 FUNCTION RHOVPROJMOM(theta,phi,nord)
C
C Select the projected component (LOS=1, POSR=2, POST=3)
C and return rho times the nord-th projected velocity moment at the
C point with polar coordinates (theta,phi) in radians, and radius 1
C (in dimensionless units).
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /projcase/ iproj
C
      IF (iproj.EQ.1) THEN
        RHOVPROJMOM = RHOVLOSMOM(theta,phi,nord)
      ELSE IF (iproj.EQ.2) THEN
        RHOVPROJMOM = RHOVPOSRMOM(theta,phi,nord)
      ELSE IF (iproj.EQ.3) THEN
        RHOVPROJMOM = RHOVPOSTMOM(theta,phi,nord)
      ELSE
        STOP 'iproj wrong value'
      END IF
C
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The subroutine VIRIALTENSORS calculates the virial tensors, using
C the fact that the intrinsic moments are power series in e^2 SIN^2(theta).
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE VIRIALTENSORS (potx,potz,xkinr,xkinth,xkinph,
     &                          xkinRc,xkinz,xkinx)
C
C Calculate the tensors for potential and kinemtic energy, and write
C the results to the screen.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /param/ gamma, beta, q, alpha, eta
C
      e2 = 1.0D0 - (q*q)
      e2 = MAX(e2,1.0D-20)
C
      potx   = 0.0D0
      potz   = 0.0D0
      prefln = PREFACLN(0,0,0) 
      DO k=0,KMAXARR(0,0,0)
        dk     = DBLE(k)
        addxln = prefln + COFARR(k,0,0,0) + (dk*LOG(e2)) +
     &           betaln(2.0D0+dk,0.5D0) 
        addzln = prefln + COFARR(k,0,0,0) + (dk*LOG(e2)) +
     &           betaln(1.0D0+dk,1.5D0) 
        potx = potx - (0.5D0*EXPP(addxln))
        potz = potz - EXPP(addzln)
      END DO
C
      xkinr  = 0.0D0
      xkinr1 = 0.0D0
      xkinr2 = 0.0D0
      prefln = PREFACLN(2,0,0)
      DO k=0,KMAXARR(2,0,0)
        dk     = DBLE(k)
        addln  = prefln + COFARR(k,2,0,0) + (dk*LOG(e2)) +
     &           betaln(1.0D0+dk,0.5D0) 
        addln1 = prefln + COFARR(k,2,0,0) + (dk*LOG(e2)) +
     &           betaln(2.0D0+dk,0.5D0) 
        addln2 = prefln + COFARR(k,2,0,0) + (dk*LOG(e2)) +
     &           betaln(1.0D0+dk,1.5D0) 
        xkinr  = xkinr + EXPP(addln)
        xkinr1 = xkinr1 + EXPP(addln1)
        xkinr2 = xkinr2 + EXPP(addln2)
      END DO
C
      xkinth  = 0.0D0
      xkinth1 = 0.0D0
      xkinth2 = 0.0D0
      prefln = PREFACLN(0,2,0)
      DO k=0,KMAXARR(0,2,0)
        dk     = DBLE(k)
        addln  = prefln + COFARR(k,0,2,0) + (dk*LOG(e2)) +
     &           betaln(1.0D0+dk,0.5D0) 
        addln1 = prefln + COFARR(k,0,2,0) + (dk*LOG(e2)) +
     &           betaln(2.0D0+dk,0.5D0) 
        addln2 = prefln + COFARR(k,0,2,0) + (dk*LOG(e2)) +
     &           betaln(1.0D0+dk,1.5D0) 
        xkinth  = xkinth  + EXPP(addln)
        xkinth1 = xkinth1 + EXPP(addln1)
        xkinth2 = xkinth2 + EXPP(addln2)
      END DO
C
      xkinph = 0.0D0
      prefln = PREFACLN(0,0,2)
      DO k=0,KMAXARR(0,0,2)
        dk    = DBLE(k)
        addln = prefln + COFARR(k,0,0,2) + (dk*LOG(e2)) +
     &          betaln(1.0D0+dk,0.5D0) 
        xkinph = xkinph + EXPP(addln)
      END DO
C
      xkinRc = xkinr1 + xkinth2
      xkinz  = xkinr2 + xkinth1
      xkinx  = 0.5D0 * (xkinRc+xkinph)
C
      END


      REAL*8 FUNCTION RHOVELMOMTHAV(lr,lth,lph)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The function RHOVELMOMTHAV calculates rho times the
C intrinsic velocity moment of a given order, integrated over
C a spherical shell of unit radius.
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (pi=3.14159265358979D0)
C
      COMMON /howint/ iint
C
      COMMON /curorders/ lrcur,lthcur,lphcur
C
C Common block that holds copies of the orders with which this 
C function is called.
C      
      EXTERNAL TOINTMER,MIDPNT
C      
CCCCCCCCCCCCCCCCCCCC
C
      lrcur  = lr
      lthcur = lth
      lphcur = lph
C      
      IF (iint.EQ.0) THEN
        CALL QROMO (TOINTMER,0.0D0,0.5D0*pi,SS,MIDPNT)
      ELSE
        CALL QGAUSLEG (TOINTMER,0.0D0,0.5D0*pi,SS)
      END IF
C
      RHOVELMOMTHAV = SS
C 
      END


      REAL*8 FUNCTION TOINTMER (theta)
C
C The function that must be integrated to get the integral over a spherical
C shell of rho times the intrinsic velocity moment of the orders given by
C /curorders/.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /curorders/ lrcur,lthcur,lphcur
C
      TOINTMER = SIN(theta) * RHOVELMOM(theta,lrcur,lthcur,lphcur)
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
C A regularized vanderMonde matrix solver
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE REG_VANDER (x,w,q,n,xlam)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Solves a VanderMonde matrix with regularization. 
C Arguments are as in de subroutine VANDER from Numerical Recipes.
C The input variable xlam is the regularization parameter. 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      INTEGER n,NMAX
      DIMENSION q(n),w(n),x(n)
C
      PARAMETER (NMAX=101)
C
      DIMENSION A(NMAX,NMAX), Anorm(NMAX,NMAX), H(NMAX,NMAX),
     &          Areg(NMAX,NMAX), bnorm(NMAX), indx(NMAX)
C
CCCCCCCCCCCCCCCCCCCC
C
C Set up the matrix that occurs in the matrix equation that is to be solved.
C
      DO i=1,n
        DO j=1,n
          A(i,j) = x(j)**(i-1)
        END DO
      END DO
C
C Calculate the matrix Atr*A that occurs in the normal equations.
C
      DO i=1,n
        DO j=1,n
          Anorm(i,j) = 0.0D0
          DO k=1,n
            Anorm(i,j) = Anorm(i,j) + (A(k,i)*A(k,j))
          END DO
        END DO
      END DO
C
C Calculate the right-hand side vector Atr*q that occurs in the normal 
C equations.
C
      DO i=1,n
        bnorm(i) = 0.0D0
        DO k=1,n
          bnorm(i) = bnorm(i) + (A(k,i)*q(k))
        END DO
      END DO
C
C Set up the desired regularization matrix H
C
      CALL REGMAT2 (H,n,NMAX)
C
C Calculate the equipartition xlambda
C
      Anormtr = 0.0D0
      Htr = 0.0D0
      DO i=1,n
        Anormtr = Anormtr + Anorm(i,i)
        Htr = Htr + H(i,i)
      END DO
C
      xlameq = Anormtr / Htr
C
C Construct the regularized matrix
C 
      DO i=1,n
        DO j=1,n
          Areg(i,j) = Anorm(i,j) + (xlam*xlameq*H(i,j))
        END DO
      END DO
C
C The system to be solved is now: Areg*w = bnorm. The solution is found 
C with LU decomposition.
C
      CALL LUDCMP (Areg,n,NMAX,indx,d)
      CALL LUBKSB (Areg,n,NMAX,indx,bnorm)
C
C The answer returned in bnorm is assigned to the vector w
C
      DO i=1,n
        w(i) = bnorm(i)
      END DO
C
      END


      SUBROUTINE REGMAT0 (H,n,np)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Set up a zeroth order regularization matrix
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION H(np,np)
C      
      IF (n.LE.2) STOP 'Matrix too small in REGMAT0'
C
      DO i=1,n
        DO j=1,n
          H(i,j) = 0.0D0
        END DO
      END DO
C
      H(1,1) = 1.0D0
      H(1,2) = -1.0D0
C
      DO i=2,n-1
        H(i,i-1) = -1.0D0
        H(i,i)   = 2.0D0
        H(i,i+1) = -1.0D0
      END DO
C
      H(n,n-1) = -1.0D0
      H(n,n)   = 1.0D0
C
      END


      SUBROUTINE REGMAT1 (H,n,np)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Set up a first order regularization matrix
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION H(np,np)
C      
      IF (n.LE.4) STOP 'Matrix too small in REGMAT1'
C
      DO i=1,n
        DO j=1,n
          H(i,j) = 0.0D0
        END DO
      END DO
C
      H(1,1) = 1.0D0
      H(1,2) = -2.0D0
      H(1,3) = 1.0D0
C
      H(2,1) = -2.0D0
      H(2,2) = 5.0D0
      H(2,3) = -4.0D0
      H(2,4) = 1.0D0
C
      DO i=3,n-2
        H(i,i-2) = 1.0D0
        H(i,i-1) = -4.0D0
        H(i,i)   = 6.0D0
        H(i,i+1) = -4.0D0
        H(i,i+2) = 1.0D0
      END DO
C
      H(n-1,n-3) = 1.0D0
      H(n-1,n-2) = -4.0D0
      H(n-1,n-1) = 5.0D0
      H(n-1,n)   = -2.0D0
C
      H(n,n-2) = 1.0D0
      H(n,n-1) = -2.0D0
      H(n,n)   = 1.0D0
C
      END


      SUBROUTINE REGMAT2 (H,n,np)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Set up a second order regularization matrix
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      DIMENSION H(np,np)
C      
      IF (n.LE.6) STOP 'Matrix too small in REGMAT2'
C
      DO i=1,n
        DO j=1,n
          H(i,j) = 0.0D0
        END DO
      END DO
C
      H(1,1) = 1.0D0
      H(1,2) = -3.0D0
      H(1,3) = 3.0D0
      H(1,4) = -1.0D0
C
      H(2,1) = -3.0D0
      H(2,2) = 10.0D0
      H(2,3) = -12.0D0
      H(2,4) = 6.0D0
      H(2,5) = -1.0D0
C
      H(3,1) = 3.0D0
      H(3,2) = -12.0D0
      H(3,3) = 19.0D0
      H(3,4) = -15.0D0
      H(3,5) = 6.0D0
      H(3,6) = -1.0D0
C
      DO i=4,n-3
        H(i,i-3) = -1.0D0
        H(i,i-2) = 6.0D0
        H(i,i-1) = -15.0D0
        H(i,i)   = 20.0D0
        H(i,i+1) = -15.0D0
        H(i,i+2) = 6.0D0
        H(i,i+3) = -1.0D0
      END DO
C
      H(n-2,n-5) = -1.0D0
      H(n-2,n-4) = 6.0D0
      H(n-2,n-3) = -15.0D0
      H(n-2,n-2) = 19.0D0
      H(n-2,n-1) = -12.0D0
      H(n-2,n)   = 3.0D0
C
      H(n-1,n-4) = -1.0D0
      H(n-1,n-3) = 6.0D0
      H(n-1,n-2) = -12.0D0
      H(n-1,n-1) = 10.0D0
      H(n-1,n)   = -3.0D0
C
      H(n,n-3) = -1.0D0
      H(n,n-2) = 3.0D0
      H(n,n-1) = -3.0D0
      H(n,n)   = 1.0D0
C
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


      SUBROUTINE qgausleg(func,a,b,ss)
C
C Modified version of qgaus. Common block /gleg/ must have been 
C filled previously
C
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 a,b,ss,func
      EXTERNAL func
      COMMON /gleg/ qx(300),qw(300),nGL
      ss = 0.0D0
      DO i=1,nGL
        ss = ss + (qw(i)*func(a+(qx(i)*(b-a))))
      END DO
      ss = ss*(b-a)
      END


      SUBROUTINE gauleg(x1,x2,x,w,n)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      PARAMETER (EPS=3.0D-14)
      INTEGER i,j,m
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END


      SUBROUTINE qromo(func,a,b,ss,choose)
C
C Modified to be double precision, and to receive eps from common block
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      COMMON /romoeps/ eps
      EXTERNAL func,choose
      PARAMETER (JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.0D0
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.0D0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0D0
11    continue
      STOP 'too many steps in qromo'
      END


      SUBROUTINE trapzd(func,a,b,s,n)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5D0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5D0*del
        sum=0.0D0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5D0*(s+(b-a)*sum/tnm)
      endif
      return
      END


      SUBROUTINE polint(xa,ya,n,x,y,dy)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0D0)STOP 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END


      SUBROUTINE midpnt(func,a,b,s,n)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5D0*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.0D0*tnm)
        ddel=del+del
        x=a+0.5D0*del
        sum=0.0D0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.0D0
      endif
      return
      END


      SUBROUTINE vander(x,w,q,n)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n,NMAX
      DIMENSION q(n),w(n),x(n)
      PARAMETER (NMAX=101)
      INTEGER i,j,k,k1
      DOUBLE PRECISION b,s,t,xx,c(NMAX)
      if(n.eq.1)then
        w(1)=q(1)
      else
        do 11 i=1,n
          c(i)=0.d0
11      continue
        c(n)=-x(1)
        do 13 i=2,n
          xx=-x(i)
          do 12 j=n+1-i,n-1
            c(j)=c(j)+xx*c(j+1)
12        continue
          c(n)=c(n)+xx
13      continue
        do 15 i=1,n
          xx=x(i)
          t=1.d0
          b=1.d0
          s=q(n)
          k=n
          do 14 j=2,n
            k1=k-1
            b=c(k)+xx*b
            s=s+q(k1)*b
            t=xx*t+b
            k=k1
14        continue
          w(i)=s/t
15      continue
      endif
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Numerical Recipes routines, transformed to double precision, for
C Regularized VanderMonde matrix solution.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ludcmp(a,n,np,indx,d)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0D-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.0D0
      do 12 i=1,n
        aamax=0.0D0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.0D0) STOP 'singular matrix in ludcmp'
        vv(i)=1.0D0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.0D0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.0D0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.0D0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.0D0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C All routines below are obsolete. However, they were kept because they
C might turn out to be useful one day. 
C
C The function densradius defines the relation r(theta) for an isodensity
C surface. The SUBROUTINE CONTOUR_ANALYSIS finds the best-fitting ellipse
C and does a Fourier analysis. In doing so it uses the Numerical Recipes
C routines cosft1, four1 and realft, which were kept at single precision.
C
C The calling sequence from the main program is:
C     EXTERNAL DENSRADIUS
C     CALL CONTOUR_ANALYSIS (DENSRADIUS)
C
C In our model the density is always perfectly spheroidal, so the results
C are rather uninteresting. However, CONTOUR_ANALYSIS can be used for any
C function r(theta), so it could be used, e.g., to determine the shape
C of the contours of the projected velocity dispersion on the sky.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION DENSRADIUS (theta)
C
C The radius r of a point at polar angle theta, on the isodensity
C surface with rho = 1.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      COMMON /param/ gamma, beta, q, alpha, eta
C
      DENSRADIUS = RHOVELMOM(theta,0,0,0)**(1.0/gamma)
C
      END


      SUBROUTINE CONTOUR_ANALYSIS (RADIUS)
C
C Given a function RADIUS, that returns a radius as function of polar angle 
C theta, do a shape analysis of the contour defined by this function. Four-fold
C symmetry is assumed.
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (pi=3.14159265358979D0)
C
      EXTERNAL RADIUS
C
      DIMENSION radarr(10000)
      REAL radreal(10000)
C
C Estimate the axis ratio by taking the ratio of points on the
C major and minor axis
C
      qest = RADIUS(0.0D0)/RADIUS(pi/2.0D0)
C
C We now write the radius to an array and Fourier analyse it. 
C Rather than using the conventional polar coordinates (r,theta), we
C use elliptical coordinates, such that 
C   R = s sin(tau)
C   z = s qcoord cos(tau)
C
      N2pow  = 256
      qcoord = qest
C
25    DO i=1,N2pow+1
C
C Choose tau, spaced linearly between 0 and pi/2.
C
        tau = ((pi/2.0D0)/(DBLE(N2pow)))*DBLE(i-1)
C
C Calculate the correspinding theta
C
        theta = ATAN(TAN(tau)/qcoord)
        st = SIN(theta)
        ct = COS(theta)
C
C Calculate the corresponding radius
C
        radthet = RADIUS(theta)
C
C Calculate and store the corresponding elliptical radius
C
        radarr(i)  = radthet * 
     &     SQRT((st**2.0D0)+((ct/qcoord)**2.0D0))
        radreal(i) = REAL(radarr(i))
C
      END DO
C
C Get the Fourier cosine decomposition of the array radarr. If the 
C contour is a perfect ellipse, the only the lowest order (constant)
C term in the transfrom is non-zero. If the cos (2 tau) is non-zero, then 
C we should choose a different qcoord, and iterate.
C 
      CALL COSFT1(radreal,N2pow)
C
C radarr now contains the Fourier cosine series, which we properly normalize,
C and make double precision again.
C
      semimajor = DBLE(radreal(1))/DBLE(N2pow)
      radarr(1) = 1.0D0
      DO i=2,N2pow+1
        radarr(i) = 2.0D0*DBLE(radreal(i))/(DBLE(N2pow)*semimajor)
      END DO
C
C No verbose output on the screen please
C
      iverb = 0
      IF (iverb.EQ.1) THEN
        WRITE (*,'(A5,F12.8)') ' ', semimajor 
        DO i=1,5
          WRITE (*,'(I5,F12.8)') 2*(i-1),radarr(i)
        END DO
        WRITE (*,*) ' '
      END IF
C
C If the COS(2 tau) coefficient is too big, then change qcoord. 
C The formula for doing this is obtained by writing the equations
C for a pure ellipse with axial ratio q, in elliptical coordinates 
C with axial ratio q'. After Fourier expanding the equation s'(tau'),
C the cos(2 tau') coefficient is to lowest order directly related to (q/q').
C
      IF (ABS(radarr(2)).GE.1.0D-8) THEN
        qcoord = qcoord * (1.0D0 + (2.0D0*radarr(2)))
        qcoord = MAX(0.0D0,MIN(qcoord,1.0D0))
        GOTO 25
      END IF
C
C Write results to the screen
C
      WRITE (*,'(A45,F15.8)') 'Ratio major/minor of contour',qest
      WRITE (*,'(A45,F15.8)') 'axial ratio best fitting ellipse',
     &         qcoord
      WRITE (*,'(A45,F15.8)') 'semi-major axis best fitting ellipse',
     &         semimajor
      WRITE (*,'(A45,F15.8)') 'cos(4 theta) coefficient',
     &         radarr(3)/radarr(1)
      WRITE (*,'(A45,F15.8)') 'cos(6 theta) coefficient',
     &         radarr(4)/radarr(1)
      WRITE (*,'(A45,F15.8)') 'cos(8 theta) coefficient',
     &         radarr(5)/radarr(1)
C
      END


      SUBROUTINE realft(data,n,isign)
      INTEGER isign,n
      REAL data(n)
CU    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      END


      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END


      SUBROUTINE cosft1(y,n)
      INTEGER n
      REAL y(n+1)
CU    USES realft
      INTEGER j
      REAL sum,y1,y2
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/n
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      sum=0.5*(y(1)-y(n+1))
      y(1)=0.5*(y(1)+y(n+1))
      do 11 j=1,n/2-1
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=0.5*(y(j+1)+y(n-j+1))
        y2=(y(j+1)-y(n-j+1))
        y(j+1)=y1-wi*y2
        y(n-j+1)=y1+wi*y2
        sum=sum+wr*y2
11    continue
      call realft(y,n,+1)
      y(n+1)=y(2)
      y(2)=sum
      do 12 j=4,n,2
        sum=sum+y(j)
        y(j)=sum
12    continue
      return
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Copies of routines already given above, but with slightly modified names,
C to avoid recursive calling when evaluating double integrals.      
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE qgauslegB(func,a,b,ss)
C
C Modified version of qgaus. Common block /gleg/ must have been 
C filled previously
C
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 a,b,ss,func
      EXTERNAL func
      COMMON /gleg/ qx(300),qw(300),nGL
      ss = 0.0D0
      DO i=1,nGL
        ss = ss + (qw(i)*func(a+(qx(i)*(b-a))))
      END DO
      ss = ss*(b-a)
      END


      SUBROUTINE qromoB(func,a,b,ss,choose)
C
C Modified to be double precision, and to receive eps from common block
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      COMMON /romoeps/ eps
      EXTERNAL func,choose
      PARAMETER (JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polintB
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.0D0
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polintB(h(j-KM),s(j-KM),K,0.0D0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.0D0
11    continue
      STOP 'too many steps in qromo'
      END


      SUBROUTINE polintB(xa,ya,n,x,y,dy)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0D0)STOP 'failure in polintB'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END


      SUBROUTINE midpntB(func,a,b,s,n)
C
C Modified to be double precision
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5D0*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.0D0*tnm)
        ddel=del+del
        x=a+0.5D0*del
        sum=0.0D0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.0D0
      endif
      return
      END

      INTEGER FUNCTION LENSTR(str)
C
C Return effective length of a character string (trailing blanks removed).
C Returns 0 for an all-blank string.
C
      IMPLICIT REAL*8 (a-h,o-z)
      CHARACTER*(*) str
      INTEGER i
      LENSTR = LEN(str)
      DO i=LEN(str),1,-1
        IF (str(i:i).NE.' ') THEN
          LENSTR = i
          RETURN
        END IF
      END DO
      LENSTR = 0
      RETURN
      END
