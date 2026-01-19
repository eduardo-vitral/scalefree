program mcdraw_gh
!----------------------------------------------------------------------
! MCdraw_gh.f90
!
! Scale-free axisymmetric mock generator (MC positions + GH velocity sampling)
!
! Behavior:
! - Draws (r,theta,phi) exactly as in legacy MCdraw5 (same rejection sampling).
! - Computes intrinsic 2nd moments via RHOVELMOM (legacy routines).
! - Delegates velocity sampling (vr,vth,vph) to Python helper gh_velocity_sampler.py
!   using (mean,sigma,h3,h4) per component (currently h3=h4=0 by default).
! - Emits mock to STDOUT:
!     # kind=mock
!     # columns: x y z vx vy vz
!   and summary stats:
!     # kind=mock_stats
!     # columns: binney_beta rms_method1 rms_method2 mean_vphi
! - Deletes temporary exchange files (best-effort) via Fortran CLOSE(...,STATUS='DELETE')
!   and reports if deletion fails.
!
! Build:
!   gfortran -O2 -std=legacy -c mcdraw_legacy.f -o mcdraw_legacy.o
!   gfortran -O2 -std=f2008 MCdraw_gh.f90 mcdraw_legacy.o -o MCdraw_gh.e
!----------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  real(dp), parameter :: pi = 3.14159265358979_dp

  ! --- legacy common blocks expected by mcdraw_legacy.f ---
  integer :: ipot, icase
  real(dp) :: gamma, beta, q, alpha, eta
  common /potent/ ipot
  common /DFcase/ icase
  common /param/  gamma, beta, q, alpha, eta

  ! --- external routines from legacy file ---
  real(dp), external :: ran1, rhovelmom
  external :: fillarrays

  ! --- model / MC parameters ---
  integer :: nsamp, idum
  real(dp) :: delta, ri, ro

  ! --- per-sample quantities ---
  integer :: j
  real(dp) :: lr, thcur, phcur, rcur
  real(dp) :: costh, sinth, cosph, sinph
  real(dp) :: rho, vr2, vth2, vph2, vphav
  real(dp) :: xtmp, ytmp, ztmp

  ! --- variables for MCdraw5-matching theta sampler ---
  real(dp) :: st, th, ct, xi, xxi, fmass

  ! --- GH params ---
  real(dp) :: vr_mean, vr_sig, vr_h3, vr_h4
  real(dp) :: vth_mean, vth_sig, vth_h3, vth_h4
  real(dp) :: vph_mean, vph_sig, vph_h3, vph_h4

  ! --- arrays ---
  real(dp), allocatable :: x(:), y(:), z(:)
  real(dp), allocatable :: vr(:), vth(:), vph(:)
  real(dp), allocatable :: vx(:), vy(:), vz(:)

  ! --- exchange files and python ---
  character(len=512) :: in_path, out_path, pyexe
  integer :: uin, uout, stat

  ! --- deletion helpers ---
  integer :: udel, ios_del

  ! --- summary stats ---
  real(dp) :: vr2sum, vt2sum, v2sum, vp1sum, betabin

  write(*,*) 'Please answer the following questions about'
  write(*,*) 'the parameters of the model:'
  write(*,*) ' '

  write(*,*) 'Kepler (1) or Logarithmic (2) Potential ?'
  read(*,*) ipot
  write(*,*) ' '

  delta = 2.0_dp - real(ipot, dp)

  write(*,*) 'Power-law slope gamma of the mass density'
  read(*,*) gamma
  write(*,*) ' '

  write(*,*) 'Intrinsic axial ratio q of the mass density'
  read(*,*) q
  write(*,*) ' '

  write(*,*) 'Case I (1) or Case II (2) DF ?'
  read(*,*) icase
  write(*,*) ' '

  write(*,*) 'Anisotropy parameter beta of the DF'
  read(*,*) beta
  write(*,*) ' '

  write(*,*) 'Odd part parameter s for the DF'
  read(*,*) eta
  write(*,*) ' '

  write(*,*) 'Odd part parameter t for the DF'
  read(*,*) alpha
  write(*,*) ' '

  write(*,*) 'Number of particles to be generated'
  read(*,*) nsamp
  write(*,*) ' '

  write(*,*) 'Seed for the random sequence'
  read(*,*) idum
  if (idum > 0) idum = -idum
  write(*,*) ' '

  write(*,*) 'Inner radius for Monte-Carlo drawings'
  read(*,*) ri
  write(*,*) ' '

  write(*,*) 'Outer radius for Monte-Carlo drawings'
  read(*,*) ro
  write(*,*) ' '

  ! Precompute velocity moment series coefficients (legacy)
  call fillarrays(20)

  allocate(x(nsamp), y(nsamp), z(nsamp))
  allocate(vr(nsamp), vth(nsamp), vph(nsamp))
  allocate(vx(nsamp), vy(nsamp), vz(nsamp))

  ! Temporary exchange files (internal plumbing; user does not provide output path)
  in_path  = 'mcdraw_gh_input.txt'
  out_path = 'mcdraw_gh_output.txt'

  call get_environment_variable('PYTHON', pyexe, status=stat)
  if (stat /= 0 .or. len_trim(pyexe) == 0) pyexe = 'python'

  open(newunit=uin, file=trim(in_path), status='replace', action='write')
  write(uin,'(A)') '# kind=mcdraw_gh_input'
  write(uin,'(A)') '# columns: r theta phi x y z  ' // &
       'vr_mean vr_sig vr_h3 vr_h4  vth_mean vth_sig vth_h3 vth_h4  ' // &
       'vph_mean vph_sig vph_h3 vph_h4'

  ! Draw positions and compute GH params
  do j=1, nsamp

     ! --- radius (same as legacy) ---
     if (abs(gamma - 3.0_dp) > 1.0e-12_dp) then
        rcur = ( (ran1(idum) * (ro**(3.0_dp-gamma) - ri**(3.0_dp-gamma))) + &
                 ri**(3.0_dp-gamma) ) ** (1.0_dp/(3.0_dp-gamma))
     else
        rcur = ri * exp( ran1(idum) * log(ro/ri) )
     end if
     lr = log(rcur)

     ! --- theta rejection sampling (MATCH MCdraw5 exactly) ---
     do
        st = ran1(idum)
        th = asin(st)
        ct = cos(th)

        fmass = ((ct*ct) + ((st/q)*(st/q)))**(-0.5_dp*gamma)
        xi = ran1(idum)
        if (xi <= fmass) exit
     end do

     xxi = ran1(idum)
     if (xxi <= 0.5_dp) th = -1.0_dp * th

     thcur = (pi/2.0_dp) - th

     ! --- phi ---
     phcur = 2.0_dp*pi*ran1(idum)

     sinth = sin(thcur)
     costh = cos(thcur)
     cosph = cos(phcur)
     sinph = sin(phcur)

     xtmp = rcur * sinth * cosph
     ytmp = rcur * sinth * sinph
     ztmp = rcur * costh

     x(j) = xtmp
     y(j) = ytmp
     z(j) = ztmp

     ! Intrinsic moments at theta
     rho   = rhovelmom(thcur,0,0,0)

     vr2   = rhovelmom(thcur,2,0,0)/rho
     vth2  = rhovelmom(thcur,0,2,0)/rho
     vph2  = rhovelmom(thcur,0,0,2)/rho
     vphav = rhovelmom(thcur,0,0,1)/rho

     ! Scale-free scaling with radius
     vr2   = vr2   * exp(-delta*lr)
     vth2  = vth2  * exp(-delta*lr)
     vph2  = vph2  * exp(-delta*lr)
     vphav = vphav * exp(-0.5_dp*delta*lr)

     ! Defaults: Gaussian (h3=h4=0); keep hook for future
     vr_mean  = 0.0_dp
     vr_sig   = sqrt(max(0.0_dp, vr2))
     vr_h3    = 0.0_dp
     vr_h4    = 0.0_dp

     vth_mean = 0.0_dp
     vth_sig  = sqrt(max(0.0_dp, vth2))
     vth_h3   = 0.0_dp
     vth_h4   = 0.0_dp

     vph_mean = vphav
     vph_sig  = sqrt(max(0.0_dp, vph2 - vphav*vphav))
     vph_h3   = 0.0_dp
     vph_h4   = 0.0_dp

     call gh_moments_hook(thcur, lr, &
          vr_mean, vr_sig, vr_h3, vr_h4, &
          vth_mean, vth_sig, vth_h3, vth_h4, &
          vph_mean, vph_sig, vph_h3, vph_h4)

     write(uin,'(18ES25.16)') rcur, thcur, phcur, xtmp, ytmp, ztmp, &
          vr_mean, vr_sig, vr_h3, vr_h4, &
          vth_mean, vth_sig, vth_h3, vth_h4, &
          vph_mean, vph_sig, vph_h3, vph_h4

  end do

  close(uin)

  ! Call python helper (expects gh_velocity_sampler.py in working dir)
  call execute_command_line(trim(pyexe)//' gh_velocity_sampler.py ' // &
       trim(in_path)//' ' // trim(out_path), exitstat=stat)
  if (stat /= 0) then
     write(*,*) 'ERROR: Python helper failed (exit=', stat, ').'
     write(*,*) 'Ensure balrogo is installed and gh_velocity_sampler.py is available.'
     stop 2
  end if

  ! Read vr/vth/vph back
  open(newunit=uout, file=trim(out_path), status='old', action='read')
  call skip_header(uout)
  do j=1, nsamp
     read(uout,*,end=900) vr(j), vth(j), vph(j)
  end do
900 continue
  close(uout)

  ! --- cleanup temporary exchange files (portable, best-effort) ---
  open(newunit=udel, file=trim(in_path), status='old', action='read', iostat=ios_del)
  if (ios_del == 0) then
     close(udel, status='delete')
  else
     write(*,*) '# note: could not open for delete: ', trim(in_path), ' iostat=', ios_del
  end if

  open(newunit=udel, file=trim(out_path), status='old', action='read', iostat=ios_del)
  if (ios_del == 0) then
     close(udel, status='delete')
  else
     write(*,*) '# note: could not open for delete: ', trim(out_path), ' iostat=', ios_del
  end if

  ! Convert to Cartesian velocities
  do j=1, nsamp
     rcur = sqrt(x(j)**2 + y(j)**2 + z(j)**2)
     if (rcur <= 0.0_dp) then
        vx(j)=0.0_dp; vy(j)=0.0_dp; vz(j)=0.0_dp
     else
        thcur = acos(z(j)/rcur)
        phcur = atan2(y(j), x(j))
        costh = cos(thcur); sinth = sin(thcur)
        cosph = cos(phcur); sinph = sin(phcur)

        vx(j) = vr(j)*sinth*cosph + vth(j)*costh*cosph - vph(j)*sinph
        vy(j) = vr(j)*sinth*sinph + vth(j)*costh*sinph + vph(j)*cosph
        vz(j) = vr(j)*costh       - vth(j)*sinth
     end if
  end do

  ! Summary statistics (Binney beta etc.)
  vr2sum = 0.0_dp
  vt2sum = 0.0_dp
  v2sum  = 0.0_dp
  vp1sum = 0.0_dp

  do j=1, nsamp
     vr2sum = vr2sum + vr(j)**2
     vt2sum = vt2sum + vth(j)**2 + vph(j)**2
     v2sum  = v2sum  + vx(j)**2 + vy(j)**2 + vz(j)**2
     vp1sum = vp1sum + vph(j)
  end do

  betabin = 1.0_dp - (vt2sum/(2.0_dp*vr2sum))

  ! Emit results to STDOUT
  write(*,'(A)') '# kind=mock'
  write(*,'(A)') '# columns: x y z vx vy vz'
  do j=1, nsamp
     write(*,'(6ES25.16)') x(j), y(j), z(j), vx(j), vy(j), vz(j)
  end do

  write(*,'(A)') '# kind=mock_stats'
  write(*,'(A)') '# columns: binney_beta rms_method1 rms_method2 mean_vphi'
  write(*,'(4ES25.16)') betabin, sqrt(((vr2sum+vt2sum)/real(nsamp,dp))), &
       sqrt((v2sum/real(nsamp,dp))), vp1sum/real(nsamp,dp)

contains

  subroutine skip_header(u)
    integer, intent(in) :: u
    character(len=1024) :: line
    integer :: ios
    do
       read(u,'(A)',iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       if (line(1:1) /= '#') then
          backspace(u)
          exit
       end if
    end do
  end subroutine skip_header

  subroutine gh_moments_hook(thcur, lr, &
       vr_mean, vr_sig, vr_h3, vr_h4, &
       vth_mean, vth_sig, vth_h3, vth_h4, &
       vph_mean, vph_sig, vph_h3, vph_h4)
    implicit none
    real(dp), intent(in) :: thcur, lr
    real(dp), intent(inout) :: vr_mean, vr_sig, vr_h3, vr_h4
    real(dp), intent(inout) :: vth_mean, vth_sig, vth_h3, vth_h4
    real(dp), intent(inout) :: vph_mean, vph_sig, vph_h3, vph_h4
    ! no-op for now
  end subroutine gh_moments_hook

end program mcdraw_gh
