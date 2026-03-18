
! TODO -- Put allocation subroutines in these modules.
      MODULE outputs
      implicit none 

      ! These were initialized with fixed size in difxcalc in c2poly.i
      ! Now, they will be declared allocatable and allocated 
      ! using the subroutine allocoutarrays in dinit.f
      Real*8, allocatable, dimension(:,:,:,:) :: Delay_f, Rate_f,         &
     &    Ubase_f, Vbase_f, Wbase_f
      Real*8, allocatable, dimension(:,:,:,:,:) :: Atmdryd_f, Atmdryr_f,  &
     &    Atmwetd_f, Atmwetr_f, El_f, Az_f, StaX_f, StaY_f, StaZ_f
      Real*8, allocatable, dimension(:,:,:,:,:,:) :: Partials_f
      Integer*4, allocatable, dimension(:,:) :: Iymdhms_f
     
      Character*20 Xsource

      Integer*4 Numsite, Numbaseline, NumPhCenter

      end module outputs

      MODULE srcmod
      implicit none
!    To replace the common block STRCM in cmxsr11.i      
!      Integer*4 MAX_ARC_SRC
!      Parameter(MAX_ARC_SRC=300)

       ! TODO Fixing business with LNSTAR --> This is only used for
       ! mark3 mode

      Real*8, allocatable, dimension(:,:) :: RADEC, P_motion, PRcorr
      Real*8, allocatable, dimension(:) :: D_psec
      Real*8     CD, CRA, SD, SRA
      CHARACTER(len=20), allocatable, dimension(:) :: SrcName
      Integer*4  Pmotion, Dpsec
      Integer(4) NUMSTR
      Integer*2, allocatable, dimension(:, :) :: LNSTAR
      Integer*4, allocatable, dimension(:) :: PhCntr

      Integer*2 i1dum

      end module srcmod

      MODULE datafiles
      implicit none

      INCLUDE 'param11.i'

      end module datafiles

      MODULE ephcom
      implicit none
!
!     Module to hold externally-provided (Python-side) ephemeris data,
!     allowing PEP to bypass the Fortran binary DE file reading.
!
!     Control flag: .true. = use Python-provided ephemeris data
      logical :: use_ext_ephem = .false.
!
!     Maximum number of time steps per 2-minute epoch.
!     Typical value is 6 (for 24s intervals over 2 min).
!     150 covers d_interval as small as ~0.8s.
      Integer*4, parameter :: MAX_EPH_STEPS = 150
!
!     Step counter, reset from Python before each adrivr call.
      Integer*4 :: eph_step_idx = 1
!
!     Pre-computed PEP outputs for each time step.
!     All in meters, m/s, m/s^2, J2000.0 frame.
!
!     Barycentric Earth position, velocity, acceleration
      Real*8 :: ext_earth(3,3,MAX_EPH_STEPS)
!     Geocentric Sun position, velocity
      Real*8 :: ext_sun(3,2,MAX_EPH_STEPS)
!     Geocentric Moon position, velocity
      Real*8 :: ext_xmoon(3,2,MAX_EPH_STEPS)
!     Barycentric planet positions, velocities (7 planets)
      Real*8 :: ext_splanet(3,2,7,MAX_EPH_STEPS)
!     Geocentric planet positions, velocities (7 planets)
      Real*8 :: ext_gplanet(3,2,7,MAX_EPH_STEPS)
!     Barycentric Sun position, velocity
      Real*8 :: ext_sunb(3,2,MAX_EPH_STEPS)
!     Barycentric Moon position, velocity
      Real*8 :: ext_moonb(3,2,MAX_EPH_STEPS)
!
      end module ephcom
