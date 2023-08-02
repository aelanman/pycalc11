
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
      Integer*2, allocatable, dimension(:, :) :: LNSTAR
      Integer*4, allocatable, dimension(:) :: PhCntr

      Integer*2  NUMSTR, i1dum

      end module srcmod

      MODULE datafiles
      implicit none

      INCLUDE 'param11.i'

      end module datafiles

!      MODULE stations
!      implicit none
!
!      end module stations
