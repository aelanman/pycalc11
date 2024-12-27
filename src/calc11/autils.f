! Originally called dstrt.f
! This file has been cut down to just necessary utility functions.


!**********************************************************************
      REAL*8 FUNCTION JDY2K (IYEAR, IMONTH, IDAY)
      Implicit None
!
!     Function JDY2K: Function to convert year, month, day to full Julian
!     day. The year can be either a four-digit year or a two-digit year.
!
!     If a 4-digit year, this function is good from 1 March 1900 to 31
!     December 2099.
!
!     If a 2-digit year, this function is good from 1 January 1970 to
!     31 December 2069. If year is 70 - 99, 1900 is added. If year is
!     00 - 69, 2000 is added.
!
!     Programmer:
!      98.07.23  D. Gordon  Function written from code in cutcu.f
!
      Integer*4 IYEAR, IMONTH, IDAY, IY, IM, ID
!
       IY = IYEAR
       IM = IMONTH
       ID = IDAY
!
       If (IY .ge. 70 .and. IY .le. 99) Then
        IY = IY + 1900
        Go To 100
       Endif
!
       If (IY .ge. 0 .and. IY .le. 69) Then
        IY = IY + 2000
        Go To 100
       Endif
!
       If (IY .gt.1900 .and. IY .le. 2099) Then
        Go To 100
       Endif
!
!     Year out of range if we get here
       Print *, ' JDY2K, Year out of Range, Stopping! ', IY
       Stop
!
 100   Continue
!
        JDY2K = 367.D0*IY - (7 * ( IY + (IM+9)/12) )/4 + &
     &          (275*IM)/9 + ID + 1721013.5D0
!
!      Write(6,1000) IYEAR, IMONTH, IDAY, IY, IM, ID, JDY2k
 1000  Format(/,'Function JDY2K: ',/,' Input, Modified Y,M,D: ', &
     &        2x,3I5,5x,3I5,/,' JDY2K ', F15.2)
      Return
      End
!
!**********************************************************************
      SUBROUTINE YMDJL(FJD, IYEAR, IMONTH, IDAY)
      Implicit none
!
!     Convert Julian date at midnight to month, day and year.
!     Copied from PEP plus an additional .5 day to convert
!     from Julian date at midnight to the PEP Julian day number.
!
!     INPUT Variables:
      REAL*8 FJD
!
! FJD = Julian date at midnight
!
!     OUTPUT Variables:
      INTEGER*4 IMONTH,IDAY,IYEAR,ITIME
!
! IDAY = Day of the month
! IMONTH = Month of the year (1-12)
! IYEAR = Year of the century (0-99)
!       => Changed to full 4-digit year
! ITIME - Number of centuries since 1900
!
!     LOCAL VARIABLES
      INTEGER*2 MDN(13),nyr,ic,inyr,i
      Real*8 XJD
      DATA MDN/0,31,59,90,120,151,181,212,243,273,304,334,365/
!
! IC - Number of centuries sice 0 January 1600
! INYR - Leap year days since 0 January 1600
! MDN - Array containing day of year at end of each month
! NYR - Number of years since 0 January 1600
! XJD - Days since 0 January 1600
!
!     HISTORY
!   WHO   WHEN     WHAT
!   DSR  7608??    Created
!   DG  2012.04.26 Integer*4 and 4-digit year mods.
!
!     YMDJL PROGRAM STRUCTURE
!
      XJD = FJD - 2305447.0D0 + 0.5D0
      NYR = XJD/365.
   16 IC = NYR/100
!     DAYS DUE TO LEAP YEARS
      INYR = XJD - NYR*365.0
      IDAY = INYR - (NYR-1)/4 + (NYR + 99)/100 - (NYR + 399)/400 - 1
      IF(IC .NE.0) GO TO 20
      IF(NYR.NE.0) GO TO 20
      IDAY = IDAY + 1
   20 IF(IDAY .GT. 0) GO TO 23
      NYR = NYR - 1
      GO TO 16
!**** IYEAR (O THRU 99) YEAR OF THE CENTURY
   23 IYEAR = NYR - IC * 100
      ITIME = IC - 3
      NYR = IYEAR
      IF(NYR .NE. 0) GO TO 27
      IF(MOD(IC,INT2(4)) .NE. 0) GO TO 34
      GO TO 30
   27 IF(MOD(NYR,INT2(4)) .NE. 0) GO TO 34
!
! Replacing ancient code:
!  30 IF(IDAY - 60) 34,39,32
   30 Continue
      If((IDAY-60) .lt. 0) Go to 34
      If((IDAY-60) .eq. 0) Go to 39
      If((IDAY-60) .gt. 0) Go to 32
!  
   32 IDAY = IDAY - 1
   34 DO 36 I=2,13
          IF(IDAY .LE. MDN(I)) GO TO 40
   36 CONTINUE
   39 IMONTH = 2
      IDAY = 29
      GO TO 45
   40 IMONTH = I - 1
      IDAY = IDAY - MDN(IMONTH)
   45 CONTINUE
! Modified - 4 digit year
       IYEAR = IYEAR + 1900 + ITIME*100
!
      RETURN
      END
!
!**********************************************************************
      SUBROUTINE dGet_input(Kjob)
      use srcmod
      IMPLICIT None
!
!
! 3.2.2 COMMON BLOCKS USED -
!
      INCLUDE 'cmxut11.i'
!            Variables 'to':
!              1. Xintv(2)  - First and last Julian Date in the data base
!              1. UT1IF(4)  - The final UT1 information array. This array
!                             contains respectively: 1) The Julian date of the
!                             first tabular point, 2) The increment in days of
!                             the tabular points, 3) The number of tabular
!                             points, 4) The units of the UT1 tabular array per
!                             second. (days, days, unitless, sec/table unit)
!              2. UT1PT(20) - The tabular values of 'TAI minus UT1'.
!                             (table units)
!              3. ISHRTFL   - The short period tidal terms flag, (unitless).
!                             = 1 --> UT1 table coming from input database is
!                             true UT1, (that is, fortnightly tidal terms have
!                             not been removed, as in the IRIS or IERS series).
!                             = -1 --> UT1 table coming from input database is
!                             UT1R, (that is, the Yoder fortnightly tidal terms
!                             HAVE been removed as in Bulletin B).
!                             = -2 --> UT1 table coming from input database is
!                             UT1S, (the S tidal terms HAVE been removed).
!              4. Leap_fix   - Used in external input mode. .True. means
!                              correct the input EOP series for accumluated
!                              leap seconds. .False. means do not correct.
!              5. UT1type    - UT1 data type: 'UT1-TAI ' or 'UT1-UTC '.
!                              For ''UT1-UTC ', leap second corrections
!                              must be made.
!              6. EOP_time_scale - EOP table time scale, allowed values:
!                              'TAI     ', 'TCG     ', 'TDB     ',
!                              'TDT     ', 'UTC     ', 'UNDEF   '.
!                              Assumed default if not present => TDB
!
      INCLUDE 'cmwob11.i'
!            Variables 'to':
!              1. WOBIF(3)  -  The wobble information array. Contains
!                              respectively: 1) The Julian date of the first
!                              tabular point, 2) The increment in days of the
!                              tabular points, 3) The number of tabular points.
!                              (days, days, unitless)
!              2. XYWOB(2,20)- The wobble tabular points for the polar motion
!                              (wobble) X & Y offsets. (milliarcsec)
!                              (Note: Based on old BIH conventions, offsets
!                              are assumed to be left-handed.)
!
      INCLUDE 'cmxst11.i'
!      Variables from:
!       1. Max_Stat             -  Maximum number of stations that can be 
!                                  processed.
!      Variables to:
!*      1. CFRAD(Max_Stat)      -  THE SITE SPHERICAL EARTH RADII.  (M)
!*      2. PLAT(3,Max_Stat)     -  THE PARTIAL DERIVATIVES OF THE SITE CRUST
!*                                 FIXED VECTOR COMPONENTS WITH RESPECT TO THE
!*                                 GEODETIC LATITUDES. (M/RAD)
!*      3. PLON(3,Max_Stat)     -  THE PARTIAL DERIVATIVES OF THE SITE CRUST
!*                                 FIXED VECTOR COMPONENTS WITH RESPECT TO THE
!*                                 EAST LONGITUDES. (M/RAD)
!       4. SITAXO(Max_Stat)     -  THE SITE ANTENNA AXIS OFFSETS. (M)
!       5. SITOAM(11,Max_Stat)  -  THE SITE VERTICAL OCEAN LOADING AMPLITUDES.
!                                  (M)
!       6. SITOPH(11,Max_Stat)  -  THE SITE VERTICAL OCEAN LOADING PHASES.
!                                  (RAD)
!       7. SITXYZ(3,Max_Stat)   -  THE SITE CRUST FIXED X, Y, & Z
!                                  COORDINATES. (M, M, M )
!*      8. SNRM(3,Max_Stat)     -  THE X, Y, AND Z COMPONENTS OF THE SITE
!*                                 NORMAL UNIT VECTORS. (UNITLESS)
!*      9. SITZEN(Max_Stat)     -  THE ZENITH ELECTRICAL PATH LENGTH AT EACH
!*                                 OBSERVATION SITE. (SEC)
!*     10. TCROT(3,3,Max_Stat)  -  THE ROTATION MATRICES WHICH ROTATE THE
!*                                 TOPOCENTRIC REFERENCE SYSTEM TO THE CRUST
!*                                 FIXED REFERENCE SYSTEM FOR EACH SITE.
!*                                 (UNITLESS)
!*     11. XLAT(Max_Stat)       -  THE SITE GEODETIC LATITUDES. (RAD)
!*     12. XLON(Max_Stat)       -  THE SITE EAST LONGITUDES. (RAD)
!      13. KTYPE(Max_Stat)      -  THE SITE ANTENNA AXIS TYPES. (UNITLESS)
!*     14. NLAST(2)             -  THE INTEGER VARIABLE WHICH DETERMINES IF
!*                                 THE BASELINE ID HAS CHANGED FROM ONE
!*                                 OBSERVATION TO THE NEXT.
!*                                 (NOTE: THE SITE GEOMETRY NEED NOT BE
!*                                 RELOADED FOR EACH OBSERVATION IF THE
!*                                 BASELINE ID DOES NOT CHANGE. NLAST MUST BE
!*                                 INITIALIZED TO ZERO IN THE INITIALIZATION
!*                                 SECTION AND PASSED TO THE GEOMETRY SECTION
!*                                 SO THAT IT WILL HAVE ZERO VALUES UNTIL
!*                                 AFTER THE FIRST OBSERVATION IS PROCESSED.)
!*     15. NUMSIT               -  THE NUMBER OF SITES IN THE SITE CATALOG.
!      16. LNSITE(4,Max_Stat)   -  THE EIGHT CHARACTER SITE NAMES OF THE
!                                  SITES IN THE SITE CATALOG. (ALPHAMERIC)
!      17. SITHOA(11,2,Max_Stat) - THE SITE HORIZONTAL OCEAN LOADING
!                                  AMPLITUDES. (M)
!      18. SITHOP(11,2,Max_Stat) - THE SITE HORIZONTAL OCEAN LOADING PHASES.
!                                  (RAD)
!*     19. HEIGHT(Max_Stat)     -  Height above the geoid. (meters)
!*     20. RTROT(3,3,Max_Stat)  -  The rotation matrices which rotate the
!*                                 'radial-transverse' reference system to the
!*                                 crust fixed reference system for each site.
!*                                 (Unitless). The 'radial-transverse' ref.
!*                                 system is nearly identical to the
!*                                 topocentric system. 'Up' is in the radial
!*                                 direction from the center of the Earth;
!*                                 'East' is East; and 'North' is perpendicular
!*                                 to the radial in the north direction.
!*     21. GLAT(Max_Stat)       -  The geocentric latitude at each site. (rad)
!      22. Zero_site            -  The site number of the site at the
!                                  geocenter, if there is one in this data
!                                  set. For correlator usage.
!      23. Dbtilt(Max_Stat,2)   -  Antenna fixed axis tilts, in arc-minutes.
!                                  For alt-az mounts, 1 => East tilt,
!                                  2 => North tilt.
!*     24. Rotilt(3,3,Max_Stat) -  Rotation matrices representing the antenna
!*                                 fixed axis tilts, in the local topocentric
!*                                 station frame. X = Up, Y = East, Z = North.
!      25. OPTL6(6,Max_stat)    -  The site ocean pole tide loading 
!                                  coefficients, interpolated from the Desai
!                                  lat/lon table.
!
!      INCLUDE 'cmxsr11.i'   --- COMMENTED OUT AEL 2021/12/7
!       VARIABLES 'TO':
!         1. LNSTAR(10,MAX_ARC_SRC)- THE ALPHANUMERIC CHARACTER NAMES
!                                    OF THE STARS IN THE STAR CATALOG.
!                                    Now up to 20 characters in difx mode.
!         2. NUMSTR                - THE NUMBER OF STARS IN THE STAR
!                                    CATALOG.
!         3. RADEC(2,MAX_ARC_SRC)  - THE RIGHT ASCENSIONS AND DECLINATIONS
!                                    OF THE STARS IN THE STAR CATALOG.
!                                    (RAD, RAD)
!         4. P_motion(3,Max_arc_src)-The RA and Dec proper motions and
!                                    appropriate epoch for stars in the
!                                    star catalog. (arc-sec/yr, arc-sec/yr,
!                                    year (i.e. 1995.5, etc.))
!         5. D_psec(Max_arc_src)   - Distances, if known, for stars in the
!                                    star catalog. Zero => unknown.
!                                    (parsecs)
!         6. PRcorr(2,Max_arc_src) - Proper motion corrections in RA and
!                                    Declination for each star. (Radians)
!
!
      Character*200 Buf1
!     CHARACTER*128 calc_file_name
      CHARACTER*26 Sou26, Sp26
      CHARACTER*12 Site12, Scan12
      Real*8    JStart, JStop, StartMJD, EopTag(20), TAIUTC(20),        &
     &          UT1UTC(20), ADJUSTL, X_Sp
      Integer*4 AntNum, SrcNum,          NumEOP, EopNum, SpNum
      Integer*4 I, J, Unit1, IOS, IX, IK, KJ, I_Sp, L, Kjob
      Integer*4 get4unit
      Character*6 Mount6 
      Save      Unit1
!
      INCLUDE 'd_input.i'
!     COMMON/Calc_input/ SrcName, ScanID, Sites, Axis, JobID, NumScans, &
!    &       StartYr, StartMo, StartDay, StartHr, StartMin, StartSec,   &
!    &       ScanStrt, ScanDur
!     Character*10 SrcName(100), ScanID
!     Character*8 Sites(Max_stat)
!     Integer*4 JobID, NumScans
!     Integer*4 StartYr, StartMo, StartDay, StartHr, StartMin, StartSec
!     Integer*4 ScanStrt, ScanDur
!     Character*4 Axis(Max_Stat) 
!
!     COMMON/NFO/ SpTag, SpPos, SpVel, NumSpace, NumRows, SpName 
!     Character*10 SpName(10)
!     Real*8 SpTag(100,10), SpPos(100,3,10), SpVel(100,3,10)
!     Integer*4 NumSpace, Numrows(10)
!
!!    Character*8 SrcName(10)
      ! TODO Need to replace these to use allocated arrays
      ! LNSTAR isn't actually used in difx mode.
!      Character*20 SrcName(MAX_ARC_SRC)
!      Equivalence (LNSTAR(1,1), SrcName(1))
!
!  Need to initialize stuff:
!
       NumSit = 0
       NumStr = 0
       NumEOP = 0
       NumSpace = 0
       SpFrame = 'ECJ2'
       Do IK=1, Max_Stat
         SITAXO(IK) = 0.D0
         Do KJ=1,3
          SITXYZ(KJ,IK) = 0.D0
         Enddo
       Enddo
!
       If (Kjob .eq. 1) Unit1 = get4unit()
       OPEN(Unit1,FILE=calc_file_name ,STATUS='OLD', IOSTAT= IOS)
       IF(IOS.ne.0) Write(6,'("Open Error for file: ",A128)') calc_file_name
!
 100  Continue
!
      Read(Unit1,'(A200)',end=200) Buf1
!
      If (Buf1(1:7) .eq.'JOB ID:') Read(Buf1(20:40),*) JobID 
      If (Buf1(1:15).eq.'JOB START TIME:') Read(Buf1(16:40),*) JStart 
      If (Buf1(1:14).eq.'JOB STOP TIME:') Read(Buf1(16:40),*) JStop 
      If (Buf1(1:10).eq.'START MJD:') Read(Buf1(20:40),*) StartMJD
      If (Buf1(1:11).eq.'START YEAR:') Read(Buf1(20:40),*) StartYr
      If (Buf1(1:12).eq.'START MONTH:') Read(Buf1(20:40),*) StartMo
      If (Buf1(1:10).eq.'START DAY:') Read(Buf1(20:40),*) StartDay
      If (Buf1(1:11).eq.'START HOUR:') Read(Buf1(20:40),*) StartHr
      If (Buf1(1:13).eq.'START MINUTE:') Read(Buf1(20:40),*) StartMin
      If (Buf1(1:13).eq.'START SECOND:') Read(Buf1(20:40),*) StartSec
! 
      If (Buf1(1:15).eq.'NUM TELESCOPES:') Read(Buf1(20:40),*) NUMSIT
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv 
      If (Buf1(1:9).eq.'TELESCOPE') Then
!  Don't fill site #1. That will be the Geocenter.

       IX = INDEX(Buf1,'NAME:')
       If (IX .gt. 0) Then
        Read(Buf1(10:IX-1),*) AntNum
        Site12 = Buf1((IX+5):(IX+16)) 
        Site12 = ADJUSTL(Site12)
        Sites(AntNum+2) = Site12(1:8)
       Endif 
!
       IX = INDEX(Buf1,'MOUNT:')
       If (IX .gt. 0) Then
        Read(Buf1(10:IX-1),*) AntNum
        Mount6 = Buf1((IX+6):(IX+11)) 
        Mount6 = ADJUSTL(Mount6)
        Axis(AntNum+2) = Mount6(1:4)
       Endif 
!
       IX = INDEX(Buf1,'OFFSET (m):')
       If (IX .gt. 0) Then
        Read(Buf1(10:IX-1),*) AntNum
        Read(Buf1((IX+11):(IX+22)),*) SITAXO(AntNum+2) 
       Endif 
!
       IX = INDEX(Buf1,'X (m):')
       If (IX .gt. 0) Then
        Read(Buf1(10:IX-1),*) AntNum
        Read(Buf1((IX+7):(IX+26)),*) SITXYZ(1,AntNum+2) 
       Endif 
!
       IX = INDEX(Buf1,'Y (m):')
       If (IX .gt. 0) Then
        Read(Buf1(10:IX-1),*) AntNum
        Read(Buf1((IX+7):(IX+26)),*) SITXYZ(2,AntNum+2) 
       Endif 
!
       IX = INDEX(Buf1,'Z (m):')
       If (IX .gt. 0) Then
        Read(Buf1(10:IX-1),*) AntNum
        Read(Buf1((IX+7):(IX+26)),*) SITXYZ(3,AntNum+2) 
       Endif 
!
      Endif    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!
      If (Buf1(1:12).eq.'NUM SOURCES:') Read(Buf1(13:30),*)  NUMSTR 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv 
      If (Buf1(1:6).eq.'SOURCE') Then
       IX = INDEX(Buf1,'NAME:')
       If (IX .gt. 0) Then
        Read(Buf1(8:IX-1),*) SrcNum
        Sou26 = Buf1((IX+5):(IX+30)) 
        Sou26 = ADJUSTL(Sou26)
        SrcName(SrcNum+1) = Sou26(1:20)
        ! Encode source names in LNSTAR without EQUIVALENCE
        do j=1, 10
           LNSTAR(j,SrcNum+1) = TRANSFER(Sou26((j-1)*2+1:(j-1)*2+2),    &
          & int(0,2))
        enddo
       Endif 
       IX = INDEX(Buf1,'RA:')
       If (IX .gt. 0) Then
        Read(Buf1(7:IX-1),*) SrcNum
        Read(Buf1((IX+3):(IX+33)),*) RADEC(1,SrcNum+1) 
       Endif 
       IX = INDEX(Buf1,'DEC:')
       If (IX .gt. 0) Then
        Read(Buf1(7:IX-1),*) SrcNum
        Read(Buf1((IX+4):(IX+33)),*) RADEC(2,SrcNum+1) 
       Endif 
      Endif 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!
      If (Buf1(1:10).eq.'NUM SCANS:') Read(Buf1(11:30),*) NumScans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Scan information now is read in by subroutine dScan.
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv 
!     If (Buf1(1:4).eq.'SCAN') Then
!      IX = INDEX(Buf1,'IDENTIFIER:')
!      If (IX .gt. 0) Then
!       Read(Buf1(5:IX-1),*) ScanNum
!       Scan12 = Buf1((IX+11):(IX+22)) 
!       Scan12 = ADJUSTL(Scan12)
!       ScanID = Scan12(1:10)
!      Endif 
!      IX = INDEX(Buf1,'START (S):')
!      If (IX .gt. 0) Then
!       Read(Buf1(5:IX-1),*) ScanNum
!       Read(Buf1((IX+10):(IX+17)),*) ScanStrt
!      Endif 
!      IX = INDEX(Buf1,'DUR (S):')
!      If (IX .gt. 0) Then
!       Read(Buf1(5:IX-1),*) ScanNum
!       Read(Buf1((IX+8):(IX+17)),*) ScanDur
!      Endif 
!      IX = INDEX(Buf1,'NUM PHS CTRS:')
!      If (IX .gt. 0) Then
!       Read(Buf1(5:IX-1),*) ScanNum
!       Read(Buf1((IX+13):(IX+16)),*) NumPhCntr
!      Endif 
!     Endif 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!
      If (Buf1(1:9).eq.'NUM EOPS:') Read(Buf1(10:30),*) NumEOP
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv 
      If (Buf1(1:3).eq.'EOP') Then
       IX = INDEX(Buf1,'TIME (mjd):')
       If (IX .gt. 0) Then
        Read(Buf1(5:IX-1),*) EopNum
        Read(Buf1((IX+11):(IX+20)),*) EopTag(EopNum+1) 
       Endif 
       IX = INDEX(Buf1,'TAI_UTC (sec):')
       If (IX .gt. 0) Then
        Read(Buf1(5:IX-1),*) EopNum
        Read(Buf1((IX+14):(IX+20)),*) TAIUTC(EopNum+1) 
       Endif 
       IX = INDEX(Buf1,'UT1_UTC (sec):')
       If (IX .gt. 0) Then
        Read(Buf1(5:IX-1),*) EopNum
        Read(Buf1((IX+14):(IX+25)),*) UT1UTC(EopNum+1) 
       Endif 
       IX = INDEX(Buf1,'XPOLE (arcsec):')
       If (IX .gt. 0) Then
        Read(Buf1(5:IX-1),*) EopNum
        Read(Buf1((IX+15):(IX+25)),*) XYWOB(1,EopNum+1) 
       Endif 
       IX = INDEX(Buf1,'YPOLE (arcsec):')
       If (IX .gt. 0) Then
        Read(Buf1(5:IX-1),*) EopNum
        Read(Buf1((IX+15):(IX+25)),*) XYWOB(2,EopNum+1) 
       Endif 
      Endif 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv 
      If (Buf1(1:15).eq.'NUM SPACECRAFT:') Read(Buf1(16:30),*) NumSpace
       If (Buf1(1:6).eq.'FRAME:') Read(Buf1(21:24),*) SpFrame 
      If (Buf1(1:10).eq.'SPACECRAFT') Then
       IX = INDEX(Buf1,'NAME:')
       If (IX .gt. 0) Then
        Read(Buf1(11:IX-1),*) SpNum
!       Read(Buf1((IX+5):(IX+30)),*) SpName(SpNum+1) 
        Sp26 = Buf1((IX+5):(IX+30)) 
        Sp26 = ADJUSTL(Sp26)
        SpName(SpNum+1) = Sp26(1:20)
!     write(6,*) 'NumSpace,SpNum,SpName(SpNum+1))', NumSpace,SpNum,SpName(SpNum+1)
       Endif 
       IX = INDEX(Buf1,'ROWS:')
       If (IX .gt. 0) Then
        Read(Buf1(11:IX-1),*) SpNum 
        Read(Buf1((IX+5):(IX+14)),*) Numrows(SpNum+1) 
!  Check for too many spacecraft positions.
           If (Numrows(SpNum+1) .gt. NF_row) Then
            Write(6,1031) Numrows(SpNum+1), NF_row
 1031       Format(/,' Number of spacecraft position rows = ',I4,/,     &
     &            ' The current limit is ',I4,/,' Increase the value ', &
     &            'of NF_row in d_input.i ',/)
            Stop 
           Endif
          Do I = 1, Numrows(SpNum+1)
            Read(Unit1,'(A200)',end=200) Buf1
             If (Buf1(1:6).eq.'FRAME:') Then
               Read(Buf1(21:24),*) SpFrame
               Go to 190
             Endif
            IX = INDEX(Buf1,':')
            Read(Buf1(IX+1:200),*) SpTag(I,SpNum+1),SpPos(I,1,SpNum+1), &
     &       SpPos(I,2,SpNum+1),SpPos(I,3,SpNum+1), SpVel(I,1,SpNum+1), &
     &       SpVel(I,2,SpNum+1),SpVel(I,3,SpNum+1)
!             Shift time tags by t_offset if requested.
            If (SpOffset .eq. 'Offset  ')                               &
     &             SpTag(I,SpNum+1) = SpTag(I,SpNum+1) + t_offset
 190         Continue
          Enddo
       Endif 
      Endif 
!
!
       Go to 100
!
  200 Continue
!
      Close (Unit1)
!
! Setup WOBIF and UT1IF arrays:
      WOBIF(1) = EOPTag(1) + 2400000.5D0
!     WOBIF(2) = 1.0D0 
      WOBIF(2) = EOPTag(2) - EOPTag(1)
      WOBIF(3) = NumEop
!
      UT1IF(1) = EOPTag(1) + 2400000.5D0
!     UT1IF(2) = 1.0D0
      UT1IF(2) = EOPTag(2) - EOPTag(1)
      UT1IF(3) = NumEop
      UT1IF(4) = 1.0D0 
!
       Xintv(1) = JStart + 2400000.5D0
       Xintv(2) = JStop  + 2400000.5D0
!
!  Save leap seconds
       Xleap_sec = TAIUTC(1)
!
!  If multiple scans, get the total duration time.
!!!!   If(NumScans .gt. 1) ScanDur = ScanStrt + ScanDur
!
! Change XY-wobble to milli-arc-seconds
      Do J = 1, NumEOP
       UT1PT(J) = TAIUTC(J) - UT1UTC(J)
       XYWOB(1,J) = XYWOB(1,J) * 1.D3
       XYWOB(2,J) = XYWOB(2,J) * 1.D3
      Enddo
!
      Numsrc = NUMSTR
      Return
      End
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FixEpoch(JTAG, TAG_SEC)
      IMPLICIT None
!
      Real*8    TAG_SEC
      Integer*4 JTAG(5), Imin
      Integer IMNTHS(12)
      DATA IMNTHS /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
!      Write (6,*) 'FixEpoch: ', JTAG, TAG_SEC
!
! Check for leap year
      IF(MOD(JTAG(1),4) .EQ. 0) IMNTHS(2) = 29
!
      IF (TAG_SEC .ge. 59.9999999) Then  ! Increment minutes
       Imin = TAG_SEC/60.D0 + .00001
       TAG_SEC = TAG_SEC - Imin*60.0D0
       JTAG(5) = JTAG(5) + Imin  ! minutes
        If (JTAG(5) .ge. 60) Then     ! Increment hours
         JTAG(5) = JTAG(5) - 60
         JTAG(4) = JTAG(4) + 1  ! hours
          If (JTAG(4) .ge. 24) Then     ! Increment days
           JTAG(4) = JTAG(4) - 24
           JTAG(3) = JTAG(3) + 1   ! days
            If (JTAG(3) .gt. IMNTHS(JTAG(2))) Then ! Increment month
             JTAG(3) = JTAG(3) - IMNTHS(JTAG(2))    ! days, should be 1
             JTAG(2) = JTAG(2) + 1   ! month
             If (JTAG(2) .gt. 12) Then               ! Increment year
              JTAG(2) = 1
              JTAG(1) = JTAG(1) + 1
             Else
              Return
             Endif                                   ! Increment year
            Else
             Return
            Endif                                  ! Increment month
          Else
           Return
          Endif                         ! Increment days
        Else
         Return
        Endif                         ! Increment hours
      ENDIF                              ! Increment minutes
!      Write (6,*) 'FixEpoch: ', JTAG, TAG_SEC
!
      Return
      End
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FixEpoch2(JTAG, TAG_SEC)
      IMPLICIT None
!
      Real*8    TAG_SEC, Xmin, Xhr
      Integer*4 JTAG(5), Imin, Ihr
      Integer IMNTHS(12)
      DATA IMNTHS /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
!      Write (6,*) 'FixEpoch2/Before: ', JTAG, TAG_SEC
!
! Check for leap year
      IF(MOD(JTAG(1),4) .EQ. 0) IMNTHS(2) = 29
!
      IF (TAG_SEC .ge. 59.9999999) Then  ! Reset seconds and increment minutes
       Xmin = TAG_SEC/60.D0 + .00001
       Imin = Xmin
       TAG_SEC = TAG_SEC - Imin*60.0D0
       JTAG(5) = JTAG(5) + Imin  ! minutes
      ENDIF
!
       If (JTAG(5) .ge. 60) Then     ! Reset minutes ana increment hours
        Xmin = JTAG(5)
        Xhr  = Xmin/60.D0
        Ihr  = Xhr
        JTAG(5) = JTAG(5) - Ihr * 60  ! minutes
        JTAG(4) = JTAG(4) + Ihr 
       Endif
!
          If (JTAG(4) .ge. 24) Then     ! Reset hours and increment days
           JTAG(4) = JTAG(4) - 24
           JTAG(3) = JTAG(3) + 1   ! days
            If (JTAG(3) .gt. IMNTHS(JTAG(2))) Then ! Increment month
             JTAG(3) = JTAG(3) - IMNTHS(JTAG(2))    ! days, should be 1
             JTAG(2) = JTAG(2) + 1   ! month
             If (JTAG(2) .gt. 12) Then               ! Increment year
              JTAG(2) = 1
              JTAG(1) = JTAG(1) + 1
             Else
              Return
             Endif                                   ! Increment year
            Else
             Return
            Endif                                  ! Increment month
          Else
           Return
          Endif                         ! Increment days
!
!      Write (6,*) 'FixEpoch2/After:  ', JTAG, TAG_SEC
!
      Return
      End
