module bessel
  use misc
  !
  ! This module is taken from MCJYNA.F90 in zhangjin.zip (http://jblevins.org/mirror/amiller/) and slightly modified.
  ! See `Computation of Special Functions' by Zhang & Jin.
  !
  implicit none

  interface besh0_1
     module procedure dbesh0_1
     module procedure cbesh0_1
  end interface besh0_1

  interface besj
     module procedure besj_real
     module procedure besj_complex
  end interface besj

  interface besy
     module procedure besy_real
     module procedure besy_complex
  end interface besy

  interface besh
     module procedure besh_real
     module procedure besh_complex
  end interface besh

  interface besh_2
     module procedure besh_2_real
     module procedure besh_2_complex
  end interface besh_2

  interface besj_d
     module procedure besj_d_real
     module procedure besj_d_complex
  end interface besj_d

  interface besy_d
     module procedure besy_d_real
     module procedure besy_d_complex
  end interface besy_d

  interface besh_d
     module procedure besh_d_real
     module procedure besh_d_complex
  end interface besh_d

  interface besh_2_d
     module procedure besh_2_d_real
     module procedure besh_2_d_complex
  end interface besh_2_d

  interface
     ! gsl/gsl_sf_bessel.h
     function gsl_sf_bessel_K0(x) result(out) bind(c,name="gsl_sf_bessel_K0")
       use iso_c_binding
       real(c_double),value :: x
       real(c_double) :: out
     end function gsl_sf_bessel_K0

     function gsl_sf_bessel_K1(x) result(out) bind(c,name="gsl_sf_bessel_K1")
       use iso_c_binding
       real(c_double),value :: x
       real(c_double) :: out
     end function gsl_sf_bessel_K1

     function gsl_sf_bessel_Kn_array(nmin, nmax, x, result_array) result(out) bind(c,name="gsl_sf_bessel_Kn_array")
       use iso_c_binding
       integer(c_int),value :: nmin
       integer(c_int),value :: nmax
       real(c_double),value :: x
       type(c_ptr),value :: result_array
       integer(c_int) :: out
     end function gsl_sf_bessel_Kn_array
  end interface

  INTEGER, PARAMETER, private  :: dp = SELECTED_REAL_KIND(12, 60)

  REAL (dp), PARAMETER, private  :: a_1(12) = (/ -.703125D-01, .112152099609375D+00, -  &
       .5725014209747314D+00, .6074042001273483D+01, -  &
       .1100171402692467D+03, .3038090510922384D+04, -  &
       .1188384262567832D+06, .6252951493434797D+07, -  &
       .4259392165047669D+09, .3646840080706556D+11, -  &
       .3833534661393944D+13, .4854014686852901D+15 /)
  REAL (dp), PARAMETER, private  :: b_1(12) = (/ .732421875D-01, -.2271080017089844D+00,  &
       .1727727502584457D+01, -.2438052969955606D+02,  &
       .5513358961220206D+03, -.1825775547429318D+05,  &
       .8328593040162893D+06, -.5006958953198893D+08,  &
       .3836255180230433D+10, -.3649010818849833D+12,  &
       .4218971570284096D+14, -.5827244631566907D+16 /)
  REAL (dp), PARAMETER, private  :: a1_1(12) = (/ .1171875D+00, -.144195556640625D+00,  &
       .6765925884246826D+00, -.6883914268109947D+01,  &
       .1215978918765359D+03, -.3302272294480852D+04,  &
       .1276412726461746D+06, -.6656367718817688D+07,  &
       .4502786003050393D+09, -.3833857520742790D+11,  &
       .4011838599133198D+13, -.5060568503314727D+15 /)
  REAL (dp), PARAMETER, private  :: b1_1(12) = (/ -.1025390625D+00, .2775764465332031D+00, -  &
       .1993531733751297D+01, .2724882731126854D+02, -  &
       .6038440767050702D+03, .1971837591223663D+05, -  &
       .8902978767070678D+06, .5310411010968522D+08, -  &
       .4043620325107754D+10, .3827011346598605D+12, -  &
       .4406481417852278D+14, .6065091351222699D+16 /)
  
  REAL (dp), PARAMETER, private  :: a_2(12) = (/ 0.125D0, 7.03125D-2, 7.32421875D-2,  &
       1.1215209960938D-1, 2.2710800170898D-1, 5.7250142097473D-1,   &
       1.7277275025845D0, 6.0740420012735D0, 2.4380529699556D01,   &
       1.1001714026925D02, 5.5133589612202D02, 3.0380905109224D03 /)
  REAL (dp), PARAMETER, private  :: b_2(12) = (/ -0.375D0, -1.171875D-1, -1.025390625D-1,  &
       -1.4419555664063D-1, -2.7757644653320D-1, -6.7659258842468D-1,  &
       -1.9935317337513D0, -6.8839142681099D0, -2.7248827311269D01,  &
       -1.2159789187654D02, -6.0384407670507D02, -3.3022722944809D03 /)
  REAL (dp), PARAMETER, private  :: a1_2(8) = (/ 0.125D0, 0.2109375D0, 1.0986328125D0,  &
       1.1775970458984D01, 2.1461706161499D02, 5.9511522710323D03,  &
       2.3347645606175D05, 1.2312234987631D07 /)
  REAL (dp), PARAMETER, private  :: pi = 3.141592653589793_dp, el = 0.5772156649015329_dp

  real(dp),parameter, private :: rp2 = 2.0_dp / pi
  complex(dp),parameter, private :: ci = (0.0_dp,1.0_dp)
contains

  SUBROUTINE cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)

    !       ===========================================================
    !       Purpose: Compute complex Bessel functions J0(z), J1(z)
    !                Y0(z), Y1(z), and their derivatives
    !       Input :  z --- Complex argument
    !       Output:  CBJ0 --- J0(z)
    !                CDJ0 --- J0'(z)
    !                CBJ1 --- J1(z)
    !                CDJ1 --- J1'(z)
    !                CBY0 --- Y0(z)
    !                CDY0 --- Y0'(z)
    !                CBY1 --- Y1(z)
    !                CDY1 --- Y1'(z)
    !       ===========================================================

    COMPLEX (dp), INTENT(IN)   :: z
    COMPLEX (dp), INTENT(OUT)  :: cbj0
    COMPLEX (dp), INTENT(OUT)  :: cdj0
    COMPLEX (dp), INTENT(OUT)  :: cbj1
    COMPLEX (dp), INTENT(OUT)  :: cdj1
    COMPLEX (dp), INTENT(OUT)  :: cby0
    COMPLEX (dp), INTENT(OUT)  :: cdy0
    COMPLEX (dp), INTENT(OUT)  :: cby1
    COMPLEX (dp), INTENT(OUT)  :: cdy1

    
    COMPLEX (dp)  :: cp, cp0, cp1, cq0, cq1, cr, cs, ct1, ct2, cu, z1, z2
    REAL (dp)     :: a0, w0, w1
    INTEGER       :: k, k0
    
    a0 = ABS(z)
    z2 = z * z
    z1 = z
    IF (a0 == 0.0_dp) THEN
       cbj0 = (1.0_dp,0.0_dp)
       cbj1 = (0.0_dp,0.0_dp)
       cdj0 = (0.0_dp,0.0_dp)
       cdj1 = (0.5_dp,0.0_dp)
       cby0 = -(1.0D300,0.0_dp)
       cby1 = -(1.0D300,0.0_dp)
       cdy0 = (1.0D300,0.0_dp)
       cdy1 = (1.0D300,0.0_dp)
       RETURN
    END IF
    IF (REAL(z) < 0.0) z1 = -z
    IF (a0 <= 12.0) THEN
       cbj0 = (1.0_dp,0.0_dp)
       cr = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          cr = -0.25_dp * cr * z2 / (k*k)
          cbj0 = cbj0 + cr
          IF (ABS(cr/cbj0) < 1.0D-15) EXIT
       END DO

       cbj1 = (1.0_dp,0.0_dp)
       cr = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          cr = -0.25_dp * cr * z2 / (k*(k+1))
          cbj1 = cbj1 + cr
          IF (ABS(cr/cbj1) < 1.0D-15) EXIT
       END DO

       cbj1 = 0.5_dp * z1 * cbj1
       w0 = 0.0_dp
       cr = (1.0_dp,0.0_dp)
       cs = (0.0_dp,0.0_dp)
       DO  k = 1, 40
          w0 = w0 + 1.0_dp / k
          cr = -0.25_dp * cr / (k*k) * z2
          cp = cr * w0
          cs = cs + cp
          IF (ABS(cp/cs) < 1.0D-15) EXIT
       END DO
       cby0 = rp2 * (LOG(z1/2.0_dp) + el) * cbj0 - rp2 * cs
       w1 = 0.0_dp
       cr = (1.0_dp,0.0_dp)
       cs = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          w1 = w1 + 1.0_dp / k
          cr = -0.25_dp * cr / (k*(k+1)) * z2
          cp = cr * (2.0_dp*w1 + 1.0_dp/(k+1))
          cs = cs + cp
          IF (ABS(cp/cs) < 1.0D-15) EXIT
       END DO
       cby1 = rp2 * ((LOG(z1/2.0_dp) + el)*cbj1 - 1.0_dp/z1 - .25_dp*z1*cs)
    ELSE
       k0 = 12
       IF (a0 >= 35.0) k0 = 10
       IF (a0 >= 50.0) k0 = 8
       ct1 = z1 - 0.25_dp * pi
       cp0 = (1.0_dp,0.0_dp)
       DO  k = 1, k0
          cp0 = cp0 + a_1(k) * z1 ** (-2*k)
       END DO
       cq0 = -0.125_dp / z1
       DO  k = 1, k0
          cq0 = cq0 + b_1(k) * z1 ** (-2*k-1)
       END DO
       cu = SQRT(rp2/z1)
       cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
       cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
       ct2 = z1 - 0.75_dp * pi
       cp1 = (1.0_dp,0.0_dp)
       DO  k = 1, k0
          cp1 = cp1 + a1_1(k) * z1 ** (-2*k)
       END DO
       cq1 = 0.375_dp / z1
       DO  k = 1, k0
          cq1 = cq1 + b1_1(k) * z1 ** (-2*k-1)
       END DO
       cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
       cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
    END IF
    IF (REAL(z) < 0.0) THEN
       IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0_dp * ci * cbj0
       IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0_dp * ci * cbj0
       IF (AIMAG(z) < 0.0) cby1 = -(cby1 - 2.0_dp*ci*cbj1)
       IF (AIMAG(z) > 0.0) cby1 = -(cby1 + 2.0_dp*ci*cbj1)
       cbj1 = -cbj1
    END IF
    cdj0 = -cbj1
    cdj1 = cbj0 - 1.0_dp / z * cbj1
    cdy0 = -cby1
    cdy1 = cby0 - 1.0_dp / z * cby1
    RETURN
  END SUBROUTINE cjy01

  SUBROUTINE cjy02(z, cbj0, cbj1, cby0, cby1)

    !       ===========================================================
    !       Purpose: Compute complex Bessel functions J0(z), J1(z)
    !                Y0(z), Y1(z)
    !       Input :  z --- Complex argument
    !       Output:  CBJ0 --- J0(z)
    !                CBJ1 --- J1(z)
    !                CBY0 --- Y0(z)
    !                CBY1 --- Y1(z)
    !       ===========================================================

    COMPLEX (dp), INTENT(IN)   :: z
    COMPLEX (dp), INTENT(OUT)  :: cbj0
    COMPLEX (dp), INTENT(OUT)  :: cbj1
    COMPLEX (dp), INTENT(OUT)  :: cby0
    COMPLEX (dp), INTENT(OUT)  :: cby1

    
    COMPLEX (dp)  :: cp, cp0, cp1, cq0, cq1, cr, cs, ct1, ct2, cu, z1, z2
    REAL (dp)     :: a0, w0, w1
    INTEGER       :: k, k0
    
    a0 = ABS(z)
    z2 = z * z
    z1 = z
    IF (a0 == 0.0_dp) THEN
       cbj0 = (1.0_dp,0.0_dp)
       cbj1 = (0.0_dp,0.0_dp)
       cby0 = -(1.0D300,0.0_dp)
       cby1 = -(1.0D300,0.0_dp)
       RETURN
    END IF
    IF (REAL(z) < 0.0) z1 = -z
    IF (a0 <= 12.0) THEN
       cbj0 = (1.0_dp,0.0_dp)
       cr = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          cr = -0.25_dp * cr * z2 / (k*k)
          cbj0 = cbj0 + cr
          IF (ABS(cr/cbj0) < 1.0D-15) EXIT
       END DO

       cbj1 = (1.0_dp,0.0_dp)
       cr = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          cr = -0.25_dp * cr * z2 / (k*(k+1))
          cbj1 = cbj1 + cr
          IF (ABS(cr/cbj1) < 1.0D-15) EXIT
       END DO
       cbj1 = 0.5_dp * z1 * cbj1
       
       w0 = 0.0_dp
       cr = (1.0_dp,0.0_dp)
       cs = (0.0_dp,0.0_dp)
       DO  k = 1, 40
          w0 = w0 + 1.0_dp / k
          cr = -0.25_dp * cr / (k*k) * z2
          cp = cr * w0
          cs = cs + cp
          IF (ABS(cp/cs) < 1.0D-15) EXIT
       END DO
       cby0 = rp2 * (LOG(z1/2.0_dp) + el) * cbj0 - rp2 * cs
       
       w1 = 0.0_dp
       cr = (1.0_dp,0.0_dp)
       cs = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          w1 = w1 + 1.0_dp / k
          cr = -0.25_dp * cr / (k*(k+1)) * z2
          cp = cr * (2.0_dp*w1 + 1.0_dp/(k+1))
          cs = cs + cp
          IF (ABS(cp/cs) < 1.0D-15) EXIT
       END DO
       cby1 = rp2 * ((LOG(z1/2.0_dp) + el)*cbj1 - 1.0_dp/z1 - .25_dp*z1*cs)
    ELSE
       k0 = 12
       IF (a0 >= 35.0) k0 = 10
       IF (a0 >= 50.0) k0 = 8
       ct1 = z1 - 0.25_dp * pi
       cp0 = (1.0_dp,0.0_dp)
       DO  k = 1, k0
          cp0 = cp0 + a_1(k) * z1 ** (-2*k)
       END DO
       cq0 = -0.125_dp / z1
       DO  k = 1, k0
          cq0 = cq0 + b_1(k) * z1 ** (-2*k-1)
       END DO
       cu = SQRT(rp2/z1)
       cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
       cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
       ct2 = z1 - 0.75_dp * pi
       cp1 = (1.0_dp,0.0_dp)
       DO  k = 1, k0
          cp1 = cp1 + a1_1(k) * z1 ** (-2*k)
       END DO
       cq1 = 0.375_dp / z1
       DO  k = 1, k0
          cq1 = cq1 + b1_1(k) * z1 ** (-2*k-1)
       END DO
       cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
       cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
    END IF
    IF (REAL(z) < 0.0) THEN
       IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0_dp * ci * cbj0
       IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0_dp * ci * cbj0
       IF (AIMAG(z) < 0.0) cby1 = -(cby1 - 2.0_dp*ci*cbj1)
       IF (AIMAG(z) > 0.0) cby1 = -(cby1 + 2.0_dp*ci*cbj1)
       cbj1 = -cbj1
    END IF
    RETURN
  END SUBROUTINE cjy02

  ! from MIK01A.F90
  ! ABS((bk0-ww)/bk0) < 1.0D-15 でwwが初期化されていないbug
  ! k > 1 .and. ABS((bk0-ww)/bk0) < 1.0D-15 に変更
  SUBROUTINE ik01a(x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)

    !       =========================================================
    !       Purpose: Compute modified Bessel functions I0(x), I1(1),
    !                K0(x) and K1(x), and their derivatives
    !       Input :  x   --- Argument ( x ò 0 )
    !       Output:  BI0 --- I0(x)
    !                DI0 --- I0'(x)
    !                BI1 --- I1(x)
    !                DI1 --- I1'(x)
    !                BK0 --- K0(x)
    !                DK0 --- K0'(x)
    !                BK1 --- K1(x)
    !                DK1 --- K1'(x)
    !       =========================================================

    REAL (dp), INTENT(IN)   :: x
    REAL (dp), INTENT(OUT)  :: bi0
    REAL (dp), INTENT(OUT)  :: di0
    REAL (dp), INTENT(OUT)  :: bi1
    REAL (dp), INTENT(OUT)  :: di1
    REAL (dp), INTENT(OUT)  :: bk0
    REAL (dp), INTENT(OUT)  :: dk0
    REAL (dp), INTENT(OUT)  :: bk1
    REAL (dp), INTENT(OUT)  :: dk1
    
    REAL (dp)  :: ca, cb, ct, r, w0, ww, x2, xr, xr2
    INTEGER    :: k, k0

    x2 = x * x
    IF (x == 0.0D0) THEN
       bi0 = 1.0D0
       bi1 = 0.0D0
       bk0 = 1.0D+300
       bk1 = 1.0D+300
       di0 = 0.0D0
       di1 = 0.5D0
       dk0 = -1.0D+300
       dk1 = -1.0D+300
       RETURN
    ELSE IF (x <= 18.0D0) THEN
       bi0 = 1.0D0
       r = 1.0D0
       DO  k = 1, 50
          r = 0.25D0 * r * x2 / (k*k)
          bi0 = bi0 + r
          IF (ABS(r/bi0) < 1.0D-15) EXIT
       END DO
       bi1 = 1.0D0
       r = 1.0D0
       DO  k = 1, 50
          r = 0.25D0 * r * x2 / (k*(k+1))
          bi1 = bi1 + r
          IF (ABS(r/bi1) < 1.0D-15) EXIT
       END DO
       bi1 = 0.5D0 * x * bi1
    ELSE
       k0 = 12
       IF (x >= 35.0) k0 = 9
       IF (x >= 50.0) k0 = 7
       ca = EXP(x) / SQRT(2.0D0*pi*x)
       bi0 = 1.0D0
       xr = 1.0D0 / x
       DO  k = 1, k0
          bi0 = bi0 + a_2(k) * xr ** k
       END DO
       bi0 = ca * bi0
       bi1 = 1.0D0
       DO  k = 1, k0
          bi1 = bi1 + b_2(k) * xr ** k
       END DO
       bi1 = ca * bi1
    END IF
    IF (x <= 9.0D0) THEN
       ct = -(LOG(x/2.0D0) + el)
       bk0 = 0.0D0
       w0 = 0.0D0
       r = 1.0D0
       DO  k = 1, 50
          w0 = w0 + 1.0D0 / k
          r = 0.25D0 * r / (k*k) * x2
          bk0 = bk0 + r * (w0+ct)
          IF (k > 1 .and. ABS((bk0-ww)/bk0) < 1.0D-15) EXIT
          ww = bk0
       END DO
       bk0 = bk0 + ct
    ELSE
       cb = 0.5D0 / x
       xr2 = 1.0D0 / x2
       bk0 = 1.0D0
       DO  k = 1, 8
          bk0 = bk0 + a1_2(k) * xr2 ** k
       END DO
       bk0 = cb * bk0 / bi0
    END IF
    bk1 = (1.0D0/x-bi1*bk0) / bi0
    di0 = bi1
    di1 = bi0 - bi1 / x
    dk0 = -bk1
    dk1 = -bk0 - bk1 / x
    RETURN
  END SUBROUTINE ik01a
  
  ! 0次のHankel関数
  complex(dp) function cbesh0(z)
    complex(dp),intent(in) :: z
    complex(dp) :: cbj0, cby0
    complex(dp) :: cp, cp0, cq0, cr, cs, ct1, cu, z1, z2
    real(dp) :: a0, w0
    integer :: k, k0

    a0 = ABS(z)
    z2 = z * z
    z1 = z
    IF (a0 == 0.0_dp) THEN
       cbj0 = (1.0_dp,0.0_dp)
       RETURN
    END IF
    IF (REAL(z) < 0.0) z1 = -z
    IF (a0 <= 12.0) THEN
       cbj0 = (1.0_dp,0.0_dp)
       cr = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          cr = -0.25_dp * cr * z2 / (k*k)
          cbj0 = cbj0 + cr
          IF (ABS(cr/cbj0) < 1.0D-15) EXIT
       END DO

       w0 = 0.0_dp
       cr = (1.0_dp,0.0_dp)
       cs = (0.0_dp,0.0_dp)
       DO  k = 1, 40
          w0 = w0 + 1.0_dp / k
          cr = -0.25_dp * cr / (k*k) * z2
          cp = cr * w0
          cs = cs + cp
          IF (ABS(cp/cs) < 1.0D-15) EXIT
       END DO
       cby0 = rp2 * (LOG(z1/2.0_dp) + el) * cbj0 - rp2 * cs

    else
       k0 = 12
       IF (a0 >= 35.0) k0 = 10
       IF (a0 >= 50.0) k0 = 8
       ct1 = z1 - 0.25_dp * pi
       cp0 = (1.0_dp,0.0_dp)
       DO  k = 1, k0
          cp0 = cp0 + a_1(k) * z1 ** (-2*k)
       END DO
       cq0 = -0.125_dp / z1
       DO  k = 1, k0
          cq0 = cq0 + b_1(k) * z1 ** (-2*k-1)
       END DO
       cu = SQRT(rp2/z1)
       cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
       cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
    end IF

    IF (REAL(z) < 0.0) THEN
       IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0_dp * ci * cbj0
       IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0_dp * ci * cbj0
    END IF
    
    cbesh0 = cbj0+ci*cby0

  end function cbesh0

  complex(dp) function cbesh1(z)
    complex(dp),intent(in) :: z
    complex(dp) :: cbj1, cby1
    COMPLEX (dp)  :: cp, cp1, cq1, cr, cs, ct2, cu, z1, z2
    REAL (dp)     :: a0, w1
    INTEGER       :: k, k0

    a0 = ABS(z)
    z2 = z * z
    z1 = z
    IF (a0 == 0.0_dp) THEN
       cbj1 = (0.0_dp,0.0_dp)
       cby1 = -(1.0D300,0.0_dp)
       RETURN
    END IF
    IF (REAL(z) < 0.0) z1 = -z
    IF (a0 <= 12.0) THEN
       cbj1 = (1.0_dp,0.0_dp)
       cr = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          cr = -0.25_dp * cr * z2 / (k*(k+1))
          cbj1 = cbj1 + cr
          IF (ABS(cr/cbj1) < 1.0D-15) EXIT
       END DO
       cbj1 = 0.5_dp * z1 * cbj1

       w1 = 0.0_dp
       cr = (1.0_dp,0.0_dp)
       cs = (1.0_dp,0.0_dp)
       DO  k = 1, 40
          w1 = w1 + 1.0_dp / k
          cr = -0.25_dp * cr / (k*(k+1)) * z2
          cp = cr * (2.0_dp*w1 + 1.0_dp/(k+1))
          cs = cs + cp
          IF (ABS(cp/cs) < 1.0D-15) EXIT
       END DO
       cby1 = rp2 * ((LOG(z1/2.0_dp) + el)*cbj1 - 1.0_dp/z1 - .25_dp*z1*cs)

    else
       k0 = 12
       IF (a0 >= 35.0) k0 = 10
       IF (a0 >= 50.0) k0 = 8
       cu = SQRT(rp2/z1)
       
       ct2 = z1 - 0.75_dp * pi
       cp1 = (1.0_dp,0.0_dp)
       DO  k = 1, k0
          cp1 = cp1 + a1_1(k) * z1 ** (-2*k)
       END DO
       cq1 = 0.375_dp / z1
       DO  k = 1, k0
          cq1 = cq1 + b1_1(k) * z1 ** (-2*k-1)
       END DO
       cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
       cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
    end IF

    IF (REAL(z) < 0.0) THEN
       IF (AIMAG(z) < 0.0) cby1 = -(cby1 - 2.0_dp*ci*cbj1)
       IF (AIMAG(z) > 0.0) cby1 = -(cby1 + 2.0_dp*ci*cbj1)
       cbj1 = -cbj1
    END IF

    cbesh1 = cbj1+ci*cby1
    
  end function cbesh1

  ! Hankel function of first kind and order 0 and 1
  ! cbesh0 + cbesh1 よりも若干はやい?
  subroutine cbesh(z, h0, h1)
    complex(dp),intent(in) :: z
    complex(dp),intent(out) :: h0, h1
    complex(dp)  :: cbj0, cbj1, cby0, cby1

    call cjy02(z, cbj0, cbj1, cby0, cby1)

    h0 = cbj0+ci*cby0
    h1 = cbj1+ci*cby1

    return
    
  end subroutine cbesh

  ! n1~n2次のbessel関数を計算
  subroutine cbesjn(z, n1, n2, out)    
    complex(8),intent(in) :: z
    integer,intent(in) :: n1, n2
    complex(8),intent(out) :: out(n1:n2)

    real(8) :: cyr(n1:n2), cyi(n1:n2)
    integer :: nz, ierr
    integer :: n

    if(n1 >= 0) then  
       call zbesj(real(z), aimag(z), n1*1.d0, 1, n2-n1+1, cyr, cyi, nz, ierr)
       if(ierr /= 0) then
          write(*,*) "Error: ierr=", ierr
          call assert(.false.)
          ! stop
       end if
       
    else
       if(n2 < 0) then          
          ! n1 ~ n2次は -n2 ~ -n1次から計算
          call zbesj(real(z), aimag(z), -n2*1.d0, 1, -n1+n2+1, cyr, cyi, nz, ierr)
          if(ierr /= 0) then
             write(*,*) "Error: ierr=", ierr
             call assert(.false.)
             ! stop
          end if
          ! 順番を反転
          cyr = cyr(n2:n1:-1)
          cyi = cyi(n2:n1:-1)
          ! (-1)^nを掛ける
          do n=n1,n2
             cyr(n) = cyr(n)*(-1)**n
             cyi(n) = cyi(n)*(-1)**n
          end do          
       else
          ! n1 ~ -1次は 1 ~ -n1次から計算
          call zbesj(real(z), aimag(z), 1*1.d0, 1, -n1, cyr(n1), cyi(n1), nz, ierr)
          if(ierr /= 0) then
             write(*,*) "Error: ierr=", ierr
             call assert(.false.)
             ! stop
          end if
          ! 順番を反転
          cyr(n1:-1) = cyr(-1:n1:-1)
          cyi(n1:-1) = cyi(-1:n1:-1)
          ! (-1)^nを掛ける
          do n=n1,-1
             cyr(n) = cyr(n)*(-1)**n
             cyi(n) = cyi(n)*(-1)**n
          end do

          ! 0 ~ n2次は普通に計算
          call zbesj(real(z), aimag(z), 0*1.d0, 1, n2+1, cyr(0), cyi(0), nz, ierr)
          if(ierr /= 0) then
             write(*,*) "Error: ierr=", ierr
             call assert(.false.)
             ! stop
          end if          
       end if
    end if

    out(:) = dcmplx(cyr(:),cyi(:))
    
  end subroutine cbesjn

  ! n1~n2次のbessel関数を計算
  subroutine cbesyn(z, n1, n2, out)    
    complex(8),intent(in) :: z
    integer,intent(in) :: n1, n2
    complex(8),intent(out) :: out(n1:n2)

    real(8) :: cyr(n1:n2), cyi(n1:n2)
    integer :: nz, ierr
    integer :: n

    real(8) :: CWRKR(n1:n2), CWRKI(n1:n2) ! work array

    if(n1 >= 0) then  
       call zbesy(real(z), aimag(z), n1*1.d0, 1, n2-n1+1, cyr, cyi, nz, CWRKR, CWRKI, ierr)
       if(ierr /= 0) then
          write(*,*) "Error: ierr=", ierr
          call assert(.false.)
          ! stop
       end if
       
    else
       if(n2 < 0) then          
          ! n1 ~ n2次は -n2 ~ -n1次から計算
          call zbesy(real(z), aimag(z), -n2*1.d0, 1, -n1+n2+1, cyr, cyi, nz, CWRKR, CWRKI, ierr)
          if(ierr /= 0) then
             write(*,*) "Error: ierr=", ierr
             call assert(.false.)
             ! stop
          end if
          ! 順番を反転
          cyr = cyr(n2:n1:-1)
          cyi = cyi(n2:n1:-1)
          ! (-1)^nを掛ける
          do n=n1,n2
             cyr(n) = cyr(n)*(-1)**n
             cyi(n) = cyi(n)*(-1)**n
          end do          
       else
          ! n1 ~ -1次は 1 ~ -n1次から計算
          call zbesy(real(z), aimag(z), 1*1.d0, 1, -n1, cyr(n1), cyi(n1), nz, CWRKR, CWRKI, ierr)
          if(ierr /= 0) then
             write(*,*) "Error: ierr=", ierr
             call assert(.false.)
             ! stop
          end if
          ! 順番を反転
          cyr(n1:-1) = cyr(-1:n1:-1)
          cyi(n1:-1) = cyi(-1:n1:-1)
          ! (-1)^nを掛ける
          do n=n1,-1
             cyr(n) = cyr(n)*(-1)**n
             cyi(n) = cyi(n)*(-1)**n
          end do

          ! 0 ~ n2次は普通に計算
          call zbesy(real(z), aimag(z), 0*1.d0, 1, n2+1, cyr(0), cyi(0), nz, CWRKR, CWRKI, ierr)
          if(ierr /= 0) then
             write(*,*) "Error: ierr=", ierr
             call assert(.false.)
             ! stop
          end if          
       end if
    end if

    out(:) = dcmplx(cyr(:),cyi(:))
    
  end subroutine cbesyn

  ! n1~n2次の第一種hankel関数を計算
  subroutine cbeshn(z, n1, n2, out)    
    complex(8),intent(in) :: z
    integer,intent(in) :: n1, n2
    complex(8),intent(out) :: out(n1:n2)

    complex(8) :: j(n1:n2), y(n1:n2)
    complex(8),parameter :: ione = (0.d0,1.d0)

    call cbesjn(z, n1, n2, j)
    call cbesyn(z, n1, n2, y)

    out = j + ione*y
    
  end subroutine cbeshn

  ! n1~n2次の第二種hankel関数を計算
  subroutine cbeshn_2(z, n1, n2, out)    
    complex(8),intent(in) :: z
    integer,intent(in) :: n1, n2
    complex(8),intent(out) :: out(n1:n2)

    complex(8) :: j(n1:n2), y(n1:n2)
    complex(8),parameter :: ione = (0.d0,1.d0)

    call cbesjn(z, n1, n2, j)
    call cbesyn(z, n1, n2, y)

    out = j - ione*y
    
  end subroutine cbeshn_2

  ! Hankel function of first kind and order 0 and 1
  ! cbesh0 + cbesh1 よりも若干はやい?
  subroutine cbesh0_1(z, out)
    complex(dp),intent(in) :: z
    complex(dp),intent(out) :: out(0:1)
    complex(dp)  :: cbj0, cbj1, cby0, cby1

    call cjy02(z, cbj0, cbj1, cby0, cby1)

    out(0) = cbj0+ci*cby0
    out(1) = cbj1+ci*cby1

    return
    
  end subroutine cbesh0_1

  subroutine dbesh0_1(x, out)
    real(8),intent(in) :: x
    complex(8),intent(out) :: out(0:1)
    
    out(0) = bessel_j0(x) + ci*bessel_y0(x)
    out(1) = bessel_j1(x) + ci*bessel_y1(x)
    
  end subroutine dbesh0_1

  subroutine besk0_1(x, out)
    use iso_c_binding
    real(8),intent(in):: x
    real(8),intent(out),target :: out(0:1)
    real(8) :: bi0
    real(8) :: di0
    real(8) :: bi1
    real(8) :: di1
    real(8) :: bk0
    real(8) :: dk0
    real(8) :: bk1
    real(8) :: dk1
    integer :: err

    ! call ik01a(x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)
    ! out(1) = bk0
    ! out(2) = bk1
    !out(1) = gsl_sf_bessel_K0(x)
    !out(2) = gsl_sf_bessel_K1(x)
    err = gsl_sf_bessel_Kn_array(0,1,x,c_loc(out(0)))
#ifdef DEBUG
    if(err /= 0) then
       write(*,*) "err=",err,"at gsl_sf_bessel_Kn_array"
       ! stop
    end if
#endif
    
    return
  end subroutine besk0_1

  complex(8) function bessel_hn(n,x)
    integer,intent(in) :: n
    real(8),intent(in) :: x

    bessel_hn = dcmplx(bessel_jn(n,x), bessel_yn(n,x))
  end function bessel_hn

  real(8) function besj_real(n, x)
    integer,intent(in) :: n
    real(8),intent(in) :: x

    besj_real = bessel_jn(n, x)
    
  end function besj_real

  complex(8) function besj_complex(n, z)
    integer,intent(in) :: n
    complex(8),intent(in) :: z
    
    complex(8) :: j(1)! , y(1)
    complex(8),parameter :: ione = (0.d0,1.d0)

    call cbesjn(z, n, n, j)
    ! call cbesyn(z, n1, n2, y)

    besj_complex = j(1)
    
  end function besj_complex

  real(8) function besy_real(n, x)
    integer,intent(in) :: n
    real(8),intent(in) :: x

    besy_real = bessel_yn(n, x)
    
  end function besy_real

  complex(8) function besy_complex(n, z)
    integer,intent(in) :: n
    complex(8),intent(in) :: z
    
    complex(8) :: j(1)! , y(1)
    complex(8),parameter :: ione = (0.d0,1.d0)

    call cbesyn(z, n, n, j)
    ! call cbesyn(z, n1, n2, y)

    besy_complex = j(1)
    
  end function besy_complex

  complex(8) function besh_real(n, x)
    integer,intent(in) :: n
    real(8),intent(in) :: x

    besh_real = dcmplx(bessel_jn(n,x), bessel_yn(n,x))
  end function besh_real

  complex(8) function besh_complex(n, z)
    integer,intent(in) :: n
    complex(8),intent(in) :: z
    
    complex(8) :: j(1) , y(1)
    complex(8),parameter :: ione = (0.d0,1.d0)

    call cbesjn(z, n, n, j)
    call cbesyn(z, n, n, y)

    besh_complex = j(1) + ione*y(1)
    
  end function besh_complex

  complex(8) function besh_2_real(n, x)
    integer,intent(in) :: n
    real(8),intent(in) :: x

    besh_2_real = dcmplx(bessel_jn(n,x), -bessel_yn(n,x))
  end function besh_2_real

  complex(8) function besh_2_complex(n, z)
    integer,intent(in) :: n
    complex(8),intent(in) :: z
    
    complex(8) :: j(1) , y(1)
    complex(8),parameter :: ione = (0.d0,1.d0)

    call cbesjn(z, n, n, j)
    call cbesyn(z, n, n, y)

    besh_2_complex = j(1) - ione*y(1)
    
  end function besh_2_complex

  ! Jnの微分
  real(8) function besj_d_real(n, z)
    integer :: n
    real(8) :: z

    besj_d_real = n*besj(n,z)/z - besj(n+1,z)
    
  end function besj_d_real

  ! Ynの微分
  real(8) function besy_d_real(n, z)
    integer :: n
    real(8) :: z

    besy_d_real = n*besy(n,z)/z - besy(n+1,z)
    
  end function besy_d_real

  ! Hnの微分
  complex(8) function besh_d_real(n, z)
    integer :: n
    real(8) :: z

    besh_d_real = n*besh(n,z)/z - besh(n+1,z)
    
  end function besh_d_real

  ! Jnの微分
  complex(8) function besj_d_complex(n, z)
    integer :: n
    complex(8) :: z

    besj_d_complex = n*besj(n,z)/z - besj(n+1,z)
    
  end function besj_d_complex

  ! Ynの微分
  complex(8) function besy_d_complex(n, z)
    integer :: n
    complex(8) :: z

    besy_d_complex = n*besy(n,z)/z - besy(n+1,z)
    
  end function besy_d_complex

  ! Hnの微分
  complex(8) function besh_d_complex(n, z)
    integer :: n
    complex(8) :: z

    besh_d_complex = n*besh(n,z)/z - besh(n+1,z)
    
  end function besh_d_complex

  ! Hnの微分
  complex(8) function besh_2_d_real(n, z)
    integer :: n
    real(8) :: z

    besh_2_d_real = n*besh_2(n,z)/z - besh_2(n+1,z)
    
  end function besh_2_d_real

  complex(8) function besh_2_d_complex(n, z)
    integer :: n
    complex(8) :: z

    besh_2_d_complex = n*besh_2(n,z)/z - besh_2(n+1,z)
    
  end function besh_2_d_complex
end module bessel
