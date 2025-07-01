module gamma_func_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  public :: up_inc_gam
  private :: gser, gcf

contains

  function up_inc_gam(a,x) result(inc_gam)
    ! returns the upper incomplete gamma function (a,x)
    implicit none

    real(dp), intent(in) :: a, x
    real(dp) :: inc_gam, gln, gamser, gamcf

    if ((x < 0.0_dp) .or. (a <= 0.0_dp)) then
      print*,'Invalid arguments to gammq: ', x, a
      stop
    end if

    if (x < a + 1.0_dp) then
      ! use series representation
      call gser(gamser,a,x,gln)
      inc_gam = (1.0_dp - gamser) * exp(gln) ! and take its complement
    else
      ! use continued fraction representation
      call gcf(gamcf,a,x,gln)
      inc_gam = gamcf * exp(gln)
    end if

  end function up_inc_gam

  subroutine gser(gamser, a, x, gln)
    implicit none
    ! Returns the normalised lower incomplete gamma function and log gamma

    real(dp), intent(in) :: a, x
    real(dp), intent(out) :: gamser, gln

    integer, parameter :: itmax = 100
    real(dp), parameter :: eps = 3.0e-7_dp

    integer :: n, ierr
    real(dp) :: ap, sum, del

    gln = log_gamma(a)

    if (x <= 0.0_dp) then
      if (x < 0.0_dp) then
        print*, 'Invalid value of x in gser: ', x
        stop 
      end if
      gamser = 0.0_dp
      return
    end if

    ap = a
    sum = 1.0_dp/a
    del = sum
    ierr = 1
    do n = 1, itmax
      ap = ap + 1.0_dp
      del = del * x/ap
      sum = sum + del
      if (abs(del) < abs(sum)*eps) then
        ierr = 0
        exit
      end if
    end do

    if (ierr == 1) then
      print*, 'not enough iterations in gser for large a: ', a, itmax
      stop
    end if

    gamser = sum * exp(-x+a*log(x) - gln)

  end subroutine gser

  subroutine gcf(gamcf, a, x, gln)
    implicit none
    ! Returns the normalised upper incomplete gamma function and log gamma

    real(dp), intent(in) :: a, x
    real(dp), intent(out) :: gamcf, gln

    integer, parameter :: itmax = 100
    real(dp), parameter :: eps = 3.0e-7_dp

    integer :: n, ierr
    real(dp) :: gold, a0, a1, b0, b1, fac, an, anf, g, ana

    gln = log_gamma(a)

    gold = 0.0_dp
    a0 = 1.0_dp
    a1 = x
    b0 = 0.0_dp
    b1 = 1.0_dp
    fac = 1.0_dp

    ierr = 1
    do n = 1, itmax
      an = real(n,dp)
      ana = an - a
      a0 = (a1 + a0*ana)*fac
      b0 = (b1 + b0*ana)*fac
      anf = an*fac
      a1 = x*a0 + anf*a1
      b1 = x*b0 + anf*b1
      if (a1 /= 0.0_dp) then
        fac = 1.0_dp/a1
        g = b1*fac
        if (abs((g-gold)/g) < eps) then
          ierr = 0
          exit
        else
          gold = g
        end if
      end if
    end do

    if (ierr == 1) then
      print*, 'not enough iterations in gcf for large a: ', a, itmax
      stop
    end if

    gamcf = exp(-x+a*log(x) - gln)*g

  end subroutine gcf

end module gamma_func_mod