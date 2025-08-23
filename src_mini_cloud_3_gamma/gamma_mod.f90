module gamma_func_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: TOL  = 1.0e-9_dp
  real(dp), parameter :: TINY = 1.0e-300_dp  ! guard for divisions
  integer,  parameter :: ITMAX = 100000

  public  :: up_inc_gam           ! unregularized upper Γ(a,x)
  private ::  gamma_p, gamma_q, series_p, cf_q

contains

  function gamma_p(a, x) result(P)
    !! Regularized lower incomplete gamma P(a,x)
    real(dp), intent(in) :: a, x
    real(dp) :: P, lgam_a, ax

    if (a <= 0.0_dp .or. x < 0.0_dp) then
     stop
    end if
    if (x == 0.0_dp) then
      P = 0.0_dp; return
    end if

    lgam_a = log_gamma(a)
    ax     = a*log(x) - x - lgam_a

    if (x < a + 1.0_dp) then
      P = series_p(a, x, ax)
    else
      P = 1.0_dp - cf_q(a, x, ax)
    end if

    ! Clamp to [0,1] in case of round-off
    if (P < 0.0_dp) P = 0.0_dp
    if (P > 1.0_dp) P = 1.0_dp
  end function gamma_p

  function gamma_q(a, x) result(Q)
    !! Regularized upper incomplete gamma Q(a,x)
    real(dp), intent(in) :: a, x
    real(dp) :: Q, lgam_a, ax

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      stop
    end if
    if (x == 0.0_dp) then
      Q = 1.0_dp; return
    end if

    lgam_a = log_gamma(a)
    ax     = a*log(x) - x - lgam_a

    if (x < a + 1.0_dp) then
      Q = 1.0_dp - series_p(a, x, ax)
    else
      Q = cf_q(a, x, ax)
    end if

    if (Q < 0.0_dp) Q = 0.0_dp
    if (Q > 1.0_dp) Q = 1.0_dp

  end function gamma_q

  function up_inc_gam(a, x) result(Gu)
    !! Unregularized upper incomplete gamma Γ(a,x) = Γ(a) * Q(a,x)
    real(dp), intent(in) :: a, x
    real(dp) :: Gu
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      stop
    end if
    Gu = exp(log_gamma(a)) * gamma_q(a, x)
  end function up_inc_gam

  !---------------- internal kernels ----------------

  pure function series_p(a, x, ax) result(P)
    !! Power series for P(a,x), using ax = a*log(x) - x - log_gamma(a)
    real(dp), intent(in) :: a, x, ax
    real(dp) :: P
    real(dp) :: term, sum, an
    integer  :: n

    term = 1.0_dp / a
    sum  = term
    n = 1
    do
      an   = a + real(n, dp)
      term = term * (x / an)
      sum  = sum + term
      if (abs(term) <= TOL * abs(sum)) exit
      n = n + 1
      if (n > ITMAX) exit
    end do
    P = exp(ax) * sum
  end function series_p

  pure function cf_q(a, x, ax) result(Q)
    !! Modified Lentz continued fraction for Q(a,x), using ax as above
    real(dp), intent(in) :: a, x, ax
    real(dp) :: Q
    real(dp) :: b, c, d, f, delta, an
    integer  :: n

    ! Start values for the standard CF of Q(a,x)
    b = x + 1.0_dp - a
    d = 1.0_dp / max(b, TINY)
    c = 1.0_dp / TINY
    f = d

    n = 1
    do
      an = real(n, dp) * (a - real(n, dp))
      b  = b + 2.0_dp
      d  = 1.0_dp / max(b + an*d, TINY)
      c  = max(b + an/c, TINY)
      delta = c * d
      f = f * delta
      if (abs(delta - 1.0_dp) <= TOL) exit
      n = n + 1
      if (n > ITMAX) exit
    end do

    Q = exp(ax) * f
  end function cf_q

end module gamma_func_mod
