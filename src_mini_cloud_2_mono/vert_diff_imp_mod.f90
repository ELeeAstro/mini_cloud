module vert_diff_imp_mod
  use, intrinsic :: iso_fortran_env
  implicit none


  integer, parameter :: dp = REAL64 ! Precision variable

  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: kb = 1.380649e-16_dp

  public :: vert_diff_imp
  private :: thomas

contains

  subroutine vert_diff_imp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, Kzz, nq, q, q0)
    implicit none

    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), dimension(nlay), intent(in) :: Tl, pl_in, Kzz, mu
    real(dp), dimension(nlev), intent(in) :: pe_in
    real(dp), dimension(nq), intent(in) :: q0

    real(dp), dimension(nlay,nq), intent(inout) :: q
    
    real(dp), dimension(nlay,nq) :: qc, q_new, q_em, q_in

    integer :: k, n
    real(dp) :: grav, q_min, inv_dt, scale, theta
    real(dp), dimension(nlev) :: alte, pe, K_e, rho_e, D, J_old
    real(dp), dimension(nlay) :: dz, pl, rho, altm, a, b, c, rhs, scales
    real(dp), dimension(nlay-1) :: dz_m

    theta = 0.6_dp

    q_min = 1e-99_dp

    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2

    grav = grav_in * 100.0_dp

    rho(:) = pl(:) / ((R / mu(:)) * Tl(:))

    K_e(1)   = Kzz(1)
    rho_e(1) = rho(1)
    do k = 1, nlay-1
      K_e(k+1)   = 0.5_dp * (Kzz(k) + Kzz(k+1))
      rho_e(k+1) = 0.5_dp * (rho(k) + rho(k+1))
    end do
    K_e(nlev)  = Kzz(nlay)
    rho_e(nlev) = rho(nlay)

    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do
    do k = 1, nlay
      altm(k) = 0.5_dp * (alte(k) + alte(k+1))
    end do
    do k = 1, nlay
      dz(k) = alte(k) - alte(k+1)
    end do
    do k = 1, nlay-1
      dz_m(k) = altm(k) - altm(k+1)
    end do

    D(:) = 0.0_dp
    do k = 2, nlay
      D(k) = rho_e(k) * K_e(k) / (dz_m(k-1) + 1e-300_dp)
    end do

    scales(:) = 1.0_dp / (rho(:) * dz(:))

    inv_dt = 1.0_dp / t_end

    do n = 1, nq

      J_old(1) = 0.0_dp
      do k = 2, nlev-1
        J_old(k) = -D(k) * (q(k,n) - q(k-1,n))
      end do
      J_old(nlev) = - D(nlay) * ( q0(n) - q(nlay-1,n) )   ! bottom face uses Dirichlet q0

      do k = 1, nlay-1

        a(k) = -theta * scales(k) * D(k)
        c(k) = -theta * scales(k) * D(k+1)
        b(k) = inv_dt + theta * scales(k) * (D(k) + D(k+1))

        rhs(k) = inv_dt * q(k,n) - (1.0_dp - theta) * scales(k) * (J_old(k+1) - J_old(k))
      end do

      a(nlay) = 0.0_dp
      b(nlay) = 1.0_dp
      c(nlay) = 0.0_dp
      rhs(nlay) = q0(n)

      call thomas(a, b, c, rhs, q(:,n))
      q(:,n) = max(q(:,n),q_min)

    end do

  end subroutine vert_diff_imp

  subroutine thomas(a, b, c, d, x)
    implicit none
    real(dp), intent(in)  :: a(:), b(:), c(:), d(:)
    real(dp), intent(out) :: x(:)
    integer               :: n, i
    real(dp), allocatable :: cp(:), dpv(:)
    real(dp)              :: denom, eps

    n = size(b)
    if (size(a) /= n .or. size(c) /= n .or. size(d) /= n .or. size(x) /= n) then
      error stop "thomas_solve: array size mismatch"
    end if

    allocate(cp(n), dpv(n))
    eps = 1.0e-300_dp  ! tiny to avoid division by zero

    ! Forward sweep
    if (abs(b(1)) <= eps) error stop "thomas_solve: zero pivot at i=1"
    cp(1)  = c(1) / b(1)
    dpv(1) = d(1) / b(1)

    do i = 2, n
      denom = b(i) - a(i) * cp(i-1)
      if (abs(denom) <= eps) error stop "thomas_solve: zero pivot in forward sweep"
      if (i < n) then
        cp(i) = c(i) / denom
      else
        cp(i) = 0.0_dp      ! c(n) is not used
      end if
      dpv(i) = (d(i) - a(i) * dpv(i-1)) / denom
    end do

    ! Back substitution
    x(n) = dpv(n)
    do i = n-1, 1, -1
      x(i) = dpv(i) - cp(i) * x(i+1)
    end do

    deallocate(cp, dpv)
  end subroutine thomas

end module vert_diff_imp_mod

