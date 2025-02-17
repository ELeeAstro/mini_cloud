module vert_diff_imp_mod
  use, intrinsic :: iso_fortran_env
  implicit none


  integer, parameter :: dp = REAL64 ! Precision variable

  real(dp), parameter :: CFL = 0.90_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp

  public :: vert_diff_imp

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

    integer :: k
    real(dp) :: grav
    real(dp), dimension(nlev) :: alte, pe, Kzze
    real(dp), dimension(nlay) :: delz, delz_mid, pl

    real(dp), dimension(nlay-1) :: ld, ud, c_p
    real(dp), dimension(nlay) :: md, rhs, d_p

    integer :: n_it, n, ierr, i
    real(dp) :: dt, t_now, alpha, beta, denom, dx_f, dx_b

    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2

    grav = grav_in * 100.0_dp


    !! First calculate the vertical height (cm) assuming hydrostatic equilibrium and differences
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
      delz(k) = alte(k) - alte(k+1)
    end do

    !! Find differences between layers directly
    do k = 1, nlay-1
      delz_mid(k) = delz(k)/2.0_dp + delz(k+1)/2.0_dp
    end do
    delz_mid(nlay) = delz(nlay)

    !! Find Kzz at levels
    Kzze(1) = Kzz(1)
    do k = 2, nlay
      Kzze(k) = (Kzz(k) + Kzz(k-1))/2.0_dp
    end do
    Kzze(nlev) = Kzz(nlay)


    qc(:,:) = q(:,:)

    dt = t_end


    !! Perform inversion for each tracer
    do n = 1, nq


    end do

    do n = 1, nq
      q(:,n) = max(q(:,n),1e-30_dp)
    end do

  end subroutine vert_diff_imp

end module vert_diff_imp_mod
