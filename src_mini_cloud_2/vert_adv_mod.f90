module vert_adv_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: CFL = 0.95_dp
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp

  public :: vert_adv
  private :: minmod

  contains
 
  subroutine vert_adv(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q, q0)
    implicit none


    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), dimension(nlay), intent(in) :: Tl, pl_in, vf, mu
    real(dp), dimension(nlev), intent(in) :: pe_in
    real(dp), dimension(nq), intent(in) :: q0

    real(dp), dimension(nlay,nq), intent(inout) :: q

    integer :: k
    real(dp) :: h1, h2, grav
    real(dp), dimension(nlev) :: alte, lpe, vf_e, Te, nde, pe
    real(dp), dimension(nlay) :: delz, delz_mid, pl

    real(dp), dimension(nlay,nq) :: qc
    real(dp), dimension(nlay) :: sig, c

    integer :: n_it, n
    real(dp) :: dt, t_now


    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2

    grav = grav_in * 100.0_dp


    ! Find velocity at levels
    vf_e(1) = vf(1)
    do k = 2, nlay
      vf_e(k) = (vf(k-1) + vf(k))/2.0_dp
    end do
    vf_e(nlev) = vf(nlay)

    !! First calculate the vertical height (cm) assuming hydrostatic equilibrium and differences
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (kb*Tl(k))/(mu(k) * amu * grav) * log(pe(k+1)/pe(k))
      delz(k) = alte(k) - alte(k+1)
    end do

    !! Find differences between layers directly
    do k = 1, nlay-1
      delz_mid(k) = (alte(k) + alte(k+1))/2.0_dp - (alte(k+1) + alte(k+2))/2.0_dp
    end do
    delz_mid(nlay) = delz(nlay)

    !! Find minimum timestep that allows the CFL condition
    dt = t_end
    do k = 1, nlay
      dt = min(dt,CFL*(delz_mid(k))/abs(vf_e(k+1)))
    enddo

    t_now = 0.0_dp
    n_it = 0

    do while ((t_now < t_end) .and. (n_it < 1000000))

      ! If next time step overshoots - last time step is equal tend
      if (t_now + dt > t_end) then
        dt = t_end - t_now
      end if

      !! Find the courant number
      c(:) = abs(vf_e(2:nlev)) * dt / delz_mid(:)

      do n = 1, nq

        !! Find the minmod limiter
        call minmod(nlay,q(:,n),delz_mid,sig)

        !! Perform McCormack step
        qc(:,n) = q(:,n)
        qc(:nlay-1,n) = q(:nlay-1,n) - sig(:nlay-1)*c(:nlay-1)*(q(2:nlay,n) - q(:nlay-1,n))
        q(2:nlay,n) = 0.5_dp * (q(2:nlay,n) + qc(2:nlay,n) - c(2:nlay)*(qc(2:nlay,n) - qc(:nlay-1,n)))

        q(:,n) = max(q(:,n),1e-30_dp)

      end do

      t_now = t_now + dt
      n_it = n_it + 1

      ! Apply boundary conditions
      q(nlay,:) = q0(:)

    end do
    
  end subroutine vert_adv

  subroutine minmod(nlay,q,dz,sig)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q, dz

    real(dp), dimension(nlay), intent(out) :: sig

    integer :: i
    real(dp) :: de_minus, de_plus

    sig(1) = 0.0_dp
    do i = 2, nlay-1
      de_minus = (q(i) - q(i-1)) / dz(i)
      de_plus = (q(i+1) - q(i)) / dz(i)
      if ((de_minus > 0.0_dp) .and. (de_plus > 0.0_dp)) then
        sig(i) = min(de_minus, de_plus)
      else if ((de_minus < 0.0_dp) .and. (de_plus < 0.0_dp)) then
        sig(i) = max(de_minus, de_plus)
      else
        sig(i) = 0.0_dp
      end if
    end do
    sig(nlay) = 0.0_dp

  end subroutine minmod

end module vert_adv_mod
