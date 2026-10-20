module vert_diff_exp_mod
  use, intrinsic :: iso_fortran_env
  implicit none


  integer, parameter :: dp = REAL64 ! Precision variable

  real(dp), parameter :: CFL = 0.95_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: kb = 1.380649e-16_dp

  public :: vert_diff_exp
  private :: compute_fluxes

contains

  subroutine vert_diff_exp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, Kzz, nq, q, q0)
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
    real(dp), dimension(nlev) :: alte, pe, Kzze, rhoe
    real(dp), dimension(nlay) :: delz, pl, rho
    real(dp), dimension(nlay-1) :: delz_mid

    real(dp), dimension(nlay,nq) :: k1, k2, k3, k4

    integer :: n_it, n, accept, ierr
    real(dp) :: dt, t_now, dt_max
    real(dp), dimension(nq) :: tol, err
    real(dp), parameter ::  pow = 0.2_dp, safe = 0.9_dp
    real(dp), parameter ::  atol = 1e-99_dp, rtol = 1e-3_dp

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

    !! Find Kzz at levels
    Kzze(1) = Kzz(1)
    do k = 2, nlay
      Kzze(k) = (Kzz(k) + Kzz(k-1))/2.0_dp
    end do
    Kzze(nlev) = Kzz(nlay)

    rho(:) = pl(:) / ((R / mu(:)) * Tl(:))
    !! Find rho at levels
    rhoe(1) = rho(1)
    do k = 2, nlay
      rhoe(k) = (rho(k) + rho(k-1))/2.0_dp
    end do
    rhoe(nlev) = rho(nlay)

    !! Copy tracer array
    do k = 1, nlay
      qc(k,:) = q(k,:)
    end do

    !! Prepare timestepping routine
    !! Find minimum timestep that satifies the CFL condition
    dt_max = t_end
    do k = 1, nlay
      dt_max = min(dt_max,CFL*(delz(k)**2/(2.0_dp*Kzz(k))))
    end do
    dt = dt_max

    !! Begin timestepping routine
    t_now = 0.0_dp
    n_it = 0

    do while ((t_now < t_end) .and. (n_it < 100000))

      !! If next time step overshoots - last time step is equal tend
      if ((t_now + dt) > t_end) then
        dt = t_end - t_now
      end if

      do n = 1, nq
        qc(:,n) = max(qc(:,n),1e-99_dp)
      end do

      call compute_fluxes(nlay, delz_mid(:), delz(:), Kzze(:), rhoe(:), rho(:), nq, q0(:), qc(:,:), k1(:,:))
      q_in(:,:) = qc(:,:) + 0.5_dp * dt * k1(:,:)
      q_in(:,:) = max(q_in(:,:), 1.0e-99_dp)
      call compute_fluxes(nlay, delz_mid(:), delz(:), Kzze(:), rhoe(:), rho(:), nq, q0(:), q_in(:,:), k2(:,:))
      q_in(:,:) = qc(:,:) + 0.75_dp * dt * k2(:,:)
      q_in(:,:) = max(q_in(:,:), 1.0e-99_dp)
      call compute_fluxes(nlay, delz_mid(:), delz(:), Kzze(:), rhoe(:), rho(:), nq, q0(:), q_in(:,:), k3(:,:))

      !! Update q_new 
      q_new(:,:) = qc(:,:) + dt * (2.0_dp / 9.0_dp * k1(:,:) + 3.0_dp / 9.0_dp * k2(:,:) + 4.0_dp / 9.0_dp * k3(:,:))

      q_new(:,:) = max(q_new(:,:), 1.0e-99_dp)

      !! Calculated q_embedded
      call compute_fluxes(nlay, delz_mid(:), delz(:), Kzze(:), rhoe(:), rho(:), nq, q0(:), q_new(:,:), k4(:,:))
      q_em(:,:) = qc(:,:) + dt * (7.0_dp / 24.0_dp * k1(:,:) + 1.0_dp / 4.0_dp * k2(:,:) + 1.0_dp / 3.0_dp * k3(:,:) &
        & + 1.0_dp/8.0_dp * k4(:,:))
      q_em(:,:) = max(q_em(:,:), 1.0e-99_dp)

      !! Compute error
      accept = 0
      ierr = 0
      do n = 1, nq

        err(n) = maxval(abs(q_new(:,n) - q_em(:,n)))
        tol(n) = atol + rtol * maxval(abs(q_new(:,n)))

        if (n == 1 .or. err(n)/max(tol(n), tiny(1.0_dp)) > err(ierr)/max(tol(ierr), tiny(1.0_dp))) then
          ierr = n
        end if

        if (err(n) > tol(n)) then
          accept = -1
        end if

      end do

      if (accept == 0) then
        ! Accept the step
        qc(:,:) = q_new(:,:)
        t_now = t_now + dt
        n_it = n_it + 1

        ! Adjust timestep with safety factor
        dt = safe * dt * (tol(ierr) / max(err(ierr), tiny(1.0_dp)))**pow
        dt = min(dt, dt_max)
      else
        ! Reject the step and reduce timestep
        dt = safe * dt * (tol(ierr) / max(err(ierr), tiny(1.0_dp)))**pow
      end if

      ! Escape condition for small timestep
      if (dt < 1.0e-10_dp) then
        print *, "Exiting: timestep too small: ", dt
        exit
      end if


    end do

    do k = 1, nlay
      q(k,:) = qc(k,:)
    end do

    do n = 1, nq
      q(:,n) = max(q(:,n),1e-99_dp)
    end do

  end subroutine vert_diff_exp

  subroutine compute_fluxes(nlay, delz_mid, delz, Kzze, rhoe, rho, nq, q0, q_in, flux)
    implicit none

    integer, intent(in) :: nlay, nq
    real(dp), dimension(nlay-1), intent(in) :: delz_mid
    real(dp), dimension(nlay), intent(in) :: delz, rho
    real(dp), dimension(nlay+1), intent(in) :: Kzze, rhoe
    real(dp), dimension(nq), intent(in) :: q0
    real(dp), dimension(nlay, nq), intent(in) :: q_in

    real(dp), dimension(nlay, nq), intent(out) :: flux

    integer :: k, n
    real(dp), dimension(nlay+1, nq) :: face_flux

    !! Form diffusive face fluxes.  The top face is zero-flux; q0 is a
    !! Dirichlet value located at the bottom face, half a cell from q(nlay).
    do n = 1, nq
      face_flux(1,n) = 0.0_dp
      do k = 2, nlay
        face_flux(k,n) = rhoe(k)*Kzze(k)*(q_in(k,n) - q_in(k-1,n))/delz_mid(k-1)
      end do
      face_flux(nlay+1,n) = rhoe(nlay+1)*Kzze(nlay+1)* &
                            (q0(n) - q_in(nlay,n))/(0.5_dp*delz(nlay))

      do k = 1, nlay
        flux(k,n) = (face_flux(k+1,n) - face_flux(k,n))/(rho(k)*delz(k))
      end do
    end do

  end subroutine compute_fluxes

end module vert_diff_exp_mod
