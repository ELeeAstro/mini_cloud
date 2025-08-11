module vert_adv_exp_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  private :: koren_phi
  public :: vert_adv_exp

contains

  pure real(dp) function koren_phi(r) result(phi)
    real(dp), intent(in) :: r
    phi = max(0.0_dp, min(2.0_dp*r, (1.0_dp + 2.0_dp*r)/3.0_dp, 2.0_dp))
  end function koren_phi

  subroutine vert_adv_exp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, &
                             vf, nq, q)
    ! TVD flux-limited Lax–Wendroff advection, downward-only (vf > 0).
    ! Inputs:
    !   vf(nlay,nq)  : settling speed per tracer (cm/s), positive downward
    !   q(nlay,nq)   : in/out, cell-avg tracer (mixing ratio)
    !   q_top(nq)    : Dirichlet inflow value at the top face
    ! Notes: top inflow via F(1)=rho_e(1)*v_e(1)*q_top; bottom is outflow.
    implicit none
    integer, intent(in)            :: nlay, nlev, nq
    real(dp), intent(in)           :: t_end, grav_in
    real(dp), intent(in)           :: Tl(nlay), pl_in(nlay), mu(nlay), pe_in(nlev)
    real(dp), intent(in)           :: vf(nlay,nq)
    real(dp), intent(inout)        :: q(nlay,nq)

    integer  :: i, k, n
    real(dp) :: R_gas, grav, eps, CFL
    real(dp) :: t, dt, dt_max, vmax, dzmin, num, den, Cup
    real(dp), dimension(nq) :: q_top
    real(dp), allocatable :: alte(:), altm(:), dz(:), dzm(:)
    real(dp), allocatable :: pl(:), pe(:), rho(:), rho_e(:), v_e(:), A(:), F(:)
    real(dp) :: q_up, q_dn, dq_up, dq_dn, r, phi

    CFL = 0.95_dp
    q_top(:) = 0.0_dp

    R_gas    = 8.31446261815324e7_dp
    grav = 100.0_dp * grav_in
    eps  = 1.0e-30_dp

    allocate(alte(nlev), altm(nlay), dz(nlay), dzm(nlay-1))
    allocate(pl(nlay), pe(nlev), rho(nlay), rho_e(nlev), v_e(nlev))
    allocate(A(nlev), F(nlev))

    ! Hydrostatic heights (cm) and spacings (your existing pattern)
    pl = 10.0_dp * pl_in
    pe = 10.0_dp * pe_in
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R_gas*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do
    do i = 1, nlay
      dz(i) = alte(i) - alte(i+1)
      altm(i) = 0.5_dp*(alte(i) + alte(i+1))
    end do
    do i = 1, nlay-1
      dzm(i) = altm(i) - altm(i+1)
    end do

    ! Density at centers and faces
    rho(:)    = pl(:) / ((R_gas / mu(:)) * Tl(:))
    rho_e(1)  = rho(1)
    do i = 1, nlay-1
      rho_e(i+1) = 0.5_dp*(rho(i)+rho(i+1))
    end do
    rho_e(nlev) = rho(nlay)

    ! Time step from upwind Courant on faces (downward-only => upwind cell is i-1 at face f=i)
    
    do n = 1, nq
      ! face-avg velocity for this tracer
      v_e(1) = vf(1,n)
      do i = 1, nlay-1
        v_e(i+1) = 0.5_dp*(vf(i,n) + vf(i+1,n))
      end do
      v_e(nlev) = vf(nlay,n)

      dt_max = t_end
      do k = 2, nlay
        num = rho(k-1) * dz(k-1)
        den = rho_e(k) * v_e(k) + eps
        dt_max = min(dt_max, CFL * num / den)
      end do

      t = 0.0_dp
      do while (t < t_end)
        if (t + dt_max <= t_end) then
          dt = dt_max
        else
          dt = t_end - t
        end if

        ! Build mass flux A = rho_e * v_e for this tracer
        do k = 1, nlev
          A(k) = rho_e(k) * v_e(k)
        end do

        ! Top face: inflow flux
        F(1) = A(1) * q_top(n)

        ! Interior faces f=2..nlay (between cells i-1 and i)
        do k = 2, nlay
          i = k              ! face index equals right cell index
          ! Upwind/downwind states (downward ⇒ upwind is i-1)
          q_up = q(i-1,n)
          q_dn = q(i  ,n)

          ! Upwind-consistent slope ratio r_f on nonuniform grid
          if (i-2 >= 1) then
            dq_up = (q(i-1,n) - q(i-2,n)) / (dzm(i-2) + eps)
          else
            dq_up = 0.0_dp
          end if
          dq_dn = (q(i,n)   - q(i-1,n)) / (dzm(i-1) + eps)
          r     = dq_up / (dq_dn + eps)
          phi   = koren_phi(r)   ! or provide an 'mc_phi' if you prefer

          ! Upwind Courant based on upwind cell (i-1)
          Cup = (A(k) * dt) / (rho(i-1) * dz(i-1) + eps)
          if (Cup > 1.0_dp) Cup = 1.0_dp
          if (Cup < 0.0_dp) Cup = 0.0_dp

          ! High-resolution TVD flux at face f
          F(k) = A(k)*q_up + 0.5_dp*A(k)*(1.0_dp - Cup) * phi * (q_dn - q_up)
        end do

        ! Bottom outer face not used in divergence
        F(nlev) = 0.0_dp

        ! Conservative update of q(:,n)
        do i = 1, nlay
          q(i,n) = q(i,n) - dt * (F(i+1) - F(i)) / (rho(i)*dz(i) + eps)
          if (q(i,n) < 1.0e-99_dp) q(i,n) = 1.0e-99_dp
        end do

        ! Optional: if you insist on fixed bottom value, enforce it here (not recommended for outflow)
        ! q(nlay,n) = q_bot(n)

        t = t + dt
      end do
    end do

    deallocate(alte, altm, dz, dzm, pl, pe, rho, rho_e, v_e, A, F)

  end subroutine vert_adv_exp

end module vert_adv_exp_mod
