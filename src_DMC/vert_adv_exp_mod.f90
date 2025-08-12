module vert_adv_exp_mod
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: dp = REAL64
  private :: koren_phi, adv_rhs_muscl_down_nugrid
  public  :: vert_adv_exp

contains

  subroutine vert_adv_exp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q)
    ! TVD MUSCL with RK2 (SSPRK(2,2)), non-uniform grid consistent.
    implicit none
    integer, intent(in)            :: nlay, nlev, nq
    real(dp), intent(in)           :: t_end, grav_in
    real(dp), intent(in)           :: Tl(nlay), pl_in(nlay), mu(nlay), pe_in(nlev)
    real(dp), intent(in)           :: vf(nlay,nq)
    real(dp), intent(inout)        :: q(nlay,nq)

    integer  :: i, k, n
    real(dp) :: R_gas, grav, eps, CFL
    real(dp) :: t, dt, dt_max, num, den
    real(dp), dimension(nq) :: q_top
    real(dp), allocatable :: alte(:), altm(:), dz(:), dzm(:)
    real(dp), allocatable :: pl(:), pe(:), rho(:), rho_e(:), v_e(:), A(:)
    real(dp), allocatable :: rhs(:), q_stage(:)

    CFL   = 0.95_dp
    q_top(:) = 0.0_dp
    R_gas = 8.31446261815324e7_dp
    grav  = 100.0_dp * grav_in
    eps   = 1.0e-30_dp

    allocate(alte(nlev), altm(nlay), dz(nlay), dzm(nlay-1))
    allocate(pl(nlay), pe(nlev), rho(nlay), rho_e(nlev))
    allocate(v_e(nlev), A(nlev))
    allocate(rhs(nlay), q_stage(nlay))

    ! geometry (hydrostatic heights)
    pl = 10.0_dp*pl_in; pe = 10.0_dp*pe_in
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R_gas*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do
    do i = 1, nlay
      dz(i)   = alte(i) - alte(i+1)
      altm(i) = 0.5_dp*(alte(i) + alte(i+1))
    end do
    do i = 1, nlay-1
      dzm(i) = altm(i) - altm(i+1)    ! center-to-center spacing
    end do

    ! densities at centers/faces
    rho(:)    = pl(:) / ((R_gas / mu(:)) * Tl(:))
    rho_e(1)  = rho(1)
    do i = 1, nlay-1
      rho_e(i+1) = 0.5_dp*(rho(i)+rho(i+1))
    end do
    rho_e(nlev) = rho(nlay)

    do n = 1, nq
      ! face velocities (downward/positive)
      v_e(1) = vf(1,n)
      do i = 1, nlay-1
        v_e(i+1) = 0.5_dp*(vf(i,n) + vf(i+1,n))
      end do
      v_e(nlev) = vf(nlay,n)
      do k = 1, nlev
        A(k) = rho_e(k) * v_e(k)
      end do

      ! dt from upwind Courant (face k uses upwind cell i=k-1)
      dt_max = t_end
      do k = 2, nlay
        num    = rho(k-1) * dz(k-1)
        den    = A(k) + eps
        dt_max = min(dt_max, CFL * num / den)
      end do

      t = 0.0_dp
      do while (t < t_end)
        dt = merge(dt_max, t_end - t, t + dt_max <= t_end)

        ! Stage 1
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q(:,n), q_top(n), eps, rhs)
        do i = 1, nlay
          q_stage(i) = q(i,n) + dt * rhs(i)
          if (q_stage(i) < 1.0e-99_dp) q_stage(i) = 1.0e-99_dp
        end do

        ! Stage 2
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q_stage, q_top(n), eps, rhs)
        do i = 1, nlay
          q(i,n) = 0.5_dp*q(i,n) + 0.5_dp*( q_stage(i) + dt*rhs(i) )
          if (q(i,n) < 1.0e-99_dp) q(i,n) = 1.0e-99_dp
        end do

        t = t + dt
      end do
    end do

    deallocate(alte, altm, dz, dzm, pl, pe, rho, rho_e, v_e, A, rhs, q_stage)
  end subroutine vert_adv_exp

  ! ===== Non-uniform grid MUSCL RHS (downward flow A>0) =====
  subroutine adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q, q_top, eps, rhs)
    implicit none
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: dz(nlay), dzm(nlay-1), rho(nlay), A(nlev)
    real(dp), intent(in) :: q(nlay), q_top, eps
    real(dp), intent(out):: rhs(nlay)

    integer :: i, k
    real(dp) :: gL, gR, r, phi
    real(dp), allocatable :: sigma(:), qL_face(:), qR_face(:), F(:)


    ! --- near top boundary (A>0 inflow) ---
    real(dp) :: dzhalf_top, gL_top, gR1, r1, phi1, qR1

    allocate(sigma(nlay), qL_face(nlay), qR_face(nlay), F(nlev))

    ! Limited slope sigma_i ~ dq/dz (use forward gradient for positive flow)
    sigma(1)     = 0.0_dp
    sigma(nlay)  = 0.0_dp
    do i = 2, nlay-1
      gL = (q(i)   - q(i-1)) / (dzm(i-1) + eps)   ! left one-sided gradient
      gR = (q(i+1) - q(i)  ) / (dzm(i)   + eps)   ! right one-sided gradient
      r  = gL / (gR + eps)
      phi = koren_phi(r)
      sigma(i) = phi * gR            ! pick the downwind gradient for positive flow
    end do

    ! Cell-face reconstructions using true half distances
    ! Left face value of cell i:  q_i^L = q_i - 0.5*dzm(i-1)*sigma(i)
    ! Right face value of cell i: q_i^R = q_i + 0.5*dzm(i)  *sigma(i)
    qL_face(1) = q(1)                 ! unused (no left interior face)
    do i = 2, nlay-1
      qL_face(i) = q(i) - 0.5_dp*dzm(i-1)*sigma(i)
      qR_face(i) = q(i) + 0.5_dp*dzm(i)  *sigma(i)
    end do
    ! end cells: one-sided (no reconstruction beyond domain)
    qR_face(1)    = q(1) + 0.5_dp*dzm(1)   *sigma(1)
    qL_face(nlay) = q(nlay) - 0.5_dp*dzm(nlay-1)*sigma(nlay)
    qR_face(nlay) = q(nlay)                ! unused

    dzhalf_top = 0.5_dp*dz(1)                          ! distance center(1) -> top face
    gL_top     = (q(1) - q_top) / (dzhalf_top + eps)   ! upwind-side gradient using ghost=inflow
    gR1        = (q(2) - q(1)) / (dzm(1) + eps)
    r1         = gL_top / (gR1 + eps)
    phi1       = koren_phi(r1)

    ! limited slope in cell 1 (use downwind grad scaled by limiter)
    sigma(1)   = phi1 * gR1

    ! reconstruct first cell right-face value
    qR_face(1) = q(1) + 0.5_dp*dzm(1)*sigma(1)

    ! enforce monotonicity between inflow and neighbour (Barth–Jespersen style clamp)
    qR_face(1) = max( min(qR_face(1), max(q_top, q(1), q(2))), &
                      min(q_top, q(1), q(2)) )

    ! Upwind fluxes (A>0 ⇒ use upwind cell's right face value)
    F(1) = A(1) * q_top            ! top inflow
    do k = 2, nlay
      ! face between cells (k-1) and k
      F(k) = A(k) * qR_face(k-1)
    end do
    F(nlev) = 0.0_dp               ! bottom outer face not used in divergence

    ! Conservative update RHS
    do i = 1, nlay
      rhs(i) = - (F(i+1) - F(i)) / (rho(i)*dz(i) + eps)
    end do

    deallocate(sigma, qL_face, qR_face, F)
  end subroutine adv_rhs_muscl_down_nugrid

  pure real(dp) function koren_phi(r) result(phi)
    real(dp), intent(in) :: r
    phi = max(0.0_dp, min(2.0_dp*r, (1.0_dp + 2.0_dp*r)/3.0_dp, 2.0_dp))
  end function koren_phi

end module vert_adv_exp_mod
