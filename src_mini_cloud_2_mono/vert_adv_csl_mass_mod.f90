module vert_adv_csl_mass_mod
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: R_GAS = 8.31446261815324e7_dp
  real(dp), parameter :: TINY_W = tiny(1.0_dp)

  integer, parameter, public :: CSLM_OK = 0
  integer, parameter, public :: CSLM_WARN_Q_CLIPPED = 1
  integer, parameter, public :: CSLM_WARN_V_CLIPPED = 2
  integer, parameter, public :: CSLM_ERR_DIMENSIONS = -1
  integer, parameter, public :: CSLM_ERR_TIMESTEP = -2
  integer, parameter, public :: CSLM_ERR_GRID = -3
  integer, parameter, public :: CSLM_ERR_TRAJECTORY = -4
  integer, parameter, public :: CSLM_ERR_REMAP = -5
  integer, parameter, public :: CSLM_ERR_CONSERVATION = -6
  integer, parameter, public :: CSLM_ERR_NONFINITE = -7

  private
  public :: vert_adv_csl_mass

contains

  subroutine vert_adv_csl_mass(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q, &
                          mass_out, mass_clipped, courant_max, ierr)
    !! Conservative semi-Lagrangian transport for downward settling in the
    !! hydrostatic column-mass coordinate m=(p-p_top)/g.
    !!
    !! The cell-integrated conserved quantity is dm*q, where hydrostatic layer
    !! mass dm=(p_bottom-p_top)/g is calculated directly from pressure edges.
    !! Arrival-cell faces are traced backwards through a frozen,
    !! piecewise-linear mass-flux field W=rho*v.  The old tracer mass between
    !! adjacent departure faces is remapped using a monotonically limited
    !! linear reconstruction of q in m.  Interpolating W directly preserves
    !! a constant Epstein mass flux exactly at layer faces.
    !!
    !! The top boundary has zero tracer inflow.  Tracer crossing the bottom
    !! boundary is removed and, when requested, returned in mass_out.
    !! mass_out and mass_clipped have column-mass units (dm*q).  vf is in
    !! cm s-1 and positive velocity points downward.  ierr is zero on success,
    !! a bit mask of CSLM_WARN_* values for recoverable clipping, or a negative
    !! CSLM_ERR_* value for a fatal input/remapping error.
    implicit none

    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), intent(in) :: mu(nlay), Tl(nlay), pl_in(nlay), pe_in(nlev)
    real(dp), intent(in) :: vf(nlay, nq)
    real(dp), intent(inout) :: q(nlay, nq)
    real(dp), intent(out), optional :: mass_out(nq)
    real(dp), intent(out), optional :: mass_clipped(nq)
    real(dp), intent(out), optional :: courant_max
    integer, intent(out), optional :: ierr

    real(dp) :: pe(nlev), pl(nlay), rho(nlay), dm(nlay)
    real(dp) :: medge(nlev), mmid(nlay), mdep(nlev)
    real(dp) :: wcell(nlay), wface(nlev)
    real(dp) :: slope(nlay), prefix(nlev), q_new(nlay)
    real(dp) :: cell_mass, old_mass, new_mass, grav, mass_tol
    real(dp) :: clipped_mass, tracer_courant
    integer :: status
    integer :: i, k, n

    status = CSLM_OK
    if (present(mass_out)) mass_out = 0.0_dp
    if (present(mass_clipped)) mass_clipped = 0.0_dp
    if (present(courant_max)) courant_max = 0.0_dp
    if (present(ierr)) ierr = CSLM_OK

    if (nlev /= nlay + 1 .or. nlay < 1 .or. nq < 1) then
      if (present(ierr)) then
        ierr = CSLM_ERR_DIMENSIONS
        return
      end if
      error stop 'vert_adv_csl_mass: invalid dimensions'
    end if
    if (.not. ieee_is_finite(t_end) .or. t_end < 0.0_dp) then
      if (present(ierr)) then
        ierr = CSLM_ERR_TIMESTEP
        return
      end if
      error stop 'vert_adv_csl_mass: t_end must be non-negative'
    end if
    if (any(.not. ieee_is_finite(mu)) .or. any(.not. ieee_is_finite(Tl)) .or. &
        any(.not. ieee_is_finite(pl_in)) .or. any(.not. ieee_is_finite(pe_in)) .or. &
        .not. ieee_is_finite(grav_in) .or. &
        any(mu <= 0.0_dp) .or. any(Tl <= 0.0_dp) .or. any(pl_in <= 0.0_dp) .or. &
        any(pe_in <= 0.0_dp) .or. grav_in <= 0.0_dp .or. &
        any(pe_in(2:nlev) <= pe_in(1:nlay))) then
      if (present(ierr)) then
        ierr = CSLM_ERR_GRID
        return
      end if
      error stop 'vert_adv_csl_mass: invalid thermodynamic state or pressure grid'
    end if
    if (t_end <= 0.0_dp) return

    grav = 100.0_dp * grav_in       ! m s-2 -> cm s-2
    pe = 10.0_dp * pe_in            ! Pa -> dyn cm-2
    pl = 10.0_dp * pl_in            ! Pa -> dyn cm-2

    ! Downward hydrostatic column-mass coordinate.  dm has units g cm-2.
    medge(1) = 0.0_dp
    do i = 1, nlay
      dm(i) = (pe(i+1) - pe(i))/grav
      medge(i+1) = medge(i) + dm(i)
      mmid(i) = 0.5_dp*(medge(i) + medge(i+1))
      rho(i) = pl(i)/((R_GAS/mu(i))*Tl(i))
      if (.not. ieee_is_finite(dm(i)) .or. dm(i) <= 0.0_dp .or. &
          .not. ieee_is_finite(rho(i)) .or. rho(i) <= 0.0_dp) then
        if (present(ierr)) then
          ierr = CSLM_ERR_GRID
          return
        end if
        error stop 'vert_adv_csl_mass: invalid layer mass or density'
      end if
    end do

    do n = 1, nq
      if (any(.not. ieee_is_finite(q(:,n))) .or. &
          any(.not. ieee_is_finite(vf(:,n)))) then
        if (present(ierr)) then
          ierr = CSLM_ERR_NONFINITE
          return
        end if
        error stop 'vert_adv_csl_mass: non-finite tracer or velocity'
      end if

      clipped_mass = sum(max(-q(:,n), 0.0_dp)*dm)
      if (clipped_mass > 0.0_dp) then
        q(:,n) = max(q(:,n), 0.0_dp)
        status = ior(status, CSLM_WARN_Q_CLIPPED)
      end if
      if (present(mass_clipped)) mass_clipped(n) = clipped_mass

      if (any(vf(:,n) < 0.0_dp)) status = ior(status, CSLM_WARN_V_CLIPPED)

      ! Interpolate W=rho*v itself.  This is important for Epstein settling:
      ! when v is proportional to 1/rho, a constant W remains exactly
      ! constant at every face.
      wcell = rho*max(0.0_dp, vf(:,n))
      wface(1) = wcell(1)
      do k = 2, nlay
        wface(k) = 0.5_dp*(wcell(k-1) + wcell(k))
      end do
      wface(nlev) = wcell(nlay)

      tracer_courant = 0.0_dp
      do i = 1, nlay
        tracer_courant = max(tracer_courant, &
                             max(wface(i), wface(i+1))*t_end/dm(i))
      end do
      if (present(courant_max)) courant_max = max(courant_max, tracer_courant)

      do k = 1, nlev
        mdep(k) = departure_point(k, t_end, nlay, medge, wface)
      end do
      if (any(mdep(2:nlev) < mdep(1:nlay))) then
        if (present(ierr)) then
          ierr = CSLM_ERR_TRAJECTORY
          return
        end if
        error stop 'vert_adv_csl_mass: characteristic trajectories crossed'
      end if

      ! In mass coordinates, q itself is the conserved-density field because
      ! tracer column mass is the integral of q dm.
      call limited_slopes(nlay, medge, mmid, q(:,n), slope)

      prefix(1) = 0.0_dp
      do i = 1, nlay
        prefix(i+1) = prefix(i) + dm(i)*q(i,n)
      end do
      old_mass = prefix(nlev)
      mass_tol = 100.0_dp*epsilon(1.0_dp)*max(old_mass, tiny(1.0_dp))

      new_mass = 0.0_dp
      do i = 1, nlay
        cell_mass = integrated_mass(mdep(i+1), nlay, medge, mmid, q(:,n), slope, prefix) - &
                    integrated_mass(mdep(i),   nlay, medge, mmid, q(:,n), slope, prefix)

        if (cell_mass < -mass_tol) then
          if (present(ierr)) then
            ierr = CSLM_ERR_REMAP
            return
          end if
          error stop 'vert_adv_csl_mass: negative remapped cell mass'
        end if
        ! The reconstruction and ordered departure points are positive by
        ! construction.  max only removes round-off-sized negative values.
        cell_mass = max(0.0_dp, cell_mass)
        q_new(i) = cell_mass/dm(i)
        new_mass = new_mass + cell_mass
      end do

      if (new_mass > old_mass + mass_tol) then
        if (present(ierr)) then
          ierr = CSLM_ERR_CONSERVATION
          return
        end if
        error stop 'vert_adv_csl_mass: total tracer mass increased'
      end if
      q(:,n) = q_new
      if (present(mass_out)) mass_out(n) = max(0.0_dp, old_mass - new_mass)
    end do

    if (present(ierr)) ierr = status

  end subroutine vert_adv_csl_mass


  function departure_point(iface, dt, nlay, xedge, wface) result(xdep)
    !! Trace one face backwards exactly through a continuous piecewise-linear
    !! and time-independent velocity profile.
    implicit none

    integer, intent(in) :: iface, nlay
    real(dp), intent(in) :: dt, xedge(nlay+1), wface(nlay+1)
    real(dp) :: xdep

    real(dp) :: pos, remaining, xleft, width
    real(dp) :: wleft, wpos, gradw, cross_time, wdep, wscale
    integer :: cell

    pos = xedge(iface)
    remaining = dt

    ! Above the model top, retain the top-face velocity.  The mass integral
    ! is zero there, which imposes zero tracer inflow.
    if (iface == 1) then
      xdep = pos - wface(1)*remaining
      return
    end if

    cell = iface - 1
    do while (remaining > 0.0_dp .and. cell >= 1)
      xleft = xedge(cell)
      width = xedge(cell+1) - xleft
      wleft = wface(cell)
      gradw = (wface(cell+1) - wleft)/width
      wpos = wleft + gradw*(pos - xleft)

      if (wpos <= TINY_W) then
        xdep = pos
        return
      end if

      wscale = max(abs(wpos), abs(wleft), tiny(1.0_dp))
      if (abs(gradw)*(pos - xleft) <= 100.0_dp*epsilon(1.0_dp)*wscale) then
        cross_time = (pos - xleft)/wpos
        if (remaining < cross_time) then
          xdep = pos - wpos*remaining
          return
        end if
      else if (wleft <= TINY_W) then
        ! A zero-velocity face is approached only asymptotically.
        wdep = wpos*exp(-gradw*remaining)
        xdep = xleft + (wdep - wleft)/gradw
        return
      else
        cross_time = log(wpos/wleft)/gradw
        if (remaining < cross_time) then
          wdep = wpos*exp(-gradw*remaining)
          xdep = xleft + (wdep - wleft)/gradw
          return
        end if
      end if

      remaining = max(0.0_dp, remaining - cross_time)
      pos = xleft
      cell = cell - 1
    end do

    xdep = pos - wface(1)*remaining
  end function departure_point


  subroutine limited_slopes(nlay, xedge, xmid, conc, slope)
    implicit none

    integer, intent(in) :: nlay
    real(dp), intent(in) :: xedge(nlay+1), xmid(nlay), conc(nlay)
    real(dp), intent(out) :: slope(nlay)

    real(dp) :: dl, dc, dr, mc_slope
    real(dp) :: half_width, cmin, cmax, slope_bound
    integer :: i

    slope = 0.0_dp
    do i = 2, nlay-1
      dl = (conc(i) - conc(i-1))/(xmid(i) - xmid(i-1))
      dc = (conc(i+1) - conc(i-1))/(xmid(i+1) - xmid(i-1))
      dr = (conc(i+1) - conc(i))/(xmid(i+1) - xmid(i))
      mc_slope = minmod3(2.0_dp*dl, dc, 2.0_dp*dr)

      ! The centre-based MC limiter alone is not positivity preserving on a
      ! stretched grid.  Cap it so both reconstructed edge values remain
      ! within the extrema of the neighbouring cell means.
      cmin = min(conc(i-1), conc(i), conc(i+1))
      cmax = max(conc(i-1), conc(i), conc(i+1))
      half_width = 0.5_dp*(xedge(i+1) - xedge(i))
      slope_bound = min(cmax - conc(i), conc(i) - cmin)/half_width
      slope(i) = sign(min(abs(mc_slope), max(0.0_dp, slope_bound)), mc_slope)
    end do
  end subroutine limited_slopes


  pure function minmod3(a, b, c) result(mm)
    implicit none

    real(dp), intent(in) :: a, b, c
    real(dp) :: mm

    if (a > 0.0_dp .and. b > 0.0_dp .and. c > 0.0_dp) then
      mm = min(a, b, c)
    else if (a < 0.0_dp .and. b < 0.0_dp .and. c < 0.0_dp) then
      mm = max(a, b, c)
    else
      mm = 0.0_dp
    end if
  end function minmod3


  function integrated_mass(x, nlay, xedge, xmid, conc, slope, prefix) result(mass)
    !! Integral of the reconstructed concentration from the model top to x.
    implicit none

    integer, intent(in) :: nlay
    real(dp), intent(in) :: x, xedge(nlay+1), xmid(nlay)
    real(dp), intent(in) :: conc(nlay), slope(nlay), prefix(nlay+1)
    real(dp) :: mass

    real(dp) :: offset, left_offset, local_offset
    integer :: lo, hi, mid, i

    if (x <= xedge(1)) then
      mass = 0.0_dp
      return
    else if (x >= xedge(nlay+1)) then
      mass = prefix(nlay+1)
      return
    end if

    lo = 1
    hi = nlay
    do while (lo < hi)
      mid = (lo + hi)/2
      if (x > xedge(mid+1)) then
        lo = mid + 1
      else
        hi = mid
      end if
    end do
    i = lo

    offset = x - xedge(i)
    left_offset = xedge(i) - xmid(i)
    local_offset = x - xmid(i)
    mass = prefix(i) + conc(i)*offset + &
           0.5_dp*slope(i)*(local_offset**2 - left_offset**2)

    ! Protect cumulative monotonicity from floating-point round-off.
    mass = min(prefix(i+1), max(prefix(i), mass))
  end function integrated_mass

end module vert_adv_csl_mass_mod
