module vert_diff_exp_2_mod
  use, intrinsic :: iso_fortran_env
  implicit none


  integer, parameter :: dp = REAL64 ! Precision variable

  real(dp), parameter :: CFL = 0.45_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: kb = 1.380649e-16_dp

  public :: vert_diff_exp_2

contains

  subroutine vert_diff_exp_2(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, Kzz_in, nq, q_in, q0)
    implicit none

    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), dimension(nlay), intent(in) :: Tl, pl_in, Kzz_in, mu
    real(dp), dimension(nlev), intent(in) :: pe_in
    real(dp), dimension(nq), intent(in) :: q0

    real(dp), dimension(nlay,nq), intent(inout) :: q_in
    
    real(dp), dimension(nlay,nq) :: q

    integer :: k
    real(dp) :: grav
    real(dp), dimension(nlev) :: alte_r, pe
    real(dp), dimension(nlay) :: nd_r, nd, Kzz, alt, altm
    real(dp), dimension(nlay) :: pl
    real(dp), dimension(2:nlay-1) :: hil, hir, dil1, dim1, dir1, dil2, dim2, dir2

    real(dp) :: h2, h3, hm1, hm2, d1l, d1m, d1r, dnl, dnm, dnr

    real(dp), dimension(2:nlay-1,nq) :: flx
    real(dp), dimension(nq) :: flx1, flxn


    integer :: n_it, n, accept, ierr
    real(dp) :: dt, t_now, dt_max

    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2

    grav = grav_in * 100.0_dp


    !! First calculate the vertical height (cm) assuming hydrostatic equilibrium and differences
    alte_r(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte_r(k) = alte_r(k+1) + (R*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do

    !! Get midpoints of altitudes
    do k = 1, nlay
      altm(k) = (alte_r(k) + alte_r(k+1))/2.0_dp
    end do

    !! Number density of atmosphere at layers
    nd_r(:) = pl(:)/(kb*Tl(:))  

    !! Reverse required arrays for ascending indexes
    do k = 1, nlay
      alt(k) = altm(nlay-k+1)
      nd(k) = nd_r(nlay-k+1)
      Kzz(k) = Kzz_in(nlay-k+1)
      q(k,:) = q_in(nlay-k+1,:)
    end do

    !! We now have all variables at the cell centers and can carry out differencing equations
    !! Calculate h and d values
    do k = 2, nlay-1
      hil(k) = (alt(k) - alt(k-1))
      hir(k) = (alt(k+1) - alt(k))

      dil1(k) = - (hir(k)/((hir(k) + hil(k))*hil(k)))
      dim1(k) = + ((hir(k) - hil(k))/(hil(k)*hir(k)))
      dir1(k) = + (hil(k)/((hir(k) + hil(k))*hir(k)))
      dil2(k) = + 2.0_dp/((hir(k) + hil(k))*hil(k))
      dim2(k) = - 2.0_dp/(hir(k)*hil(k))
      dir2(k) = + 2.0_dp/((hir(k) + hil(k))*hir(k))
    end do

    h2 = alt(2) - alt(1) ; h3 = alt(3) - alt(1)
    hm1 = alt(nlay) - alt(nlay-1) ; hm2 = alt(nlay) - alt(nlay-2)

    d1l =  -(h2-h3)/(h2*h3) ; d1m = h3/(h2*(h3 - h2));  d1r = -h2/(h3*(h3 - h2))
    dnl =  (hm1 + hm2)/(hm1*hm2); dnm = -hm2/(hm1*(hm2 - hm1)) ; dnr = hm1/(hm2*(hm2 - hm1))


    !! Prepare timestepping routine
    !! Find minimum timestep that satifies the CFL condition
    dt_max = t_end
    do k = 2, nlay
      dt_max = min(dt_max,CFL*((alt(k) - alt(k-1))**2/(0.5_dp*(Kzz(k) + Kzz(k-1)))))
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

      do k = 2, nlay-1
        flx(k,:) = (dil1(k)*nd(k-1)*Kzz(k-1) + dim1(k)*nd(k)*Kzz(k) + dir1(k)*nd(k+1)*Kzz(k+1)) & 
          & * (dil1(k)*q(k-1,:) + dim1(k)*q(k,:) + dir1(k)*q(k+1,:)) & 
          & + nd(k)*Kzz(k)*(dil2(k)*q(k-1,:) + dim2(k)*q(k,:) + dir2(k)*q(k+1,:))
      end do

      flx1(:) = 0.0_dp!-nd(1)*Kzz(1)*(d1l*q(1,:) + d1m*q(2,:) + d1r*q(3,:))
      flxn(:) = 0.0_dp

      do k = 2, nlay-1
        q(k,:) = q(k,:) + flx(k,:)/nd(k) * dt
      end do

      do n = 1, nq
        q(:,n) = max(q(:,n),1e-30_dp)
      end do

      q(1,:) = q0(:)
      q(nlay,:) = q(nlay-1,:)

      t_now = t_now + dt
      n_it = n_it + 1

    end do

    !! Reverse required arrays for output
    do k = 1, nlay
      q_in(k,:) = q(nlay-k+1,:)
    end do

  end subroutine vert_diff_exp_2

end module vert_diff_exp_2_mod

 