module mini_cloud_settle
  use mini_cloud_precision
  implicit none

  public :: adv_tracer_mccormack
  private :: adv_prepare, minmod

contains

  subroutine adv_tracer_mccormack(nlay, nq, tend, q, vf, h)
    implicit none

    integer, intent(in) :: nlay, nq
    real(dp), dimension(nq,nlay), intent(inout) :: q
    real(dp), dimension(nlay), intent(in) :: vf 
    real(dp), dimension(nlay+1), intent(in) :: h
    real(dp), intent(in) :: tend

    real(dp), parameter :: cfl = 0.90_dp

    real(dp), dimension(nq,nlay) :: qc
    real(dp), dimension(nlay) :: sig, c, vfc
    real(dp) :: D, tnow, dt
    integer :: i, iit, n

    vfc(:) = vf(:) / 100.0_dp

    call adv_prepare(nlay+1, h, dh)

    dt = 9.0e99_dp
    do i = 1, nlay
      D  = abs(vfc(i))
      if (D <= 0.0_dp) then
        cycle
      end if
      dt = MIN(dt,cfl*(dh(i))/D)
    enddo
    dt = dt * (1.0_dp + 1.e-12_dp)

    tnow = 0.0_dp
    iit = 1

    do while ((tnow < tend) .and. (iit < 1000000))

      ! If next time step overshoots - last time step is equal tend
      if (tnow + dt > tend) then
        dt = tend - tnow
      end if

      !! Find the courant number
      c(:) = abs(vfc(:)) * dt / dh(:)

      do n = 1, nq

        !! Find the minmod limiter
        call minmod(nlay,q(n,:),dh,sig)

        !! Perform McCormack step
        qc(n,:) = q(n,:)
        qc(n,:nlay-1) = q(n,:nlay-1) - sig(:nlay-1)*c(:nlay-1)*(q(n,2:nlay) - q(n,:nlay-1))
        q(n,2:nlay) = 0.5_dp * (q(n,2:nlay) + qc(n,2:nlay) - c(2:nlay)*(qc(n,2:nlay) - qc(n,:nlay-1)))

        q(n,:) = max(q(n,:),1e-30_dp)

      end do

      tnow = tnow + dt
      iit = iit + 1

      q(:,nlay) = 1e-30_dp

    end do

  end subroutine adv_tracer_mccormack   

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

  subroutine adv_prepare(nlev, h, dh)
    implicit none

    integer, intent(in) :: nlev
    real(dp), dimension(nlev), intent(in) :: h

    real(dp), dimension(nlev-1), intent(out) :: dh
    
    real(dp), dimension(nlev-1) :: hmid
    integer :: i

    ! Find distance between each layer center
    hmid(:) = (h(1:nlev-1) +  h(2:nlev))/2.0_dp
    do i = nlay-1,1,-1
      dh(i) =  hmid(i) - hmid(i+1)
    end do
    dh(nlay) = hmid(nlay)

  end subroutine adv_prepare

end module mini_cloud_settle