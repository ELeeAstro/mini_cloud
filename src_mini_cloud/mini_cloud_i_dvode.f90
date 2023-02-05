module mini_cloud_i_dvode
  use mini_cloud_precision
  use mini_cloud_class
  use mini_cloud_saturation
  use mini_cloud_nucleation
  use mini_cloud_chem
  use mini_cloud_rhs
  implicit none

  logical :: first_call = .True.
  !$omp threadprivate (first_call)

  public ::  mini_cloud_dvode, RHS_update, jac_dummy

contains

  subroutine mini_cloud_dvode(n_dust, T_in, P_in, t_end, sp, k, k3, VMR, VMR0)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T_in, P_in, t_end
    character(len=11), dimension(n_dust), intent(in) :: sp

    real(dp), dimension(4), intent(inout) :: k
    real(dp), dimension(n_dust), intent(inout) :: k3
    real(dp), dimension(n_dust), intent(inout) :: VMR, VMR0

    integer :: ncall

    ! Time controls
    real(dp) :: t_begin, t_now

    ! dvode variables
    integer :: neq 

    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork, rtol, atol
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim
    real(dp) :: rpar
    integer :: ipar

    integer :: n
    real(dp) :: a_check, Ar_check, V_check, total_k3s

    if (first_call .eqv. .True.) then
      call mini_cloud_init(n_dust, sp)
      first_call = .False.
    end if

    neq = 4 + n_dust + n_dust
    allocate(y(neq))

    ! Prepare all species for the solver
    k(:) = max(k(:),1e-30_dp)
    k3(:) = max(k3(:),1e-30_dp)
    VMR(:) = max(VMR(:),1e-30_dp)
    VMR(:) = min(VMR0(:),VMR(:))

    d_sp(:)%VMR0 =  VMR0(:)

    ! Work variables
    T = T_in
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T)  ! Find initial number density [cm-3] of atmosphere

    ! Check if species are at seed particle limits 
    a_check = k(2)/k(1)
    Ar_check = fourpi * k(3)/k(1)
    V_check = fourpi3 * k(4)/k(1)

    if (k(1) > 1e-30_dp .and. &
      & (a_check <= a_seed*1.0_dp .or. Ar_check <= Ar_seed*1.0_dp .or. V_check <= V_seed*1.0_dp)) then
      k(1) = k(1)
      k(2) = a_seed * k(1)
      k(3) = Ar_seed/fourpi * k(1)
      k(4) = V_seed/fourpi3 * k(1)

      k3(1) = k(4)
      k3(2:n_dust) = 1e-30_dp
    end if



    call calc_saturation(n_dust, VMR(:))
    call calc_nucleation(n_dust, VMR(:))
    ! ! If small number of dust particles and small nucleation rate - don't integrate
    if ((sum(d_sp(:)%Js) < 1.0e-10_dp) .and. (k(1) < 1.0e-10_dp)) then
      return
    end if


    call calc_chem(n_dust, k(:), k3(:), VMR(:))
    ! If overall growth is super slow, don't bother integrating
     if (k(1) > 1.0e-10_dp .and. (sum(d_sp(:)%Js) < 1.0e-10_dp &
       & .and. abs(sum(d_sp(:)%chis)) < 1.0e-12_dp)) then
         return
    end if

    ! -----------------------------------------
    ! ***  parameters for the dvode-solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    ! Problem is stiff (usual)
    mf = 22
    rworkdim = 22 +  9*neq + 2*neq**2
    iworkdim = 30 + neq
    allocate(rtol(neq), atol(neq), rwork(rworkdim),iwork(iworkdim))

    itol = 4
    rtol(:) = 1e-10_dp           ! Relative tolerances for each scalar
    atol(:) = 1e-30_dp               ! Absolute tolerance for each scalar (floor value)

    rwork(:) = 0.0_dp
    iwork(:) = 0

    rwork(1) = 0.0_dp               ! Critical T value (don't integrate past time here)
    rwork(5) = 0.0_dp              ! Initial starting timestep (start low, will adapt in DVODE)
    rwork(6) = 1.0e-1_dp       ! Maximum timestep (for heavy evaporation ~0.1 is required)

    iwork(5) = 5               ! Max order required
    iwork(6) = 100000               ! Max number of internal steps
    iwork(7) = 1                ! Number of error messages 


    ! Prepare all species for the solver
    y(1:4) = max(k(:),1e-30_dp)
    y(5:5+n_dust-1) = max(k3(:),1e-30_dp)
    y(5+n_dust:5+n_dust+n_dust-1) = max(VMR(:),1e-30_dp)

      ! Rescale k3s again
      total_k3s = sum(y(5:5+n_dust-1))
      if (total_K3s < 1.0e-20_dp) then
        y(5) = y(4)
        y(6:5+n_dust-1) = 1.0e-30_dp
      else
        y(5:5+n_dust-1) = (y(5:5+n_dust-1) / total_k3s) * y(4)
      end if

    rpar = 0.0_dp
    ipar = n_dust

    t_begin = 0.0_dp
    t_now = t_begin

    ncall = 1

    do while((t_now < t_end) .and. ncall <= 2)

      y(1:5+n_dust-1) = y(1:5+n_dust-1)/nd_atm

      call DVODE (RHS_update, neq, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dummy, mf, rpar, ipar)

      y(1:5+n_dust-1) = y(1:5+n_dust-1) * nd_atm

       do n = 1, ipar
        y(5+ipar-1+n) = min(VMR0(n),y(5+ipar-1+n))
      end do

      ncall = ncall + 1

      ! Check if species are at seed particle limits 
      a_check = y(2)/y(1)
      Ar_check = fourpi * y(3)/y(1)
      V_check = fourpi3 * y(4)/y(1)

      if (y(1) > 1e-30_dp .and. &
        & (a_check <= a_seed*1.0_dp .or. Ar_check <= Ar_seed*1.0_dp .or. V_check <= V_seed*1.0_dp)) then
        y(1) = y(1)
        y(2) = a_seed * y(1)
        y(3) = Ar_seed/fourpi * y(1)
        y(4) = V_seed/fourpi3 * y(1)

        y(5) = y(4)
        y(6:5+n_dust-1) = 1e-30_dp
      end if

      if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        ! Some critical failure - comment out if don't care
        print*, istate, real(t_now), real(rwork(11)), int(T), real(P_cgs/1e6_dp)
        exit
      end if

      ! Rescale k3s again
      total_k3s = sum(y(5:5+n_dust-1))
      if (total_K3s < 1.0e-20_dp) then
        y(5) = y(4)
        y(6:5+n_dust-1) = 1.0e-30_dp
      else
        y(5:5+n_dust-1) = (y(5:5+n_dust-1) / total_k3s) * y(4)
      end if

    end do

          ! Rescale k3s again
      total_k3s = sum(y(5:5+n_dust-1))
      if (total_K3s < 1.0e-20_dp) then
        y(5) = y(4)
        y(6:5+n_dust-1) = 1.0e-30_dp
      else
        y(5:5+n_dust-1) = (y(5:5+n_dust-1) / total_k3s) * y(4)
      end if

    k(:) = max(y(1:4),1e-30_dp)
    k3(:) = max(y(5:5+n_dust-1),1e-30_dp)
    VMR(:) = max(y(5+n_dust:5+n_dust+n_dust-1),1e-30_dp)

    VMR(:) = min(VMR0(:),VMR(:))

    deallocate(rtol, atol, rwork, iwork, y)

  end subroutine mini_cloud_dvode

  subroutine RHS_update(neq, time, y, f, rpar, ipar)
    implicit none

    integer, intent(in) :: neq
    real(dp), intent(inout) :: time
    real(dp), dimension(neq), intent(inout) :: y
    real(dp), dimension(neq), intent(inout) :: f
    real(dp), intent(inout) :: rpar
    integer, intent(inout) :: ipar

    integer :: n
    real(dp) :: tot_k3s

    y(1:5+ipar-1) = y(1:5+ipar-1) * nd_atm

    !y(5:5+ipar-1) = max(y(5:5+ipar-1),1e-30_dp)
    y(:) = max(y(:),1.0e-30_dp)

    do n = 1, ipar
      y(5+ipar-1+n) = min(d_sp(n)%VMR0,y(5+ipar-1+n))
    end do

    ! Rescale k3s again
    tot_k3s = sum(y(5:5+ipar-1))
    if (tot_K3s < 1.0e-20_dp) then
      y(5) = y(4)
      y(6:5+ipar-1) = 1.0e-30_dp
    else
      y(5:5+ipar-1) = (y(5:5+ipar-1) / tot_k3s) * y(4)
    end if

    call calc_saturation(ipar, y(5+ipar:5+ipar+ipar-1))
    call calc_nucleation(ipar, y(5+ipar:5+ipar+ipar-1))
    call calc_chem(ipar, y(1:4), y(5:5+ipar-1), y(5+ipar:5+ipar+ipar-1))
    call calc_seed_evap(ipar, y(1:4))
    call form_RHS(ipar, y(1:4), f)

    y(1:5+ipar-1) = y(1:5+ipar-1)/nd_atm
    f(1:5+ipar-1) = f(1:5+ipar-1)/nd_atm

  end subroutine RHS_update

  pure subroutine jac_dummy (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: T
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD,NEQ), intent(inout) :: PD
    real(dp), intent(in) :: RPAR
    integer, intent(in) :: IPAR
  end subroutine jac_dummy

end module mini_cloud_i_dvode
