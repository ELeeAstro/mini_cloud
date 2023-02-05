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

  integer, parameter :: n_dust = 4
  integer :: neq, nn

  real(dp), dimension(4) :: LL
  real(dp), dimension(n_dust) :: L3s
  real(dp), dimension(n_dust) :: VMR

  public ::  mini_cloud_dvode, RHS_update, jac_dummy
  private :: n_dust, neq, nn

contains

  subroutine mini_cloud_dvode(T_in, P_in, mu_in, t_end, sp)
    implicit none

    real(dp), intent(in) :: T_in, P_in, mu_in, t_end
    character(len=11), dimension(n_dust), intent(in) :: sp

    integer :: ncall

    ! Time controls
    real(dp) :: t_begin, t_now

    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork, rtol, atol
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim
    real(dp) :: rpar
    integer :: ipar

    integer :: n
    real(dp) :: a_check, Ar_check, V_check
    real(dp) :: tot_L3s
    real(dp), dimension(n_dust) :: nd

    if (first_call .eqv. .True.) then
      call mini_cloud_init(n_dust, sp)
      first_call = .False.
    end if

    neq = 4 + n_dust + n_dust
    nn = 4 + n_dust
    allocate(y(neq))

    ! Work variables
    T = T_in
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T)  ! Find initial number density [cm-3] of atmosphere
    rho = (P_cgs* mu_in * amu)/(kb*T)
    nd(:) = VMR(:) * nd_atm

    a_check = LL(2)/LL(1)  * vol2rad
    Ar_check = LL(3)/LL(1) * Afak
    V_check = LL(4)/LL(1)

    if (LL(1)*rho > 1e-10_dp .and. &
      & (a_check <= a_seed*1.05_dp .or. Ar_check <= Ar_seed*1.05_dp .or. V_check <= V_seed*1.05_dp)) then
      LL(1) = LL(1)
      LL(2) = (a_seed * LL(1)) / vol2rad
      LL(3) = (Ar_seed * LL(1)) / Afak
      LL(4) = V_seed * LL(1)

      L3s(1) = LL(4)
      L3s(2:n_dust) = L_Lim
    end if

    call calc_saturation(n_dust, nd(:))
    call calc_nucleation(n_dust, nd(:))
    ! ! If small number of dust particles and small nucleation rate - don't integrate
    if ((sum(d_sp(:)%Js) < 1.0e-10_dp) .and. (LL(1)*rho < 1.0e-10_dp)) then
      return
    end if


    call calc_chem(n_dust, LL(:), L3s(:), nd(:))
    ! If overall growth is super slow, don't bother integrating
     if (LL(1)*rho > 1.0e-10_dp .and. (sum(d_sp(:)%Js) < 1.0e-10_dp &
       & .and. abs(sum(d_sp(:)%chis)) < 1.0e-20_dp)) then
         return
    end if

    ! -----------------------------------------
    ! ***  parameters for the dvode-solver  ***
    ! -----------------------------------------

    itask = 4
    istate = 1
    iopt = 1

    ! Problem is stiff (usual)
    mf = 22
    rworkdim = 22 +  9*neq + 2*neq**2
    iworkdim = 30 + neq
    allocate(rtol(neq), atol(neq), rwork(rworkdim),iwork(iworkdim))

    itol = 4
    rtol(:) = 1e-12_dp           ! Relative tolerances for each scalar
    atol(:) = 1e-30_dp               ! Absolute tolerance for each scalar (floor value)

    rwork(:) = 0.0_dp
    iwork(:) = 0

    rwork(1) = t_end               ! Critical T value (don't integrate past time here)
    rwork(5) = 1.0e-30_dp              ! Initial starting timestep (start low, will adapt in DVODE)
    rwork(6) = 1.0e-2_dp       ! Maximum timestep (for heavy evaporation ~0.1 is required)

    iwork(5) = 5               ! Max order required
    iwork(6) = 100000               ! Max number of internal steps
    iwork(7) = 1                ! Number of error messages 


    ! Prepare all species for the solver
    y(1:4) = max(LL(:),1e-30_dp)
    y(5:nn) = max(L3s(:),1e-30_dp)
    y(nn+1:neq) = max(VMR(:),1e-30_dp)

    rpar = 0.0_dp
    ipar = n_dust

    t_begin = 0.0_dp
    t_now = t_begin

    ncall = 1

    do while((t_now < t_end) .and. ncall <= 3)

      LL(:) = max(LL(:),1e-30_dp)
      L3s(:) = max(L3s(:),1e-30_dp)
      VMR(:) = max(VMR(:),1e-30_dp)

      y(1:4) = max(LL(:),1e-30_dp)
      y(5:nn) = max(L3s(:),1e-30_dp)
      y(nn+1:neq) = max(VMR(:),1e-30_dp)  
      
      call DVODE (RHS_update, neq, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dummy, mf, rpar, ipar)


      LL(:) = y(1:4)
      L3s(:) = y(5:nn)
      VMR(:) = y(nn+1:neq)

      LL(:) = max(LL(:),1e-30_dp)
      L3s(:) = max(L3s(:),1e-30_dp)
      VMR(:) = max(VMR(:),1e-30_dp)

      nd(:) = VMR(:) * nd_atm

      ncall = ncall + 1

      call calc_saturation(n_dust, nd(:))
      call calc_nucleation(n_dust, nd(:))
      call calc_chem(n_dust, LL(:), L3s(:), nd(:))
      call calc_seed_evap(n_dust, LL(:))

       a_check = LL(2)/LL(1) * vol2rad
       Ar_check = LL(3)/LL(1) * Afak
       V_check = LL(4)/LL(1)

       if (LL(1)*rho > 1e-10_dp .and. &
         & (a_check <= a_seed*1.05_dp .or. Ar_check <= Ar_seed**1.05_dp .or. V_check <= V_seed**1.05_dp)) then
         LL(1) = LL(1)
         LL(2) = (a_seed * LL(1)) / vol2rad
         LL(3) = (Ar_seed * LL(1)) / Afak
         LL(4) = V_seed * LL(1)

         L3s(1) = LL(4)
         L3s(2:n_dust) = L_Lim

         if (d_sp(1)%sevap < 0.0_dp) then
           print*, 'Seed evap'
           LL(:) = L_lim
           L3s(:) = L_lim
           t_now = t_end
           cycle
         end if
        end if


      if (istate == -1) then
        istate = 2
        !rwork(5) = 1.0e-30_dp
      else if (istate < -1) then
        ! Some critical failure - comment out if don't care
        print*, istate, real(t_now), real(rwork(11)), int(T), real(P_cgs/1e6_dp)
        exit
      end if
  
    end do

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
    real(dp) :: tot_L3s
    real(dp), dimension(ipar) :: nd 

    ! Solution vector from DVODE - return to local variables
    LL(:) = y(1:4)
    L3s(:) = y(5:nn)
    VMR(:) = y(nn+1:neq)

    ! Limiters here
    LL(:) = max(LL(:), L_lim)
    L3s(:) = max(L3s(:), L_lim)
    VMR(:) = max(VMR(:), L_lim)

    nd(:) = VMR(:) * nd_atm

    ! Rescale L3s
    tot_L3s = sum(L3s(:))
    if (tot_L3s < 1.0e-20_dp) then
       L3s(1) = LL(4)
       L3s(2:) = L_lim
    else
      L3s(:) = (L3s(:) / tot_L3s) * LL(4)
    end if

    call calc_saturation(ipar, nd(:))
    call calc_nucleation(ipar, nd(:))
    call calc_chem(ipar, LL(:), L3s(:), nd(:))
    call calc_seed_evap(ipar, LL(:))
    call form_RHS(ipar, LL(:), f)

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
