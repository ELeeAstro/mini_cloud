module mini_dihrt_i_ros4
  use mini_dihrt_precision
  use mini_dihrt_class
  use mini_dihrt_saturation
  use mini_dihrt_nucleation
  use mini_dihrt_chem
  use mini_dihrt_rhs
  implicit none

  logical :: first_call = .True.
  integer :: ipar


  public ::  mini_dihrt_ros4, RHS_update, jac_dummy, mas_dummy, solout

contains

  subroutine mini_dihrt_ros4(n_dust, T_in, P_in, t_end, sp, k, k3, VMR)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T_in, P_in, t_end
    character(len=11), dimension(n_dust), intent(in) :: sp

    real(dp), dimension(4), intent(inout) :: k
    real(dp), dimension(n_dust), intent(inout) :: k3
    real(dp), dimension(n_dust), intent(inout) :: VMR

    integer :: n
    integer :: ncall

    ! Time controls
    real(dp) :: t_begin, t_now, dt_init, t_old

    ! ros4 variables
    integer :: neq 
    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    real(dp) :: rtol, atol
    integer :: itol, ijac, mljac, mujac, imas, mlmas, mumas, iout, lrwork, liwork, idid
    integer :: ifcn, idfx

    real(dp) :: a_check, Ar_check, V_check

    if (first_call .eqv. .True.) then
      call mini_dihrt_init(n_dust, sp)
      first_call = .False.
    end if

    neq = 4 + n_dust + n_dust
    allocate(y(neq))

    ! Work variables
    T = T_in
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T)  ! Find initial number density [cm-3] of atmosphere

    ! Check if species are at seed particle limits 
    a_check = k(2)/k(1)
    Ar_check = fourpi * k(3)/k(1)
    V_check = fourpi3 * k(4)/k(1)

    if (k(1) > 1e-10_dp .and. &
      & (a_check <= a_seed*1.05_dp .or. Ar_check <= Ar_seed*1.05_dp .or. V_check <= V_seed*1.05_dp)) then
      k(1) = k(1)
      k(2) = a_seed * k(1)
      k(3) = Ar_seed/fourpi * k(1)
      k(4) = V_seed/fourpi3 * k(1)

      k3(1) = k(4)
      k3(2:n_dust) = 1e-30_dp
    end if


    call calc_saturation(n_dust, VMR(:))
    call calc_nucleation(n_dust, VMR(:))

    !print*, d_sp(:)%sat, d_sp(:)%Js, k(1)

    ! ! If small number of dust particles and small nucleation rate - don't integrate
    if ((sum(d_sp(:)%Js) < 1.0e-10_dp) .and. (k(1) < 1.0e-10_dp)) then
      return
    end if


    call calc_chem(n_dust, k(:), k3(:), VMR(:))
    ! If overall growth is super slow, don't bother integrating
     if (k(1) > 1.0e-10_dp .and. (sum(d_sp(:)%Js) < 1.0e-10_dp &
       & .and. abs(sum(d_sp(:)%chis)) < 1.0e-20_dp)) then
         return
    end if

    ! -----------------------------------------
    ! ***  parameters for the ros4-solver  ***
    ! -----------------------------------------

    ifcn = 0
    idfx = 0
    rtol = 1.0e-4_dp
    atol = 1.0e-30_dp
    itol = 0
    ijac = 0
    mljac = neq
    mujac = 0
    imas = 0
    mlmas = neq
    mumas = 0
    iout = 0
    idid = 0

    ! Real work array
    lrwork = 2*neq*neq + 8*neq+5
    allocate(rwork(lrwork))
    rwork(:) = 0.0_dp

    rwork(1) = 0.0_dp!1.0e-16_dp ! Rounding unit - default 1e-16
    rwork(2) = t_end ! Max step size
    rwork(3) = 0.2_dp ! Parameter for step size selection - default 0.2
    rwork(4) = 6.0_dp ! Parameter for step size selection - deftaul 6.0
    rwork(5) = 0.1_dp ! Step rejection factor - default 0.1

    ! Integer work array
    liwork = neq + 2
    allocate(iwork(liwork))
    iwork(:) = 0

    iwork(1) = 0 ! Default step numbers - 0 = 100000
    iwork(2) = 1 ! Coefficent choice - default 1 (1,2,3)

    ! Prepare all species for the solver
    y(1:4) = max(k(:),1e-30_dp)
    y(5:5+n_dust-1) = max(k3(:),1e-30_dp)
    y(5+n_dust:5+n_dust+n_dust-1) = max(VMR(:),1e-30_dp)

    ipar = n_dust

    t_begin = 0.0_dp
    t_now = t_begin
    dt_init = 1.0e-99_dp

    ncall = 1

    do while((t_now < t_end))


      call  ROS4(neq,RHS_update,ifcn,t_now,y,t_end,dt_init, &
      &                  rtol,atol,itol, &
      &                  jac_dummy,ijac,mljac,mujac, &
      &                  dfx_dummy,idfx, &
      &                  mas_dummy,imas,mlmas,mumas, &
      &                  solout,iout, &
      &                  rwork,lrwork,iwork,liwork,idid)

      ncall = ncall + 1

      call calc_saturation(ipar, y(5+ipar:5+ipar+ipar-1))
      call calc_nucleation(ipar, y(5+ipar:5+ipar+ipar-1))
      call calc_chem(ipar, y(1:4), y(5:5+ipar-1), y(5+ipar:5+ipar+ipar-1))
      call calc_seed_evap(ipar, y(1:4))

      ! Check if species are at seed particle limits 
      a_check = y(2)/y(1)
      Ar_check = fourpi * y(3)/y(1)
      V_check = fourpi3 * y(4)/y(1)

      if (y(1) > 1e-10_dp .and. &
        & (a_check <= a_seed*1.05_dp .or. Ar_check <= Ar_seed*1.05_dp .or. V_check <= V_seed*1.05_dp)) then
        y(1) = y(1)
        y(2) = a_seed * y(1)
        y(3) = Ar_seed/fourpi * y(1)
        y(4) = V_seed/fourpi3 * y(1)

        y(5) = y(4)
        y(6:5+n_dust-1) = 1e-30_dp

        if (d_sp(1)%sevap < 0.0_dp) then
          print*, 'Seed evap'
          y(1:5+n_dust-1) = 1e-30_dp
          t_now = t_end
          cycle
        end if

      end if

      if (ncall > 10) then
        exit
      end if

    end do


    k(:) = max(y(1:4),1e-30_dp)
    k3(:) = max(y(5:5+n_dust-1),1e-30_dp)
    VMR(:) = max(y(5+n_dust:5+n_dust+n_dust-1),1e-30_dp)

    deallocate(rwork, iwork, y)

  end subroutine mini_dihrt_ros4

  subroutine RHS_update(neq, time, y, f)
    implicit none

    integer, intent(in) :: neq
    real(dp), intent(inout) :: time
    real(dp), dimension(neq), intent(inout) :: y
    real(dp), dimension(neq), intent(inout) :: f

    real(dp) :: tot_k3 

    ! Rescale L3s
    tot_k3 = sum(y(5:5+ipar-1))
    if (tot_k3 < 1.0e-20_dp) then
       y(5) = y(4)
       y(6:5+ipar-1) = 1e-30_dp
    else
      y(5:5+ipar-1) = (y(5:5+ipar-1) / tot_k3) * y(4)
    end if

    call calc_saturation(ipar, y(5+ipar:5+ipar+ipar-1))
    call calc_nucleation(ipar, y(5+ipar:5+ipar+ipar-1))
    call calc_chem(ipar, y(1:4), y(5:5+ipar-1), y(5+ipar:5+ipar+ipar-1))
    call calc_seed_evap(ipar, y(1:4))
    call form_RHS(ipar, y(1:4), y(5:5+ipar-1), y(5+ipar:5+ipar+ipar-1), f)

  end subroutine RHS_update

  subroutine jac_dummy(N,X,Y,DFY,LDFY,RPAR,IPAR)
    integer :: N,LDFY,IPAR
    double precision :: X,Y(N),DFY(LDFY,N),RPAR
  end subroutine jac_dummy

  subroutine mas_dummy(N,AM,LMAS,RPAR,IPAR)
    integer :: N, LMAS, IPAR
    double precision :: AM(LMAS,N), RPAR
  end subroutine mas_dummy

  subroutine dfx_dummy(N,X,Y,FX,RPAR,IPAR)
    integer :: N,IPAR
    double precision :: X,Y(N),FX(N),RPAR
  end subroutine dfx_dummy

  subroutine solout(NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
    integer :: NR, LRC, N, IPAR, IRTRN
    double precision :: XOLD, X, Y(N), CONT(LRC), RPAR
  end subroutine solout

end module mini_dihrt_i_ros4
