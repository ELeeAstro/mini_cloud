module mini_cloud_i_limex
  use mini_cloud_precision
  use mini_cloud_class
  use mini_cloud_saturation
  use mini_cloud_nucleation
  use mini_cloud_chem
  use mini_cloud_rhs
  implicit none

  logical :: first_call = .True.
  !$omp threadprivate (first_call)

  integer :: ipar

  public ::  mini_cloud_limex, RHS_update, jac_dummy

contains

  subroutine mini_cloud_limex(n_dust, T_in, P_in, t_end, sp, k, k3, VMR)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T_in, P_in, t_end
    character(len=11), dimension(n_dust), intent(in) :: sp

    real(dp), dimension(4), intent(inout) :: k
    real(dp), dimension(n_dust), intent(inout) :: k3
    real(dp), dimension(n_dust), intent(inout) :: VMR

    integer :: ncall

    ! Time controls
    real(dp) :: t_begin, t_now, t_min, t_step, t_step0

    ! LIMEX variables
    integer :: neq
    integer, dimension(3) :: ifail
    integer, dimension(30) :: iopt
    integer, allocatable, dimension(:) :: ipos
    real(dp), dimension(5) :: ropt
    real(dp), allocatable, dimension(:) :: y, ys
    real(dp), allocatable, dimension(:)  :: rtol, atol

    real(dp) :: a_check, Ar_check, V_check

    if (first_call .eqv. .True.) then
      call mini_cloud_init(n_dust, sp)
      first_call = .False.
    end if

    neq = 4 + n_dust + n_dust
    ipar = n_dust

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
       & .and. abs(sum(d_sp(:)%chis)) < 1.0e-30_dp)) then
         return
    end if

    ! -----------------------------------------
    ! ***  parameters for the limex-solver  ***
    ! -----------------------------------------

    t_step = 1.0e-12_dp    ! First timestep try
    t_step0 = t_step       ! Keep track of inital timestep
    t_min = 1.0e-16_dp      ! Minimum timestep

    !! Allocate limex arrays
    allocate(ipos(neq),y(neq),ys(neq),rtol(neq),atol(neq))

    ! LIMEX opt parameters
    iopt(:) = 0
    iopt(1)  = 0       ! how much output? (0=no, 1=standard, 2=more)
    iopt(2)  = 0       ! unit number for output (0:default=6)
    iopt(3)  = 0       ! solution output? (0=no)
    iopt(4)  = 0       ! unit number for solution output (0:default=6)
    iopt(5)  = 0       ! nonsigular matrix BB (0=singular, 1=nonsingular)
    iopt(6)  = 0       ! determination of initial values for FF,BB (1=yes)
    iopt(7)  = 0       ! analytic Jacobian? (0=numerical, 1=analytic)
    iopt(8)  = neq     ! Lower bandwidth of the Jacobian (Ndim=full)
    iopt(9)  = neq     ! Upper bandwidth of the Jacobian (Ndim=full)
    iopt(10) = 1       ! reuse of the Jacobian? (1=yes)
    iopt(11) = 1       ! Switch for error toleranz (0=scalars, 1=vectors)
    iopt(12) = 0       ! Return after one integration time step? (0=no)
    iopt(13) = 0 !2      ! Dense output option (0=no dense output)
    iopt(14) = 0 !2     ! The number of equidistant dense output points
    iopt(15) = 0       ! unit for dense output
    iopt(16) = 0       ! type of call (0=first call, 1=continuation call)
    iopt(17) = 0       ! bevahior at t1 (0=exactly up to t1)
    iopt(18) = 0       ! Generation of of a PostScript plot of the Jacobian

    ropt(1)  = t_end !(t_end - t_begin)/4.0_dp ! maximum allowed stepsize
    ropt(2)  = 0.0     ! maximal distance of dense outputs (if iopt(13)=3)
    ropt(3)  = 0.0     ! upper limit for t (if iopt(17)=1)

    ipos(:) = 1
    rtol(:) = 1e-3_dp
    atol(:) = 1e-30_dp
    ifail(:) = 0


    ! Prepare all species for the solver
    y(1:4) = max(k(:),1e-30_dp)
    y(5:5+n_dust-1) = max(k3(:),1e-30_dp)
    y(5+n_dust:5+n_dust+n_dust-1) = max(VMR(:),1e-30_dp)

    t_begin = 0.0_dp
    t_now = t_begin

    ncall = 1

    do while((t_now < t_end) .and. ncall <= 100)

      call calc_saturation(ipar, y(5+ipar:5+ipar+ipar-1))
      call calc_nucleation(ipar, y(5+ipar:5+ipar+ipar-1))
      call calc_chem(ipar, y(1:4), y(5:5+ipar-1), y(5+ipar:5+ipar+ipar-1))
      call calc_seed_evap(ipar, y(1:4))
      call form_RHS(ipar, y(1:4), ys)

      y(1:5+n_dust-1) = y(1:5+n_dust-1)/nd_atm
      ys(1:5+ipar-1) = ys(1:5+ipar-1)/nd_atm


      call LIMEX (neq, RHS_update, jac_dummy, t_now, t_end, y, ys, &
        & rTol, aTol, t_step, Iopt, Ropt, IPos, IFail )

      y(1:5+n_dust-1) = y(1:5+n_dust-1) * nd_atm

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

      ! Start solving problems with LIMEX integration
      select case (ifail(1))
        case (0)
          ! Do nothing and continue
        case (-46)
          ! Timestep problem (too large or too many reductions)
          t_step = t_step /10.0_dp ! Lower timestep by a factor
          iopt(10) = 0                         ! no reuse of the Jacobian
          iopt(16) = 0                         ! new initial call
          ifail(:) = 0
          if (t_step < t_min) then
            t_step = t_step0                    ! Do not go below minimum timestep
          end if

        case (-48)
          ! Timestep problem (too many iterations required)
          t_step = t_step                   ! Reset timestep to max
          iopt(10) = 0                        ! no reuse of the Jacobian
          iopt(16) = 0                        ! new initial call
          ifail(:) = 0

        case (-50)
          ! Timestep too small
          t_step = t_step0                     ! Reset timestep to initial
          iopt(10) = 0                      ! no reuse of the Jacobian
          iopt(16) = 0                         ! new initial call
          ifail(:) = 0

        case default
          ! If the error is unknown, print and stop
          print*, ifail(1), t_now, t_step, "!!! Failure unknown !!!"
          stop
        end select

    end do

    k(:) = max(y(1:4),1e-30_dp)
    k3(:) = max(y(5:5+n_dust-1),1e-30_dp)
    VMR(:) = max(y(5+n_dust:5+n_dust+n_dust-1),1e-30_dp)

    deallocate(ipos, y, ys, rtol, atol)

  end subroutine mini_cloud_limex


  subroutine RHS_update(n_in, nz, t_in, y, f, b, ir, ic, iflag)
    !! Subroutine used by LIMEX solver for intermediate timestep RHS values
    implicit none

    integer, intent(in) ::  n_in
    integer, intent(out) :: nz
    integer, intent(inout) :: iflag
    integer, dimension(n_in), intent(out) :: ir, ic
    integer :: n, d

    real(dp), intent(in) :: t_in
    real(dp), dimension(n_in),intent(out) :: b
    real(dp), dimension(n_in),intent(inout) :: y
    real(dp), dimension(n_in),intent(out) :: f

    ! Keep dimensions the same
    nz = n_in

    do n = 1, n_in
      ir(n) = n   ! non-zero entries of BB in coordinate format
      ic(n) = n
      b(n) = 1.0_dp   ! BB is and remains unity
    end do

    y(1:5+ipar-1) = y(1:5+ipar-1) * nd_atm

    y(5:5+ipar-1) = max(y(5:5+ipar-1),1.0e-30_dp)

    call calc_saturation(ipar, y(5+ipar:5+ipar+ipar-1))
    call calc_nucleation(ipar, y(5+ipar:5+ipar+ipar-1))
    call calc_chem(ipar, y(1:4), y(5:5+ipar-1), y(5+ipar:5+ipar+ipar-1))
    call calc_seed_evap(ipar, y(1:4))
    call form_RHS(ipar, y(1:4), f)

    y(1:5+ipar-1) = y(1:5+ipar-1)/nd_atm
    f(1:5+ipar-1) = f(1:5+ipar-1)/nd_atm

    ! Update iflag
    iflag = 0

  end subroutine RHS_update

  pure subroutine jac_dummy( n_in, t_in, y, ys, Jac, LDJac, ml, mu, Full_or_Band, JacInfo )
    !! A dummy Jacobian subroutine for LIMEX
    implicit none

    integer, intent(inout) :: n_in, LDJac, ml, mu, Full_or_Band, JacInfo
    double precision, intent(inout) :: t_in
    double precision, dimension(LDJac,*), intent(inout) :: Jac
    double precision, dimension(*), intent(inout) :: y , ys !(*) array defintions??!!??
    !  print*, 'In Jac'

  end subroutine jac_dummy

end module mini_cloud_i_limex
