module mini_cloud_opac
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  type o_table

    character(len=11) :: name
    character(len=50) :: fname

    integer :: n_wl, n_a
    real(dp) :: rhog
    real(dp), allocatable, dimension(:) :: wl, a
    real(dp), allocatable, dimension(:,:) :: ksca, kabs, gg

  end type o_table

  integer, parameter :: b_nwl = 400, b_na = 21

  character(len=50) :: p_2_budaj

  type(o_table), allocatable, dimension(:)  :: o_tab 
  !$omp threadprivate (o_tab)

  logical :: first_call = .True.
  !$omp threadprivate (first_call)

  private :: read_budaj_tables, read_budaj_data, &
    & first_call, locate, bezier_interp, interp_budaj_tables
  public :: opac_budaj_tables

contains

  subroutine opac_budaj_tables(n_dust, sp, T_in, mu_in, P_in, k, k3, n_wl, wl, k_ext, alb, gg)
    implicit none

    integer, intent(in) :: n_dust, n_wl
    character(len=11), dimension(n_dust), intent(in) :: sp
    real(dp), intent(in) :: T_in, mu_in, P_in
    real(dp), dimension(4), intent(in) :: k
    real(dp), dimension(n_dust), intent(in) :: k3
    real(dp), dimension(n_wl), intent(in) :: wl

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    integer :: l, n
    real(dp) :: amean, Vmean, rho, T, mu, P_cgs
    real(dp), dimension(n_dust) :: ksca_n, kabs_n, g_n, b_mix

    if (first_call .eqv. .True.) then
      p_2_budaj = 'budaj_tables/'
      call read_budaj_tables(n_dust,sp)
      first_call = .False.
    end if

    !! Start interpolation routine

    T = T_in
    mu = mu_in
    P_cgs = P_in * 10.0_dp
    rho = (P_cgs * mu * amu)/(kb * T)

    amean = log10(k(2)/k(1) * 1e4_dp)
    Vmean = fourpi3 * k(4)/k(1)
    b_mix(:) = max(k3(:)/k(4),1e-30_dp)
    b_mix(:) = min(k3(:)/k(4),1.0_dp)
    
    do l = 1, n_wl

      k_ext(l) = 0.0_dp
      alb(l) = 0.0_dp
      gg(l) = 0.0_dp
      do n = 1, n_dust
        call interp_budaj_tables(n,amean,wl(l),ksca_n(n),kabs_n(n),g_n(n))
        ksca_n(n) = ksca_n(n) * b_mix(n)
        kabs_n(n) = kabs_n(n) * b_mix(n)
        g_n(n) = g_n(n) * b_mix(n)
        k_ext(l) = k_ext(l) + (ksca_n(n) + kabs_n(n))*o_tab(n)%rhog
        alb(l) = alb(l) + ksca_n(n)*o_tab(n)%rhog
        gg(l) = gg(l) + g_n(n)*ksca_n(n)*o_tab(n)%rhog
      end do

      gg(l) = gg(l)/alb(l)

      alb(l) = alb(l)/k_ext(l)

      k_ext(l) = (k_ext(l) * Vmean * k(1))/rho

      ! Set limits on a and g
      alb(l) = max(0.0_dp, alb(l))
      alb(l) = min(0.99_dp, alb(l))

      gg(l) = max(-0.99_dp, gg(l))
      gg(l) = min(0.99_dp, gg(l))

    end do

  end subroutine opac_budaj_tables

  subroutine interp_budaj_tables(n,a,wl,ksca,kabs,g)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: a, wl

    real(dp), intent(out) :: ksca, kabs, g

    integer :: iwl1, iwl2, iwl3
    integer :: ia1, ia2, ia3
    real(dp), dimension(3) :: lwl, la, lksca, lkabs, lg
    real(dp), dimension(3) :: lkscaa, lkabsa, lga

    ! Find wavelength array triplet and check bounds
    call locate(o_tab(n)%wl(:), o_tab(n)%n_wl, wl, iwl2)

    iwl1 = iwl2 - 1
    iwl3 = iwl2 + 1

    if (iwl1 <= 0) then
      iwl1 = 1
      iwl2 = 2
      iwl3 = 3
    else if (iwl3 > o_tab(n)%n_wl) then
      iwl1 = o_tab(n)%n_wl - 2
      iwl2 = o_tab(n)%n_wl - 1
      iwl3 = o_tab(n)%n_wl
    end if

    lwl(1) = o_tab(n)%wl(iwl1)
    lwl(2) = o_tab(n)%wl(iwl2)
    lwl(3) = o_tab(n)%wl(iwl3)

    ! Check if mean size is within range
    if (a < o_tab(n)%a(1)) then
      ! Interpolate at smallest grain size

      lksca(1) = o_tab(n)%ksca(iwl1,1)
      lksca(2) = o_tab(n)%ksca(iwl2,1)
      lksca(3) = o_tab(n)%ksca(iwl3,1)
      call Bezier_interp(lwl(:), lksca(:), 3, wl, ksca)

      lkabs(1) = o_tab(n)%kabs(iwl1,1)
      lkabs(2) = o_tab(n)%kabs(iwl2,1)
      lkabs(3) = o_tab(n)%kabs(iwl3,1)
      call Bezier_interp(lwl(:), lkabs(:), 3, wl, kabs)

      lg(1) = o_tab(n)%gg(iwl1,1)
      lg(2) = o_tab(n)%gg(iwl2,1)
      lg(3) = o_tab(n)%gg(iwl3,1)
      call Bezier_interp(lwl(:), lg(:), 3, wl, g)

    else if (a > o_tab(n)%a(o_tab(n)%n_a)) then
      ! Interpolate at largest grain size

      lksca(1) = o_tab(n)%ksca(iwl1,o_tab(n)%n_a)
      lksca(2) = o_tab(n)%ksca(iwl2,o_tab(n)%n_a)
      lksca(3) = o_tab(n)%ksca(iwl3,o_tab(n)%n_a)
      call Bezier_interp(lwl(:), lksca(:), 3, wl, ksca)

      lkabs(1) = o_tab(n)%kabs(iwl1,o_tab(n)%n_a)
      lkabs(2) = o_tab(n)%kabs(iwl2,o_tab(n)%n_a)
      lkabs(3) = o_tab(n)%kabs(iwl3,o_tab(n)%n_a)
      call Bezier_interp(lwl(:), lkabs(:), 3, wl, kabs)

      lg(1) = o_tab(n)%gg(iwl1,o_tab(n)%n_a)
      lg(2) = o_tab(n)%gg(iwl2,o_tab(n)%n_a)
      lg(3) = o_tab(n)%gg(iwl3,o_tab(n)%n_a)
      call Bezier_interp(lwl(:), lg(:), 3, wl, g)

    else
      ! Size is within table grid
      ! Perform 2D Bezier interpolation by performing interpolation 4 times

      ! Find size index triplet
      call locate(o_tab(n)%a(:), o_tab(n)%n_a, a, ia2)
      ia1 = ia2 - 1
      ia3 = ia2 + 1

      la(1) = o_tab(n)%a(ia1)
      la(2) = o_tab(n)%a(ia2)
      la(3) = o_tab(n)%a(ia3)

      lksca(1) = o_tab(n)%ksca(iwl1,ia1)
      lksca(2) = o_tab(n)%ksca(iwl2,ia1)
      lksca(3) = o_tab(n)%ksca(iwl3,ia1)
      call Bezier_interp(lwl(:), lksca(:), 3, wl, lkscaa(1)) ! Result at wl, ia1
      lksca(1) = o_tab(n)%ksca(iwl1,ia2)
      lksca(2) = o_tab(n)%ksca(iwl2,ia2)
      lksca(3) = o_tab(n)%ksca(iwl3,ia2)
      call Bezier_interp(lwl(:), lksca(:), 3, wl, lkscaa(2)) ! Result at wl, ia2
      lksca(1) = o_tab(n)%ksca(iwl1,ia3)
      lksca(2) = o_tab(n)%ksca(iwl2,ia3)
      lksca(3) = o_tab(n)%ksca(iwl3,ia3)
      call Bezier_interp(lwl(:), lksca(:), 3, wl, lkscaa(3)) ! Result at wl, ia3
      call Bezier_interp(la(:), lkscaa(:), 3, a, ksca) ! Result at wl, a

      lkabs(1) = o_tab(n)%kabs(iwl1,ia1)
      lkabs(2) = o_tab(n)%kabs(iwl2,ia1)
      lkabs(3) = o_tab(n)%kabs(iwl3,ia1)
      call Bezier_interp(lwl(:), lkabs(:), 3, wl, lkabsa(1)) ! Result at wl, ia1
      lkabs(1) = o_tab(n)%kabs(iwl1,ia2)
      lkabs(2) = o_tab(n)%kabs(iwl2,ia2)
      lkabs(3) = o_tab(n)%kabs(iwl3,ia2)
      call Bezier_interp(lwl(:), lkabs(:), 3, wl, lkabsa(2)) ! Result at wl, ia2
      lkabs(1) = o_tab(n)%kabs(iwl1,ia3)
      lkabs(2) = o_tab(n)%kabs(iwl2,ia3)
      lkabs(3) = o_tab(n)%kabs(iwl3,ia3)
      call Bezier_interp(lwl(:), lkabs(:), 3, wl, lkabsa(3)) ! Result at wl, ia3
      call Bezier_interp(la(:), lkabsa(:), 3, a, kabs) ! Result at wl, a

      lg(1) = o_tab(n)%gg(iwl1,ia1)
      lg(2) = o_tab(n)%gg(iwl2,ia1)
      lg(3) = o_tab(n)%gg(iwl3,ia1)
      call Bezier_interp(lwl(:), lg(:), 3, wl, lga(1)) ! Result at wl, ia1
      lg(1) = o_tab(n)%gg(iwl1,ia2)
      lg(2) = o_tab(n)%gg(iwl2,ia2)
      lg(3) = o_tab(n)%gg(iwl3,ia2)
      call Bezier_interp(lwl(:), lg(:), 3, wl, lga(2)) ! Result at wl, ia2
      lg(1) = o_tab(n)%gg(iwl1,ia3)
      lg(2) = o_tab(n)%gg(iwl2,ia3)
      lg(3) = o_tab(n)%gg(iwl3,ia3)
      call Bezier_interp(lwl(:), lg(:), 3, wl, lga(3)) ! Result at wl, ia3
      call Bezier_interp(la(:), lga(:), 3, a, g) ! Result at wl, a

    end if

  end subroutine interp_budaj_tables

  subroutine read_budaj_tables(n_dust, sp)
    implicit none

    integer, intent(in) :: n_dust
    character(len=11), dimension(n_dust), intent(in) :: sp

    integer :: n

    allocate(o_tab(n_dust))

    o_tab(:)%n_wl = b_nwl 
    o_tab(:)%n_a = b_na

    do n = 1, n_dust

      allocate(o_tab(n)%wl(b_nwl))
      allocate(o_tab(n)%a(b_na))

      allocate(o_tab(n)%ksca(b_nwl,b_na))
      allocate(o_tab(n)%kabs(b_nwl,b_na))
      allocate(o_tab(n)%gg(b_nwl,b_na))

      select case(trim(sp(n)))
      case('CaTiO3')
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'perovskite_opac_all.txt'
        o_tab(n)%rhog = 4.1_dp        
      case('TiO2')
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'perovskite_opac_all.txt'
        o_tab(n)%rhog = 4.1_dp    
      case('Al2O3')
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'corundum_opac_all.txt'
        o_tab(n)%rhog = 2.9_dp
      case('Fe')
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'iron_opac_all.txt'
        o_tab(n)%rhog = 7.874_dp
      case('Mg2SiO4') 
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'forsterite_opac_all.txt'
        o_tab(n)%rhog = 3.0_dp
      case('MgSiO3') 
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'enstatite_opac_all.txt'
        o_tab(n)%rhog = 3.0_dp  
      case('C') 
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'carbon1000_opac_all.txt'
        o_tab(n)%rhog = 1.988_dp
      case('H2O') 
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'waterice_opac_all.txt'
        o_tab(n)%rhog = 1.0_dp
       case('NH3') 
        o_tab(n)%name = sp(n)
        o_tab(n)%fname = 'ammonia_opac_all.txt'
        o_tab(n)%rhog = 0.88_dp
      end select

      call read_budaj_data(n)

    end do

  end subroutine read_budaj_tables

  subroutine read_budaj_data(n)
    implicit none

    integer, intent(in) :: n

    integer :: u, i, j

    real(dp) :: dum1, dum2, k_sca, k_abs, gg

    open(newunit=u,file=trim(p_2_budaj)//trim(o_tab(n)%fname),action='read')

    print*, 'reading budaj table @: ', trim(p_2_budaj)//trim(o_tab(n)%fname)

    do i = 1, o_tab(n)%n_a
      do j = 1, o_tab(n)%n_wl
        read(u,*) dum1, dum2, k_sca , k_abs, gg
        o_tab(n)%ksca(j,i) = max(k_sca, 1e-30_dp)
        o_tab(n)%kabs(j,i) = max(k_abs, 1e-30_dp)
        o_tab(n)%gg(j,i) = gg
        !print*, dum1, dum2, k_sca, k_abs, gg
        if (i == 1) then
          o_tab(n)%wl(j) = dum2
        end if
      end do
      o_tab(n)%a(i) = dum1 
    end do

    close(u)

  end subroutine read_budaj_data


  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

  subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(dp), dimension(n), intent(in) :: arr
    real(dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

end module mini_cloud_opac