module cli_progress
  use, intrinsic :: iso_fortran_env, only: int32, real64, output_unit, error_unit
  implicit none
  private
  public :: progress_begin, progress_update, progress_end, spinner_tick

  integer(int32) :: t0 = 0_int32
  integer(int32) :: rate = 1_int32

contains

  subroutine progress_begin()
    integer(int32) :: c
    call system_clock(count_rate=rate)
    call system_clock(c)
    t0 = c
  end subroutine progress_begin

  subroutine progress_update(i, n, sim_time, max_ratio, width, to_stderr)
    integer(int32), intent(in) :: i, n
    integer,        intent(in), optional :: width
    logical,        intent(in), optional :: to_stderr
    real(real64),   intent(in) :: sim_time      ! current time variable
    real(real64),   intent(in) :: max_ratio     ! maxval(del(:)/t_step)
    integer :: W, hashes, dots, pct
    integer(int32) :: c, elapsed, remaining
    integer        :: u

    u = merge(error_unit, output_unit, present(to_stderr) .and. to_stderr)
    W = merge(width, 40, present(width))

    pct = int(100.0d0 * dble(i) / max(1_int32, n))
    hashes = int(dble(W) * dble(i) / max(1_int32, n))
    dots   = max(0, W - hashes)

    call system_clock(c)
    elapsed = c - t0
    if (i > 0_int32 .and. rate > 0_int32) then
      remaining = int( (dble(n - i) / max(1.0d0, dble(i))) * dble(elapsed) )
    else
      remaining = -1_int32
    end if

    write(u,'(a)',advance='no') achar(13)
    write(u,'("[",a,"] ",i3,"%",1x,"t=",f8.2," max=",es10.3,1x,"(eta ",a,")")',advance='no') &
        repeat("#", hashes)//repeat(".", dots), pct, sim_time, max_ratio, trim(eta_string(remaining))
    call flush(u)
  end subroutine progress_update

  subroutine progress_end(to_stderr)
    logical, intent(in), optional :: to_stderr
    integer :: u
    u = merge(error_unit, output_unit, present(to_stderr) .and. to_stderr)
    write(u,*)                              ! newline to finish the bar
    call flush(u)
  end subroutine progress_end

  subroutine spinner_tick(step, to_stderr)
    integer, intent(in) :: step
    logical, intent(in), optional :: to_stderr
    character(len=*), parameter :: frames(4) = (/'|','/','-','\'/)
    integer :: u
    u = merge(error_unit, output_unit, present(to_stderr) .and. to_stderr)
    write(u,'(a)',advance='no') achar(13)//frames(1+mod(step,4))//' working...'
    call flush(u)
  end subroutine spinner_tick

  pure function eta_string(ticks) result(s)
    integer(int32), intent(in) :: ticks  ! in system_clock ticks
    character(len=16) :: s
    real(real64) :: secs
    if (ticks < 0_int32 .or. rate <= 0_int32) then
      s = ' --:--'
    else
      secs = dble(ticks) / dble(rate)
      write(s,'(i2.2, ":", i2.2)') int(secs)/60, mod(int(secs),60)
    end if
  end function eta_string

end module cli_progress
