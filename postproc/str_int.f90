module str_int
  !
  ! convert a string to an integer for base 10 numbers
  !
  ! adapted from decode subroutine by
  ! Dan Nagle @ Purple Sage Computing Solutions, Inc.
  ! <http://users.erols.com/dnagle/>
  !
  ! Original code released under GPL.
  ! I've modified (removed) a lot of it.
  ! Don't know license status of this code.
  ! 
  implicit none
  private
  public decode

contains

  pure subroutine decode(i, str)

    integer, intent(out) :: i               ! output integer
    character (len = *), intent(in) :: str  ! input string

    !  local
    !  
    character (len = len(str)) :: str_buff
    integer :: jstr

    ! character set, base 10
    ! 
    integer :: base
    character (len = 10) :: ttable           

    base = 10
    ttable = '0123456789'

    !  check input
    !  
    str_buff = adjustl(str)        ! process left to right
    i = 0                          ! string not yet read

    !  scan str
    !  
    each_char: do while(str_buff /= ' ')   ! scan thru string
       jstr = index(ttable, str_buff(1:1)) ! look up character
       i = i*base + (jstr - 1)             ! add a digit
       str_buff = str_buff(2:)             ! next character
    enddo each_char

  end subroutine decode

end module str_int


