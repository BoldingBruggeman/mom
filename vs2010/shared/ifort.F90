subroutine flush_ifort(unit)
    implicit none
    integer,intent(in) :: unit
    flush unit
end subroutine flush_ifort

subroutine getenv_ifort(ename,evalue)
    use ifport
    implicit none
    character*(*) :: ename,evalue
#undef GETENV
    call GETENV(ename,evalue)
end subroutine getenv_ifort