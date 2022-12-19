subroutine fooF() bind(C, name='fooFortran')
  use iso_c_binding
  implicit none
  write(*,*) 'Fortran fooF'
end
