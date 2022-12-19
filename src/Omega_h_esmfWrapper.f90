subroutine fooF() bind(C, name='fooFortran')
  use iso_c_binding
  use ESMF
  implicit none
  type(ESMF_VM) :: vm
  integer :: rc, result

  write(*,*) 'Fortran fooF'
  call  ESMF_Initialize(vm=vm, defaultlogfilename="fooF.Log", &
                    logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
end
