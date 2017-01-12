module newton_minimizing_mod
  use healpix_types
  use quiet_utils
  implicit none


contains

  SUBROUTINE newton_min(p,gtol,step,iter,dfunc,ddfunc,ierr, rangelimit)
  ! Searches for minimum of function with Newton's method.
  ! Must supply the derivative (dfunc) and the second derivative (ddfunc)
  ! p -parameter values
  ! gtol -threshold for stopping criterium
  ! step -step length factor. If in doubt, try 1.d0
  ! iter -number of iterations needed to reach stopping criterium (max=1000)
  ! ierr not in use yet, but can be used as a error flag.
    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: iter
    REAL(dp), INTENT(IN) :: gtol, step
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
    integer(i4b), optional, intent(out)   :: ierr
    real(dp), optional, intent(in)        :: rangelimit
    integer(i4b)                          :: its, itmax=1000
    real(dp), dimension(size(p))          :: delta_p
    real(dp)                              :: healnan=-1.6375d30
    INTERFACE
       FUNCTION dfunc(p)
         USE healpix_types
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: p
         REAL(dp), DIMENSION(size(p)) :: dfunc
       END FUNCTION dfunc
       FUNCTION ddfunc(p)
         USE healpix_types
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), INTENT(IN) :: p
         REAL(dp), DIMENSION(size(p),size(p)) :: ddfunc
       END FUNCTION ddfunc
    END INTERFACE

    do its=1,itmax
       iter=its
       delta_p=matmul(ddfunc(p),dfunc(p))
       p=p - step*delta_p
       if (abs(delta_p(1)/p(1)) < gtol) return
       if (present(rangelimit)) then
          if( (p(1).LT.-rangelimit).or.(p(1).GT.rangelimit) ) then
             p(1)= 101.d0
             return
          end if
       end if
       write(*,*) its, p(1)
    end do

  end SUBROUTINE newton_min


end module newton_minimizing_mod
