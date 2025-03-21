!
!   Module for common parameters
!
MODULE mCommon
     
    IMPLICIT NONE
!
!    Parameters
!
! ----------------------------------------------------------------------
!
    integer, parameter :: rk = selected_real_kind(8)
    real(rk), parameter :: pi = 3.141592653589793238462643383_rk
    real(rk), parameter :: vSmall = 1.0E-10_rk
    logical, parameter :: debug = .FALSE.

END MODULE mCommon