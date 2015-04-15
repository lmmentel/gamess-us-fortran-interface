module VarModule
  use ISO_FORTRAN_ENV
  implicit none
!	----------------------
!	Define KIND parameters
! ----------------------
  integer,     parameter :: IK = kind(1)
  integer(IK), parameter :: LK2 = kind(2)
  integer(IK), parameter :: LK8 = kind(8)
  integer(IK), parameter :: SP = kind(1.0)
  integer(IK), parameter :: DP = kind(1.0D0)
! ----------------------------
! Specify kinds used in Gamess
! ----------------------------
  integer(IK), parameter :: I8 = 8
  integer(IK), parameter :: IPgamess = 8  !4 = 32 bit integer, 8 = 64 bit integer
  integer(IK), parameter :: DPgamess = DP
  integer(IK), parameter :: LKgamess = 8
! -------------------------------
! Define length of default string
! -------------------------------
  integer(IK), parameter :: LCHARS = 60
!	---------------------------------------------
!	Set std unit numbers to system specified ones
!	---------------------------------------------
  integer(IK), parameter :: stdin  = INPUT_UNIT
  integer(IK), parameter :: stdout = OUTPUT_UNIT 
  integer(IK), parameter :: stderr = ERROR_UNIT
!	---------
!	Constants
!	---------
  real(DP), parameter :: PI    = 3.1415926535897932384626433_DP
  real(DP), parameter :: TWOPI = PI+PI
  real(DP), parameter :: GOLD  = 1.6180339887498948482045868_DP
  real(DP), parameter :: CGOLD = 0.3819660112501051517954132_DP
! ------------------------------  
! some common format definitions
! ------------------------------  
 character(len=20), parameter :: plain   = '(a100)' 
 character(len=20), parameter :: mofmt   = '(i6,es22.12)' 
 character(len=20), parameter :: orbefmt = '(a4,f10.6)' 
 character(len=20), parameter :: gmsvfmt = '(i2,i3,5(es15.8))' 
! ------------------
! Define vector norm
! ------------------
  interface Norm2
    module procedure Norm2Real
    module procedure Norm2Complex
  end interface

  private Norm2Real, Norm2Complex

  interface matTrans
      module procedure subroutineMatrixTransform
  end interface

  private :: subroutineMatrixTransform

contains

  real(DP) function Norm2Real(vec)
    real(DP), intent(in) :: vec(:)
    Norm2Real = dot_product(vec,vec)
  end function
  real(DP) function Norm2Complex(vec)
    complex(DP), intent(in) :: vec(:)
    Norm2Complex = real(dot_product(vec,vec), kind=DP) !dot_product already takes conjugate
  end function

  function OuterProduct(vec1, vec2) result(matrix)
    real(DP), intent(in) :: vec1(:), vec2(:)
    real(DP)             :: matrix(size(vec1),size(vec2))
    matrix = spread(vec1, 2, size(vec2)) * spread(vec2, 1, size(vec1))
  end function

  elemental complex(DP) function ImExp(x)
    real(DP), intent(in) :: x

    ImExp = exp(cmplx(0.0_DP, x, kind=DP))
  end function

  elemental real(DP) function phase(re, im)
    real(DP), intent(in) :: re, im

    if (re > epsilon(1.0_DP)) then
      phase = atan(im/re)
    elseif (re < -epsilon(1.0_DP)) then
      if (im < 0.0_DP) then
        phase = atan(im/re) - PI
      else
        phase = atan(im/re) + PI
      endif
    else
      if (im > 0.0_DP) then
        phase =  0.5_DP*PI
      elseif (im < 0.0_DP) then
        phase = -0.5_DP*PI
      else
        phase =  0.0_DP
      endif
    endif
  end function

   function functionMatrixTransform(X,A,Y) result(R)
! R = X^T * A * Y
    real(dp), dimension(:,:), intent(in)  :: X,Y,A
    real(dp)                              :: R(size(A,1),size(A,2)) 
    real(dp), dimension(:,:), allocatable :: tmp 

  allocate(tmp(size(A,1),size(A,2)))

  tmp = matmul(A,Y)
  R   = matmul(transpose(X),tmp)
 
  deallocate(tmp)
  end function FunctionMatrixTransform
  
  subroutine SubroutineMatrixTransform(X,A,Y) 
! A = X^T * A * Y
    real(dp), dimension(:,:), intent(in)    :: X,Y 
    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:,:), allocatable   :: tmp 

  allocate(tmp(size(A,1),size(A,2)))

  tmp = matmul(A,Y)
  A   = matmul(transpose(X),tmp)

  deallocate(tmp)
  end subroutine SubroutineMatrixTransform

end module varModule
