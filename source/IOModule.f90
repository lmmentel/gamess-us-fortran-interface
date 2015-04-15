module ioModule
 use varModule 

 private

 interface matprint 
    module procedure print_eigenproblem
    module procedure print_vector
    module procedure print_vector_with_header
    module procedure print_matrix
    module procedure print_matrix_with_header
 end interface 

 interface print_gamess_vecs
    module procedure print_gamess_vec
 end interface

 interface printHeader
     module procedure printHeader
 end interface

 public :: matprint, print_gamess_vecs, printHeader

 contains

 subroutine printHeader(header)
  character(*), intent(in) :: header
  character(len=30)        :: headerfmt

  lspace = 7 + int((50-len(header))/2)
  write(headerfmt,'("(",i2,"x,a)")') lspace
  write(*,*)
  write(*,headerfmt) repeat('=',len(trim(header)))
  write(*,headerfmt) trim(header)
  write(*,headerfmt) repeat('=',len(trim(header)))
 end subroutine printHeader

 subroutine print_eigenproblem(M, V)
  integer(ik), parameter :: vec_per_page = 5 
  integer(ik)          :: i, j, iel, k, l
  real(dp), intent(in) :: M(:,:)
  real(dp), intent(in) :: V(:)

 if (mod(size(M,2),vec_per_page) == 0 ) then 
   iel = size(M,2)/vec_per_page
 else
   iel = size(M,2)/vec_per_page+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
     write(*,'(/7x,5(f10.5)/)') (V(j+(i-1)*10),j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
     write(*,'(/7x,5(f10.5)/)') (V(j+(i-1)*5),j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_eigenproblem

 subroutine print_vector(V)
  real(dp), intent(in) :: V(:)
  integer(ik)          :: i

  do i = 1, size(V)
    write(*,'(1x,i3,f10.5)') i, V(i)
  enddo
 end subroutine print_vector

 subroutine print_vector_with_header(V, header)
  real(dp),     intent(in) :: V(:)
  character(*), intent(in) :: header
  integer(ik)              :: i

  write(*,'(/4x,a)') repeat('=',len(trim(header))) 
  write(*,'(4x,a)') trim(header) 
  write(*,'(4x,a/)') repeat('=',len(trim(header))) 

  do i = 1, size(V)
    write(*,'(1x,i3,f10.5)') i, V(i)
  enddo
 end subroutine print_vector_with_header

 subroutine print_matrix(M)
  real(dp), intent(in) :: M(:,:)
  integer(ik)          :: i,j,iel,k,l

 if (mod(size(M,2),5) == 0 ) then 
   iel = size(M,2)/5
 else
   iel = size(M,2)/5+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_matrix

 subroutine print_matrix_with_header(M, header)
  real(dp),     intent(in) :: M(:,:)
  character(*), intent(in) :: header
  integer(ik)              :: i,j,iel,k,l
  character(len=30)        :: headerfmt

  lspace = 7 + int((50-len(header))/2)
  write(headerfmt,'("(",i2,"x,a)")') lspace  
  write(*,headerfmt) repeat('=',len(trim(header))) 
  write(*,headerfmt) trim(header) 
  write(*,headerfmt) repeat('=',len(trim(header))) 

 if (mod(size(M,2),5) == 0 ) then 
   iel = size(M,2)/5
 else
   iel = size(M,2)/5+1
 endif

 do i = 1, iel 
   if(i < iel) then
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5)
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5)
       enddo
       write(*,*)
   else   
     write(*,'(/10x,5(i5,5x))') (j+(i-1)*5,j=1,5-(iel*5-size(M,2)))
       do k = 1,size(M,1)
         write(*,'(1x,i3,3x,5f10.5)') k, (M(k,l+(i-1)*5),l=1,5-(iel*5-size(M,2)))
       enddo
   endif
 enddo    
 write(*,*)
 end subroutine print_matrix_with_header

 subroutine print_gamess_vec(coeffs)
  real(dp),           intent(in) :: coeffs(:,:)
  !character(len=100), intent(in) :: vec_file
  integer(ik) :: i, j, l, ilab, llab, nlines
  integer(ik), parameter :: gunit=16

  if (mod(size(coeffs,1),5) == 0) then
    nlines = int(size(coeffs,1)/5)
  else
    nlines = int(size(coeffs,1)/5,ik)+1
  endif

  print *, 'nlines = ', nlines

  open(gunit, file='vectors.out', status="replace")
 ! write header and then orbitals
  write(gunit,'(/a,i3,a/)') ' $guess guess=moread norb=',size(coeffs,2),' punmo=.true. prtmo=.true. $end'
  write(gunit,'(a)') ' $vec'
  do i = 1, size(coeffs,2)
    if (i >= 100) then
      ilab = mod(i, 100)
    else
      ilab = i
    endif
    do l = 1, nlines
      if (l > 1000) then
        llab = mod(l, 1000)
      else
        llab = l
      endif
      if (l < nlines) then
        write(gunit,gmsvfmt) ilab,llab,(coeffs(j+(l-1)*5,i),j=1,5)
      else
        write(gunit,gmsvfmt) ilab,llab,(coeffs(j+(l-1)*5,i),j=1,5-(l*5-size(coeffs,1)))
      endif
    enddo
  enddo
  write(gunit,'(a)') ' $end'
  close(gunit)

 end subroutine print_gamess_vec


end module ioModule
