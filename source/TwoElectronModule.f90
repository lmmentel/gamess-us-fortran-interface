module TwoElectronModule
  use VarModule

  implicit none

  private

  interface ReadInAO
      module procedure ReadInAO
  end interface

  interface ReadInMO
      module procedure ReadInMO
  end interface

  interface print4index
      module procedure print4index
  end interface

  interface getTwoElEnergy
      module procedure exactenergy
  end interface


  public :: ReadInAO, ReadInMO, print4index, getTwoElEnergy

contains

   integer(I8) function address(i,j,k,l)
   integer(I8) :: i, j, k, l, ij, kl

    ij = max(i,j)*(max(i,j)-1)/2 + min(i,j)
    kl = max(k,l)*(max(k,l)-1)/2 + min(k,l)

    address = max(ij,kl)*(max(ij,kl)-1)/2 + min(ij,kl)

  end function

 function factor(i,j,k,l) result(fac)
  integer(I8), intent(in)  :: i, j, k, l
  real(DP)                 :: fac
  if(i == j .and. k == l .and. i == k) then
    fac = 1.0_DP
  elseif(i == j .and. k == l) then
    fac = 2.0_DP
  elseif((i == k .and. j == l) .or. (i==j.and.i==k) .or. (j==k.and.j==l).or.(i==j.or.k==l)) then
    fac = 4.0_DP
  else
    fac = 8.0_DP
  endif
 end function

  subroutine ReadInAO(TwoIntAO, filename)
! read two-electron integrals over AO's
  integer(IK), parameter :: twoein = 9

    real(DP),         intent(inout)   :: TwoIntAO(:)
    character(len=*), intent(in)      :: filename

    real(DP),          allocatable  :: buffer(:)
    integer(IPgamess), allocatable  :: indexBuffer(:)

    real(DP)          :: temp

    integer(IK)       :: readStatus
    integer(IK)       :: m
    integer(IPgamess) :: i,j,k,l
    integer(IK)       :: bufLength, twoIntIndexBufSize, twoIntBufferSize
    integer(IPgamess) :: label, label1, label2
    integer(IPgamess) :: length, nintmx, labsiz, twoeao
    logical           :: largeLabels

    nintmx = 15000                               ! gamess parameter controling read buffer
    labsiz = 1_IPgamess                          ! gamess parameter controling the index range, for maxao > 256 set to 2
    twoeao = twoein

    if (labsiz /= 1_IPgamess .and. labsiz /= 2_IPgamess) then
      write(*,*) 'RdTwoIntAO:  CONFUSION IN LABSIZ! '
      stop
    endif

    largeLabels = (labsiz == 2_IPgamess)

    twoIntBufferSize = int(nintmx, kind=IK)

    twoIntIndexBufSize = twoIntBufferSize
    if (largeLabels) then
      if (IPgamess == 4) twoIntIndexBufSize = 2*twoIntIndexBufSize
    else
      if (IPgamess == 8) twoIntIndexBufSize = (twoIntIndexBufSize + 1) / 2
    endif

    open(unit=twoeao, file=trim(filename), status='old', form='unformatted')
    rewind(twoeao)

    allocate(buffer(twoIntBufferSize))
    allocate(indexBuffer(twoIntIndexBufSize))

    length = 1_IPgamess
    do while (length > 0_IPgamess)
      Read(twoeao,iostat=readStatus) length,indexBuffer,buffer
      if (readStatus /= 0) then
        if (readStatus == 2) then
          write(*,*) 'RdTwoIntAO: ENCOUNTERED UNEXPECTED END WHILE READING TWO-ELECTRON FILE'
        else
          write(*,*) 'RdTwoIntAO: ENCOUNTERED I/O ERROR WHILE READING TWO-ELECTRON FILE'
        endif
        stop
      endif

      bufLength = abs(int(length, kind=IK))
      if (bufLength > twoIntBufferSize) stop

      do m=1, bufLength
        if (IPgamess == 4) then                  ! 32-bit Gamess integers
          if (largeLabels) then       
            label1 = indexbuffer(2*m-1) 
            label2 = indexBuffer(2*m)
            i = ishft(label1, -16_IPgamess)
            j = iand( label1,  65535_IPgamess)
            k = ishft(label2, -16_IPgamess)
            l = iand( label2,  65535_IPgamess)
          else
            label = indexBuffer(m)
            i =      ishft(label, -24_IPgamess)
            j = iand(ishft(label, -16_IPgamess), 255_IPgamess)
            k = iand(ishft(label,  -8_IPgamess), 255_IPgamess)
            l = iand(      label,                255_IPgamess)
          endif  
        else                                     ! 64-bit Gamess integers
          if (largeLabels) then
            label = indexBuffer(m)
            i = int(     ishft(label,   -48_IPgamess),                  kind=IK)
            j = int(iand(ishft(label,   -32_IPgamess), 65535_IPgamess), kind=IK)
            k = int(iand(ishft(label,   -16_IPgamess), 65535_IPgamess), kind=IK)
            l = int(iand(      label,                  65535_IPgamess), kind=IK)  
          else
            if (mod(m,2) == 0) then
              label = indexBuffer(m/2)
              i = int(iand(ishft(label, -24_IPgamess ), 255_IPgamess), kind=IK)
              j = int(iand(ishft(label, -16_IPgamess ), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label,  -8_IPgamess ), 255_IPgamess), kind=IK)
              l = int(iand(      label,                 255_IPgamess), kind=IK)
            else
              label = indexBuffer(m/2+1)
              i = int(     ishft(label, -56_IPgamess),                kind=IK)
              j = int(iand(ishft(label, -48_IPgamess), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label, -40_IPgamess), 255_IPgamess), kind=IK)
              l = int(iand(ishft(label, -32_IPgamess), 255_IPgamess), kind=IK)
            endif
          endif
        endif

        temp = buffer(m)
        TwoIntAO(address(i,j,k,l)) = temp

      enddo
    enddo

    deallocate(buffer)
    deallocate(indexBuffer)
    close(twoeao)
  end subroutine

  subroutine ReadInMO(TwoIntMO, filename)
  integer(IK), parameter :: LKgamess = 8
  integer(IK), parameter :: twoein = 9

    real(DP),         intent(inout)   :: TwoIntMO(:)
    character(len=*), intent(in)      :: filename

    real(DP),          allocatable  :: buffer(:)
    integer(IPgamess), allocatable  :: indexBuffer(:)

    real(DP)          :: temp

    integer(IK)       :: readStatus
    integer(IK)       :: m
    integer(IPgamess)       :: i,j,k,l
    integer(IK)       :: bufLength, twoIntIndexBufSize, twoIntBufferSize
    integer(IPgamess) :: label, label1, label2
    integer(IPgamess) :: length, nintmx, labsiz, twoemo
    logical           :: largeLabels

    nintmx = 15000
    labsiz = 1
    twoemo = twoein

    if (labsiz /= 1_IPgamess .and. labsiz /= 2_IPgamess) then
      write(*,*) 'RdTwoIntAO:  CONFUSION IN LABSIZ! '
      stop
    endif

    largeLabels = (labsiz == 2_IPgamess)

    twoIntBufferSize = int(nintmx, kind=IK)

    twoIntIndexBufSize = twoIntBufferSize
    if (largeLabels) then
      if (IPgamess == 4) twoIntIndexBufSize = 2*twoIntIndexBufSize
    else
      if (IPgamess == 8) twoIntIndexBufSize = (twoIntIndexBufSize + 1) / 2
    endif

    open(unit=twoemo, file=trim(filename), status='old', form='unformatted')
    rewind(twoemo)
    read(twoemo)

    allocate(buffer(twoIntBufferSize))
    allocate(indexBuffer(twoIntIndexBufSize))

    length = 1_IPgamess
    do while (length > 0_IPgamess)
      Read(twoemo,iostat=readStatus) length,indexBuffer,buffer
      if (readStatus /= 0) then
        if (readStatus == 2) then
          write(*,*) 'RdTwoIntMO: ENCOUNTERED UNEXPECTED END WHILE READING TWO-ELECTRON FILE'
        else
          write(*,*) 'RdTwoIntMO: ENCOUNTERED I/O ERROR WHILE READING TWO-ELECTRON FILE'
        endif
        stop
      endif

      bufLength = abs(int(length, kind=IK))
      if (bufLength > twoIntBufferSize) stop

      do m=1, bufLength
        if (IPgamess == 4) then                  ! 32-bit Gamess integers
          if (largeLabels) then
            label1 = indexbuffer(2*m-1)
            label2 = indexBuffer(2*m)
            i = ishft(label1, -16_IPgamess)
            j = iand( label1,  65535_IPgamess)
            k = ishft(label2, -16_IPgamess)
            l = iand( label2,  65535_IPgamess)
          else
            label = indexBuffer(m)
            i =      ishft(label, -24_IPgamess)
            j = iand(ishft(label, -16_IPgamess), 255_IPgamess)
            k = iand(ishft(label,  -8_IPgamess), 255_IPgamess)
            l = iand(      label,                255_IPgamess)
          endif
        else                                     ! 64-bit Gamess integers
          if (largeLabels) then
            label = indexBuffer(m)
            i = int(     ishft(label,   -48_IPgamess),                  kind=IK)
            j = int(iand(ishft(label,   -32_IPgamess), 65535_IPgamess), kind=IK)
            k = int(iand(ishft(label,   -16_IPgamess), 65535_IPgamess), kind=IK)
            l = int(iand(      label,                  65535_IPgamess), kind=IK)
          else
            if (mod(m,2) == 0) then
              label = indexBuffer(m/2)
              i = int(iand(ishft(label, -24_IPgamess ), 255_IPgamess), kind=IK)
              j = int(iand(ishft(label, -16_IPgamess ), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label,  -8_IPgamess ), 255_IPgamess), kind=IK)
              l = int(iand(      label,                 255_IPgamess), kind=IK)
            else
              label = indexBuffer(m/2+1)
              i = int(     ishft(label, -56_IPgamess),                kind=IK)
              j = int(iand(ishft(label, -48_IPgamess), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label, -40_IPgamess), 255_IPgamess), kind=IK)
              l = int(iand(ishft(label, -32_IPgamess), 255_IPgamess), kind=IK)
            endif
          endif
        endif

        temp = buffer(m)
        TwoIntMO(address(i,j,k,l)) = temp

      enddo
    enddo

    deallocate(buffer)
    deallocate(indexBuffer)
    close(twoemo)
  end subroutine


 function exactEnergy(twoe, twordm, nb) result(energy)
  real(DP),    intent(in) :: twoe(:), twordm(:)
  integer(I8), intent(in) :: nb
  integer(I8)             :: i, j, k, l, ij, kl
  real(DP)                :: energy, cml
  character(len=30), parameter :: Frmt = '(1x,4i5,2x,3f20.14)'

  cml = 0.0_DP
  ij = 0
  do i = 1, nb
      do j = 1, i
          ij = ij + 1
          kl = 0
          do k = 1, nb
              do l = 1, k
                  kl = kl + 1
                  if (ij >= kl) then
                       cml = cml + twordm(address(i,j,k,l)) * factor(i,j,k,l) * twoe(address(i,j,k,l))
                  endif
              enddo
          enddo
      enddo
  enddo

  energy = cml
 end function exactEnergy

  subroutine print4index(fi, nb)
  real(I8), intent(in)  :: fi(:)
  integer(I8)           :: i, j, k, l, ij, kl, nb

  real(DP), parameter :: tol = 1.0d-14
  character(len=30), parameter :: frmt = '(2x,4i5,2x,f4.1,e30.16)'

  write(*,*)
   ij = 0
   do i = 1,nb
     do j = 1,i
       ij = ij+1
       kl = 0
       do k = 1,nb
         do l = 1,k
           kl = kl + 1
             if(ij >= kl .and. abs(fi(address(i,j,k,l))) > tol) write(*,frmt) i,j,k,l,factor(i,j,k,l),fi(address(i,j,k,l))
         enddo
       enddo
     enddo
   enddo

  end subroutine print4index

end module TwoElectronModule
