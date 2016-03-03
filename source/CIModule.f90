module CIModule
  use VarModule
  use DictionaryModule
  implicit none

  private

  public :: readCIVectorSize, readCIVector

contains

  subroutine readCIVectorSize(filename, nstates, nconfs)
! read CI eigenvector from NFT12 file
  integer(IK), parameter :: nft12 = 12

    integer(IPGamess),   intent(inout)   :: nstates, nconfs
    character(len=*), intent(in)      :: filename

    character(len=80) :: title, title1

    integer(IK) :: readStatus

    open(unit=nft12, file=trim(filename), status='old', form='unformatted')
    rewind(nft12)

    read(nft12) nstates,nconfs,title,title1

    write(*,*) nstates,nconfs,title,title1
    close(nft12)

  end subroutine readCIVectorSize

  subroutine readCIVector(vector, filename)
  integer(IK), parameter :: nft12 = 12
    real(DP), intent(inout) :: vector(:)
    character(len=*), intent(in)      :: filename

    integer(IK) :: readStatus

    logical :: ex, op
    character (len=11) :: nam, acc, seq, fnn
    integer :: irec, nr, length

    open(unit=nft12, file=trim(filename), status='old', form='unformatted')
    rewind(nft12)

    inquire(nft12, exist=ex, opened=op, name=nam, access=acc, &
    sequential=seq, form=fnn, recl=irec, nextrec=nr)
    write(*,*) 'recl = ', irec
    write(*,*) 'nextrec = ', nr
    read(nft12)
    read(nft12, iostat=readStatus) vector
    write(*,*) "read all the coefficients"
    close(nft12)

  end subroutine readCIVector
end module CIModule
