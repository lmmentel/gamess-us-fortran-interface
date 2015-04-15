module DictionaryModule
  use VarModule

  implicit none

  private

  type dictionaryType
      integer(IPgamess) :: idaf
      integer(IPgamess) :: irecln
      integer(IPgamess) :: irecst
      integer(IPgamess) :: ioda(950)
      integer(IPgamess) :: ifilen(950)
      integer(IPgamess) :: is
      integer(IPgamess) :: ipk
  end type

  public :: dictionaryType

  interface new
      module procedure newPrivate
  end interface
  interface delete
      module procedure deletePrivate
  end interface

  public :: new, delete

  interface readFrom
      module procedure ReadReals
      module procedure ReadRealMatrix
  end interface

  public :: readFrom

contains

 subroutine newPrivate(self, filename)
  type(dictionaryType) :: self
  character(*) :: filename

  self%idaf   = 10_IPgamess
  self%irecln = 4090_IPgamess

  open(unit=self%idaf,file=filename,access='direct',recl=8*self%irecln,form='unformatted')
  read(10,rec=1) self%irecst, self%ioda, self%ifilen, self%is, self%ipk

 end subroutine newPrivate

 subroutine deletePrivate(self)
  type(dictionaryType) :: self

  close(self%idaf)

 end subroutine deletePrivate

 subroutine ReadReals(self, section, array)
    type(dictionaryType), intent(in)  :: self
    character(*),         intent(in)  :: section
    real(DP),             intent(out) :: array(:)

!   Ordered alphabetically
    select case(trim(section))
      case ('alpha energies')
        call daread(self, array, int(size(array), kind=IPgamess), 17_IPgamess)
      case ('core hamiltonian')
        call daread(self, array, int(size(array), kind=IPgamess), 11_IPgamess)
      case ('energies')
        call daread(self, array, int(size(array), kind=IPgamess),  2_IPgamess)
      case ('kinetic energy')
        call daread(self, array, int(size(array), kind=IPgamess), 13_IPgamess)
      case ('occupation numbers')
        call daread(self, array, int(size(array), kind=IPgamess), 21_IPgamess)
      case ('overlap')
        call daread(self, array, int(size(array), kind=IPgamess), 12_IPgamess)
      case ('X dipole')
        call daread(self, array, int(size(array), kind=IPgamess), 95_IPgamess)
      case ('Y dipole')
        call daread(self, array, int(size(array), kind=IPgamess), 96_IPgamess)
      case ('Z dipole')
        call daread(self, array, int(size(array), kind=IPgamess), 97_IPgamess)
      case default
        write(stderr,*) 'No valid section for reading the dictionary file selected'
        stop
    end select

 end subroutine ReadReals

 subroutine ReadRealMatrix(self, section, matrix)
    type(DictionaryType), intent(in)  :: self
    character(*),         intent(in)  :: section
    real(DP),             intent(out) :: matrix(:,:)

    real(DP), allocatable :: vector(:)
    allocate(vector(size(matrix)))

    select case(trim(section))
      case ('guess')
        call daread(self, vector, int(size(matrix), kind=IPgamess), 265_IPgamess)
      case ('HF MOs')
        call daread(self, vector, int(size(matrix), kind=IPgamess),  15_IPgamess)
      case ('natural orbitals')
        call daread(self, vector, int(size(matrix), kind=IPgamess),  19_IPgamess)
      case ('salc')
        call daread(self, vector, int(size(matrix), kind=IPgamess),  44_IPgamess)
      case ('orthonormal symmetrized AOs')
        call daread(self, vector, int(size(matrix), kind=IPgamess),  45_IPgamess)
      case default
        write(stderr,*) 'No valid section for reading the dictionary file selected in -ReadRealMatrix-'
        stop
    end select

    matrix = reshape(vector, (/size(matrix,1), size(matrix, 2)/))

 end subroutine ReadRealMatrix

 subroutine daread(self, vector, length, nrec)
    type(DictionaryType), intent(in)  :: self
    integer(IPGamess),    intent(in)  :: length, nrec
    real(DP), dimension(length), intent(out) :: vector

    integer(IPGamess) :: n, ns, is, nsp, fi, lent, lenw

    n = self%ioda(nrec)

    if (n == -1) then
        write(*,'(/1X,"ERROR *** ATTEMPTING A BOGUS READ OF A DAF RECORD."/ &
   &              1X,"RECORD NUMBER",I5," OF LENGTH",I10, " WAS NEVER PREVIOUSLY WRITTEN.")') nrec, length
        stop "exiting"
    endif

    if (length <= 0 ) then
        write(*,'(/1X,"ERROR *** ATTEMPTING A BOGUS READ OF A DAF RECORD."/ &
   &              1X,"RECORD NUMBER",I5," OF LENGTH",I10," HAS NO LENGTH.")') nrec, length
        stop "exiting"
    endif

    if (length > self%ifilen(nrec)) then
        write(*,'(/1X,"ERROR *** ATTEMPTING A BOGUS READ OF A DAF RECORD."/      &
   &              1X,"ATTEMPTING TO READ",I10," WORDS FROM RECORD NUMBER",I5/    &
   &              1X,"BUT THIS RECORD WAS PREVIOUSLY WRITTEN WITH ONLY",         &
   &              I10," WORDS.")') length, nrec, self%ifilen(nrec)
        stop "exiting"
    endif

    IS = -self%irecln + 1
    NS = N
    LENT = length
    do while (lent >= 1)
        IS = IS + self%irecln
        fi = IS + LENT - 1
        if ((fi-IS+1) .GT. self%irecln) fi = IS + self%irecln - 1
        NSP = NS
        LENW = fi - IS + 1
        read(unit=self%idaf, rec=nsp) vector
        LENT = LENT - self%irecln
        NS = NS + 1
        N = NS
    enddo
 end subroutine daread

end module dictionaryModule
