program read_integral
!===============================================================================
! documentation of this module can be found in dmft_doc.pdf file
!===============================================================================
  use varModule
  use dictionaryModule
  use OneElectronModule
  use TwoElectronModule
  use IOModule
  implicit none
  logical               :: oneepart, twoepart, readtwoeao, readtwoemo, readtwordm
  character(len=100)    :: buffer, filename
  character(len=100)    :: dictfile, twointaofile, twointmofile, twordmfile
  type(dictionaryType)  :: dictionary
  type(energyType)      :: energy
  real(DP), allocatable :: twointao(:), twointmo(:), twordm(:)
  real(DP), allocatable :: S(:,:), H(:,:), T(:,:), X(:,:), Y(:,:), Z(:,:)
  real(DP), allocatable :: Vaomo(:,:), Vaono(:,:), Occ(:), E(:), hfocc(:)
  integer(I8) :: nmo, nao, nt, lengthao, lengthmo, nprint
  real(DP)    :: eee

  ! input processing

  namelist /input/ nao, nmo, nprint, oneepart, twoepart, readtwoeao, readtwoemo, &
 &                 readtwordm, dictfile, twointmofile, twointaofile, twordmfile

  ! set defaults

  oneepart = .false.
  twoepart = .false.
  nao = 0
  nmo = 0
  nprint = 4

  ! get the input file name from command line

  call getarg(1, buffer)
  read(buffer, *) filename

  ! read the input file contents into the namelist /input/

  open(unit=111, file=trim(filename), status='old', form='formatted', delim='apostrophe')
  read(111, nml=input)

  ! stuff read from the dictionary file

  if (oneepart) then

    call new(dictionary, trim(dictfile))

    write(*,*) dictionary%irecst
    allocate(S(nao,nao), H(nao,nao), T(nao,nao), X(nao,nao), Y(nao,nao), Z(nao,nao))
    allocate(Vaomo(nao,nmo), Vaono(nao,nmo), Occ(nmo), E(115), hfocc(nmo))

    call readS(S, dictionary)
    call readH(H, dictionary)
    call readT(T, dictionary)
    call readXdip(X, dictionary)
    call readYdip(Y, dictionary)
    call readZdip(X, dictionary)
    call readMOs(Vaomo, dictionary)
!    call readNOs(Vaono, dictionary)
!    call readOcc(Occ, dictionary)
    call readHFOcc(hfocc, dictionary)
    call readEnergies(energy, dictionary)

    call matprint(S, 'Overlap Matrix in AO')
    call matprint(T, 'Kinetic Energy Matrix in AO')
    call matprint(H, 'Core Hamiltonian in AO')
    call matprint(X, 'X Dipole Matrix in AO')
    call matprint(Y, 'Y Dipole Matrix in AO')
    call matprint(Z, 'Z Dipole Matrix in AO')
    call matprint(Vaomo, 'HF Molecular Orbitals in AO')
!    call matprint(Vaono, 'Natural Orbital in AOs')
!    call matprint(Occ, 'Occupation Numbers')
    call matprint(hfocc, 'HF Occupation Numbers/Energies')

    call print_energies(energy)
    call delete(dictionary)
    deallocate(S, H, T, X, Y, Z, Vaomo, Vaono, Occ, E, hfocc)
  endif

  ! two-electron part

  if (twoepart) then

    nt = nmo*(nmo+1)/2
    lengthmo = nt*(nt+1)/2

    nt = nao*(nao+1)/2
    lengthao = nt*(nt+1)/2

    if (readtwoeao) then
      allocate(twointao(lengthao))
      call ReadInAO(twointao, twointaofile)
      call printHeader("First few two-electron integrals in AO")
      if (nao <= nprint) then
        call print4index(twointao, nao)
      else
        call print4index(twointao, nprint)
      endif
    endif

    if (readtwoemo) then
      allocate(twointmo(lengthmo))
      call ReadInMO(twointmo, twointmofile)
      call printHeader("First few two-electron integrals in MO")
      if (nmo <= nprint) then
        call print4index(twointmo, nmo)
      else
        call print4index(twointmo, nprint)
      endif
    endif

    if (readtwordm) then
      allocate(twordm(lengthmo))
      call ReadInAO(twordm, twordmfile)
      call printHeader("Few first elements of the two-particle denisty matrix")
      if (nmo <= nprint) then
        call print4index(twordm, nmo)
      else
        call print4index(twordm, nprint)
      endif
    endif

    if (readtwoemo .and. readtwordm) then
      eee = getTwoElEnergy(twointmo, twordm, nmo)
      write(*,"(a,f15.10)") "Two electron energy \sum_ijkl(\Gamma_ijkl * <ij|kl>) : ",0.5d00*eee
    endif

    if (allocated(twointao)) deallocate(twointao)
    if (allocated(twointmo)) deallocate(twointmo)
    if (allocated(twordm)) deallocate(twordm)

  endif

end program
