module OneElectronModule
  use VarModule
  use DictionaryModule
  implicit none

  private

  type energyType
      real(DP) :: nuclearRepulsion
      real(DP) :: electronic
      real(DP) :: total
      real(DP) :: sz
      real(DP) :: szz
      real(DP) :: core
      real(DP) :: scf
      real(DP) :: electronNucelus
      real(DP) :: electronElectron
      real(DP) :: potential
      real(DP) :: kinetic
  end type

  interface print_energies
      module procedure print_energy_type
  end interface

  public :: energyType, print_energies

  interface readS
      module procedure readOverlap
  end interface
  interface readH
      module procedure readCoreHamiltonian
  end interface
  interface readT
      module procedure readKineticEnergyIntegrals
  end interface
  interface readXdip
      module procedure readXdipoleIntegrals
  end interface
  interface readYdip
      module procedure readYdipoleIntegrals
  end interface
  interface readZdip
      module procedure readZdipoleIntegrals
  end interface

  public :: readS, readH, readT, readXdip, readYdip, readZdip

  interface readMOs
      module procedure readHFMolecularOrbitals
  end interface
  interface readNOs
      module procedure readNaturalOrbitals
  end interface
  interface readHFOcc
      module procedure readHFOccupationNumbers
  end interface
  interface readOcc
      module procedure readOccupationNumbers
  end interface
  interface readEnergies
      module procedure readEnergiesCommonBlock
  end interface

  public :: readMOs, readNOs 
  public :: readOcc, readEnergies, readHFOcc

  contains

 subroutine readOverlap(S, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: S(:,:)
  real(DP), allocatable :: Stemp(:)
  integer(IK)           :: i, j, ij

  allocate(Stemp(size(S,1)*(size(S,1)+1)/2))

  call readFrom(dictionary, 'overlap', Stemp)

! compose the complete matrix
  ij = 0
  do i = 1, size(S, 1)
    do j = 1, i-1
      ij = ij + 1
      S(i,j) = Stemp(ij)
      S(j,i) = Stemp(ij)
    enddo
    ij = ij + 1
    S(i,i) = Stemp(ij)
  enddo

  deallocate(Stemp)
 end subroutine readOverlap

 subroutine readCoreHamiltonian(H, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: H(:,:)
  real(DP), allocatable :: Htemp(:)
  integer(IK)           :: i, j, ij

  allocate(Htemp(size(H,1)*(size(H,1)+1)/2))

  call readFrom(dictionary, 'core hamiltonian', Htemp)

! compose the complete matrix
  ij = 0
  do i = 1, size(H, 1)
    do j = 1, i-1
      ij = ij + 1
      H(i,j) = Htemp(ij)
      H(j,i) = Htemp(ij)
    enddo
    ij = ij + 1
    H(i,i) = Htemp(ij)
  enddo

  deallocate(Htemp)
 end subroutine readCoreHamiltonian

 subroutine readKineticEnergyIntegrals(T, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: T(:,:)
  real(DP), allocatable :: Ttemp(:)
  integer(IK)           :: i, j, ij

  allocate(Ttemp(size(T,1)*(size(T,1)+1)/2))

  call readFrom(dictionary, 'kinetic energy', Ttemp)

! compose the complete matrix
  ij = 0
  do i = 1, size(T, 1)
    do j = 1, i-1
      ij = ij + 1
      T(i,j) = Ttemp(ij)
      T(j,i) = Ttemp(ij)
    enddo
    ij = ij + 1
    T(i,i) = Ttemp(ij)
  enddo

  deallocate(Ttemp)
 end subroutine readKineticEnergyIntegrals

 subroutine readXdipoleIntegrals(X, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: X(:,:)
  real(DP), allocatable :: Xtemp(:)
  integer(IK)           :: i, j, ij

  allocate(Xtemp(size(X,1)*(size(X,1)+1)/2))

  call readFrom(dictionary, 'X dipole', Xtemp)

! compose the complete matrix
  ij = 0
  do i = 1, size(X, 1)
    do j = 1, i-1
      ij = ij + 1
      X(i,j) = Xtemp(ij)
      X(j,i) = Xtemp(ij)
    enddo
    ij = ij + 1
    X(i,i) = Xtemp(ij)
  enddo

  deallocate(Xtemp)
 end subroutine readXdipoleIntegrals

 subroutine readYdipoleIntegrals(Y, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Y(:,:)
  real(DP), allocatable :: Ytemp(:)
  integer(IK)           :: i, j, ij

  allocate(Ytemp(size(Y,1)*(size(Y,1)+1)/2))

  call readFrom(dictionary, 'Y dipole', Ytemp)

! compose the complete matrix
  ij = 0
  do i = 1, size(Y, 1)
    do j = 1, i-1
      ij = ij + 1
      Y(i,j) = Ytemp(ij)
      Y(j,i) = Ytemp(ij)
    enddo
    ij = ij + 1
    Y(i,i) = Ytemp(ij)
  enddo

  deallocate(Ytemp)
 end subroutine readYdipoleIntegrals

 subroutine readZdipoleIntegrals(Z, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Z(:,:)
  real(DP), allocatable :: Ztemp(:)
  integer(IK)           :: i, j, ij

  allocate(Ztemp(size(Z,1)*(size(Z,1)+1)/2))

  call readFrom(dictionary, 'Z dipole', Ztemp)

! compose the complete matrix
  ij = 0
  do i = 1, size(Z, 1)
    do j = 1, i-1
      ij = ij + 1
      Z(i,j) = Ztemp(ij)
      Z(j,i) = Ztemp(ij)
    enddo
    ij = ij + 1
    Z(i,i) = Ztemp(ij)
  enddo

  deallocate(Ztemp)
 end subroutine readZdipoleIntegrals

 subroutine readHFMolecularOrbitals(Vaomo, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Vaomo(:,:)

  call readFrom(dictionary, 'HF MOs', Vaomo)
 end subroutine readHFMolecularOrbitals

 subroutine readNaturalOrbitals(Vaono, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Vaono(:,:)

  call readFrom(dictionary, 'natural orbitals', Vaono)
 end subroutine readNaturalOrbitals

 subroutine readOccupationNumbers(Occ, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Occ(:)

  call readFrom(dictionary, 'occupation numbers', Occ)
 end subroutine readOccupationNumbers

 subroutine readHFOccupationNumbers(Occ, dictionary)
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Occ(:)

  call readFrom(dictionary, 'alpha energies', Occ)
 end subroutine readHFOccupationNumbers

 subroutine readEnergiesCommonBlock(self, dictionary)
  type(energyType)      :: self
  type(dictionaryType)  :: dictionary
  real(DP), allocatable :: E(:)

! extract the values form the common block and store them in the
! custom type energyType
!COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
!                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP

  allocate(E(115))
  call readFrom(dictionary, 'energies', E)
  self%nuclearRepulsion = E( 1)
  self%electronic       = E( 2)
  self%total            = E( 3)
  self%sz               = E( 4)
  self%szz              = E( 5)
  self%core             = E( 6)
  self%scf              = E( 7)
  self%electronNucelus  = E(11)
  self%electronElectron = E(12)
  self%potential        = E(13)
  self%kinetic          = E(14)
  deallocate(E)
 end subroutine readEnergiesCommonBlock

 subroutine print_energy_type(self)
  type(energyType) :: self

  write(*,'(/4x,a)') repeat('=',len(trim('Energies from dictionary file')))
  write(*,'(4x,a)') trim('Energies from dictionary file')
  write(*,'(4x,a/)') repeat('=',len(trim('Energies from dictionary file')))

  write(*,'(1x,a40," = ",f12.8)') 'Nuclear repulsion energy', self%nuclearRepulsion
  write(*,'(1x,a40," = ",f12.8)') 'Electronic energy', self%electronic
  write(*,'(1x,a40," = ",f12.8)') 'Total energy', self%total
  write(*,'(1x,a40," = ",f12.8)') 'Total spin Sz', self%sz
  write(*,'(1x,a40," = ",f12.8)') 'Square of spin ', self%szz
  write(*,'(1x,a40," = ",f12.8)') 'Frozen core energy', self%core
  write(*,'(1x,a40," = ",f12.8)') 'SCF total energy', self%scf
  write(*,'(1x,a40," = ",f12.8)') 'Electron-nucleus attraction', self%electronNucelus
  write(*,'(1x,a40," = ",f12.8)') 'Electron-electron repulsion', self%electronElectron
  write(*,'(1x,a40," = ",f12.8)') 'Potential energy', self%potential
  write(*,'(1x,a40," = ",f12.8)') 'Kinetic energy', self%kinetic
 end subroutine print_energy_type

end module OneElectronModule
