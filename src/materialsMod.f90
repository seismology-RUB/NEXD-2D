!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2019 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
!
!   This file is part of NEXD 2D.
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful, but
!   WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------
module materialsMod

    use parameterMod
    use constantsMod
    use errorMessage
    use linearSystemMod

    implicit none

    interface setupMaterials
        module procedure setupMaterialsElastic
        module procedure setupMaterialsPoroelastic
    end interface setupMaterials

    type materialVar
        real(kind=custom_real) :: vp         !P-wave velocity
        real(kind=custom_real) :: vs         !S-wave velocity
        real(kind=custom_real) :: lambda     !1st Lamé parameter
        real(kind=custom_real) :: mu         !2nd Lamé parameter
        real(kind=custom_real) :: rho        !Density
        real(kind=custom_real) :: imp_vp     !Impedance with respect to vp
        real(kind=custom_real) :: imp_vs     !Impedance with respect to vs
        real(kind=custom_real) :: qp         !Qualityfactor for P-wave speed
        real(kind=custom_real) :: qs         !Qualityfactor for S-wave speed
        real(kind=CUSTOM_REAL), dimension(nMB) :: ylambda   !Anelastic coefficient regarding the 1st Lamé parameter
        real(kind=CUSTOM_REAL), dimension(nMB) :: ymu       !Anelastic coefficient regarding the 2st Lamé parameter
        real(kind=CUSTOM_REAL), dimension(nMB) :: omegaL    !Relaxation frequencies of the Maxwellbodys
    end type materialVar

    type porous_material
        integer :: nmat
        !parameters to read in
        integer :: typ
        real(kind=custom_real) :: rhos        !solid density
        real(kind=custom_real) :: lambda      !1st Lamé parameter
        real(kind=custom_real) :: my          !2nd Lamé parameter / shear modulus
        real(kind=custom_real) :: phi         !porosity
        real(kind=custom_real) :: kappa       !permeabilty
        real(kind=custom_real) :: biot        !Biot coefficient
        real(kind=custom_real) :: invT        !inverse Tortuosity
        real(kind=custom_real) :: invN        !inverse of Biot modulus
        real(kind=custom_real) :: rho1        !fluid 1 density
        real(kind=custom_real) :: S1          !fluid 1 saturation
        real(kind=custom_real) :: K1          !fluid 1 bulk modulus
        real(kind=custom_real) :: ny1         !fluid 1 viscosity
        real(kind=custom_real) :: rho2        !fluid 2 density
        real(kind=custom_real) :: S2          !fluid 2 saturation
        real(kind=custom_real) :: K2          !fluid 2 bulk modulus
        real(kind=custom_real) :: ny2         !fluid 2 viscosity

        real(kind=custom_real) :: fitting_m   !fitting parameter for van Genuchten model
        real(kind=custom_real) :: fitting_n   !fitting parameter for van Genuchten model
        real(kind=custom_real) :: fitting_chi !fitting parameter for van Genuchten model
        real(kind=custom_real) :: Sr1         !fluid 1 residual saturation
        real(kind=custom_real) :: Sr2         !fluid 2 residual saturation
        
        !parameters to calculate
        real(kind=custom_real) :: krel1       !fluid 1 relative permeability
        real(kind=custom_real) :: krel2       !fluid 2 relative permeability
        real(kind=custom_real) :: dpcdS1      !capillary pressure term
        real(kind=custom_real) :: M           !modulus

        real(kind=custom_real) :: S1eff       !fluid 1 effective saturation
        real(kind=custom_real) :: S2eff       !fluid 2 effective saturation

        real(kind=custom_real), dimension(:,:), pointer :: A,B  !matrix A,B
        real(kind=custom_real), dimension(:,:), pointer :: AP,AM,BP,BM
        real(kind=custom_real), dimension(:,:), pointer :: E  !matrix E

        real(kind=custom_real) :: vmax        !maximum velocity calculated from matrix A
        real(kind=custom_real) :: vmin        !minimum velocity calculated from matrix A
    end type porous_material

    type materialIndizes
        integer, dimension(:), allocatable :: type           !Index of the material-type that is used in each element
        integer, dimension(:), allocatable :: pml            !Index on wether an element belongs to a pml or not
    end type materialIndizes

!    type attenuationVar
!        real(kind=custom_real) :: vp                        !P-wave velocity (unrelaxed)
!        real(kind=custom_real) :: vs                        !S-wave velocity (unrelaxed)
!        real(kind=custom_real) :: lambda                    !1st Lamé parameter (unrelaxed)
!        real(kind=custom_real) :: mu                        !2nd Lamé parameter (unrelaxed)
!
!    end type attenuationVar
    contains

    subroutine setupMaterialsElastic(par, mat, matInd, nelem, iregion, errmsg)
        !This subroutine initializes the material variables for the different cases (elastic, attenuation)
        !input
        type(parameterVar) :: par
        type(error_message) :: errmsg
        integer :: nelem
        integer, dimension(:), allocatable :: iregion
        !in/out
        type(materialIndizes) :: matInd
        type(materialVar), dimension(:), allocatable :: mat
        !local
        character(len=14) :: myname = "setupMaterials"

        call addTrace(errmsg, myname)

        ! materials values
        call readMaterialProperties(par, mat, 19, iregion, trim('cubit/matprop'), errmsg)

        ! materials array
        call readMaterialFile(par, matInd, nelem, 19, trim('cubit/mat'), errmsg)

        !Attenuation
        if (par%attenuation) then
            call prepareAnelasticCalc(par, mat, errmsg)
        end if
    end subroutine setupMaterialsElastic

    subroutine setupMaterialsPoroelastic(par, poromat, matInd, nelem, errmsg)
        !This subroutine initializes the material variables for poroelastic materials
        !input
        type(parameterVar) :: par
        type(error_message) :: errmsg
        integer :: nelem
        !in/out
        type(materialIndizes) :: matInd
        type(porous_material), dimension(:), allocatable :: poromat
        !local
        character(len=14) :: myname = "setupMaterials"

        call addTrace(errmsg, myname)

        ! materials values
        if (par%extmatprop) then
            call readPorousMaterialProperties(par, poromat, 19, trim(par%extmatpropfilename), errmsg)
        else
            call readPorousMaterialProperties(par, poromat, 19, trim('cubit/matprop'), errmsg)
        end if

        ! materials array
        call readMaterialFile(par, matInd, nelem, 19, trim('cubit/mat'), errmsg)

        !Attenuation
        !if (par%attenuation) then
        !    call prepareAnelasticCalc(par, mat, att, errmsg)
        !end if
    end subroutine setupMaterialsPoroelastic

    subroutine readMaterialProperties(par, mat, lu, iregion, filename, errmsg)
        !Subroutine to read properties of the materials from the cubit/trelis-generated file
        !input
        type(parameterVar) :: par
        type (error_message) :: errmsg
        integer, intent(in) :: lu
        character(len=*), intent(in) :: filename
        !in/output
        type(materialVar), dimension(:), allocatable :: mat
        integer, dimension(:), allocatable :: iregion              !Some parameter for the external model. Not sure if is needed anymore... TM TM
        !local
        integer :: i
        integer :: dummy
        integer :: ios
        character(len=22) :: myname = "readMaterialProperties"

        call addTrace(errmsg, myname)

        open(unit=lu,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif

        read(lu,*) par%matn
        allocate(mat(par%matn))
        allocate(iregion(par%matn))

        do i = 1, par%matn
           !vp,vs,rho,qp,qs
           read(lu,*) iregion(i), dummy, mat(i)%vp, mat(i)%vs, mat(i)%rho, mat(i)%qp ,mat(i)%qs
           !compute mu and lambda and the impedances
           mat(i)%mu = mat(i)%vs**2 * mat(i)%rho
           mat(i)%lambda = mat(i)%vp**2 * mat(i)%rho - 2*mat(i)%mu
           mat(i)%imp_vp = mat(i)%vp * mat(i)%rho
           mat(i)%imp_vs = mat(i)%vs * mat(i)%rho
        end do
        close(lu)
        deallocate(iregion)
    end subroutine readMaterialProperties

    subroutine readPorousMaterialProperties(par, poromat, lu, filename, errmsg)
        !in/output
        type(porous_material), dimension(:), allocatable :: poromat
        integer, dimension(:), allocatable :: iregion              !Some parameter for the external model. Not sure if is needed anymore... TM TM
        !input
        type (parameterVar), intent(in) :: par
        type (error_message) :: errmsg
        integer :: lu
        character(len=*) :: filename
        !local
        integer :: i!,iii,iij
        integer :: ios
        integer :: mattypes
        character(len=27) :: myname = 'readPorousMaterialProperies'
        character(len=255) :: text
        !character(len=3) :: str
        !character(len=10) :: mati
        real(kind=custom_real) :: rho,rhoast,M1,M2,Mti,Mti1,Mti2
        real(kind=custom_real) :: lambdaast1_inv,lambdaast2_inv,lambdati11_inv,lambdati12_inv,lambdati21_inv,lambdati22_inv
        real(kind=custom_real) :: vmaxA, vmaxB, vminA, vminB

        call addTrace(errmsg,myname)
        open(unit=lu, file = trim(filename), status = 'old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            return
        endif

        text = 'none'
        if (par%extmatprop) then
            do while(text /= 'BEGIN')
                read(lu,'(a)') text
            enddo
        endif
        read(lu,*) mattypes

        !!check whether material specified in 'mat_bap' is available in 'porousmaterials' or not
        !do i = 1, size(mat%i)
        !    if (mat%i(i) > mattypes) then
        !        write (mati,'(i10)') mat%i(i)
        !        call add(errmsg,2,"Material "//trim(adjustl(mati))//&
        !         " specified in 'mat_bap' is not available in 'porousmaterials'.",myname)
        !        return
        !    endif
        !enddo

        allocate(poromat(mattypes))
        allocate(iregion(mattypes))

        do i = 1, mattypes
            if (par%extmatprop) then
                read(lu,*) poromat(i)%typ, poromat(i)%rhos, poromat(i)%lambda, poromat(i)%my, poromat(i)%phi, poromat(i)%kappa,&
                 poromat(i)%biot, poromat(i)%invT, poromat(i)%invN, poromat(i)%rho1, poromat(i)%S1, poromat(i)%K1,&
                 poromat(i)%ny1, poromat(i)%rho2, poromat(i)%S2, poromat(i)%K2, poromat(i)%ny2,&
                 poromat(i)%fitting_n, poromat(i)%fitting_chi, poromat(i)%Sr1, poromat(i)%Sr2
            else
                read(lu,*) iregion(i), poromat(i)%typ, poromat(i)%rhos, poromat(i)%lambda, poromat(i)%my, poromat(i)%phi, poromat(i)%kappa,&
                 poromat(i)%biot, poromat(i)%invT, poromat(i)%invN, poromat(i)%rho1, poromat(i)%S1, poromat(i)%K1,&
                 poromat(i)%ny1, poromat(i)%rho2, poromat(i)%S2, poromat(i)%K2, poromat(i)%ny2,&
                 poromat(i)%fitting_n, poromat(i)%fitting_chi, poromat(i)%Sr1, poromat(i)%Sr2
            endif

            if (par%fluidn == 1) then
                poromat(i)%S1 = 1.0
                poromat(i)%S2 = 0.0
            endif
            if ((poromat(i)%S1 - 1.0) < epsilon(poromat(i)%S1)) then
                !It will be required to have rho2=rho1 if one fluid
                !vanishes.
                poromat(i)%rho2 = poromat(i)%rho1
                !The second fluid must not have a residual saturation.
                poromat(i)%Sr2 = 0.
            endif

            if ((poromat(i)%S1 - 1.0) < epsilon(poromat(i)%S1) .and. poromat(i)%S2 < epsilon(poromat(i)%S2)) then
                !Just to avoid a division by zero. All terms with K2 will become zero anyway in this case.
                poromat(i)%K2 = eps
            elseif (poromat(i)%S1 < epsilon(poromat(i)%S1) .and. (poromat(i)%S2 - 1.0) < epsilon(poromat(i)%S2)) then
                !Just to avoid a division by zero. All terms with K1 will become zero anyway in this case.
                poromat(i)%K1 = eps
            endif

            !check for correct input
            if ((poromat(i)%phi < 0.) .or. (poromat(i)%phi > 1.)) call add(errmsg,2,'phi must be between 0 and 1',myname)
            if (par%fluidn == 2 .and. (poromat(i)%S1 + poromat(i)%S2 - 1.0) > epsilon(poromat(i)%S1)) call add(errmsg,2,'S1+S2 must be 1',myname)
            if ((poromat(i)%biot < 0.) .or. (poromat(i)%biot > 1.)) call add(errmsg,2,'b must be between 0 and 1',myname)
            if (poromat(i)%phi > epsilon(poromat(i)%phi)) then
                if (par%calculate_tortuosity) then
                    if (poromat(i)%invT > 1. .or. poromat(i)%invT < 0.) call add(errmsg,2,'r must be between 0 and 1',myname)
                else
                    if (poromat(i)%invT > 1. .or. poromat(i)%invT < 0.) call add(errmsg,2,'1/T must be between 0 and 1',myname)
                endif
            endif
            if (poromat(i)%S1 > 1-poromat(i)%Sr2) call add(errmsg,2,'S1 <= 1-Sr2 required',myname)
            if (poromat(i)%S2 > 1-poromat(i)%Sr1) call add(errmsg,2,'S2 <= 1-Sr1 required',myname)
            if (poromat(i)%S1 < poromat(i)%Sr1) call add(errmsg,2,'S1 >= Sr1 required',myname)
            if (poromat(i)%S2 < poromat(i)%Sr2) call add(errmsg,2,'S2 >= Sr2 required',myname)

            if (poromat(i)%phi < epsilon(poromat(i)%phi)) then
                !When there is no porespace, the skeleton and the
                !material building up the skeleton have the same bulk
                !modulus. This means, that b must be 0.
                poromat(i)%biot = 0.
            elseif (poromat(i)%typ == 3 .or. poromat(i)%typ == 4 .or. poromat(i)%typ == 8) then
                !here: K^d = poromat(i)%lambda and K_s = poromat(i)%invN
                poromat(i)%biot = 1 - poromat(i)%lambda / poromat(i)%invN
            elseif (poromat(i)%typ == 5 .or. poromat(i)%typ == 6 .or. poromat(i)%typ == 9) then
                !here: K^d = poromat%lambda
                poromat(i)%biot = (poromat(i)%phi+1)/2 - &
                 sqrt((poromat(i)%phi+1)**2/4 - poromat(i)%phi - poromat(i)%lambda*poromat(i)%invN)
            endif

            if (poromat(i)%typ == 3 .or. poromat(i)%typ == 4) then
                !here: K_s = poromat(i)%invN before calculation
                poromat(i)%invN = (poromat(i)%biot-poromat(i)%phi)/poromat(i)%invN
                !now 1/N = poromat(i)%invN!
            endif

            if (par%calculate_tortuosity) then
                !here: r = poromat(i)%invT before calculation
                poromat(i)%invT = 1 / (1 - poromat(i)%invT * ( 1 - 1 / poromat(i)%phi ))
            endif

            poromat(i)%krel1 = 1.
            poromat(i)%krel2 = 1.
            poromat(i)%dpcdS1 = 1.
            if (par%fluidn == 2 .and. .not. poromat(i)%S2 < epsilon(poromat(i)%S2)) then
                poromat(i)%S1eff = (poromat(i)%S1-poromat(i)%Sr1)/(1-poromat(i)%Sr2-poromat(i)%Sr1)
                poromat(i)%S2eff = (poromat(i)%S2-poromat(i)%Sr2)/(1-poromat(i)%Sr1-poromat(i)%Sr2)
                select case (poromat(i)%typ)
                    case (1,3,5) ! van Genuchten model
                        !calculate more properties using van Genuchten (1980) model
                        poromat(i)%fitting_m = 1 - 1/poromat(i)%fitting_n
                        !calculate relative permeability
                        poromat(i)%krel1 = sqrt(poromat(i)%S1eff)*(1-(1-poromat(i)%S1eff**(1/poromat(i)%fitting_m))**&
                         poromat(i)%fitting_m)**2
                        poromat(i)%krel2 = sqrt(1-poromat(i)%S1eff)*(1-poromat(i)%S1eff**(1/poromat(i)%fitting_m))**&
                         (2*poromat(i)%fitting_m)
                        !calculate capillary pressure influence
                        poromat(i)%dpcdS1 = - ((poromat(i)%rho1*g)/((poromat(i)%fitting_n-1)*poromat(i)%fitting_chi)) * &
                         ((poromat(i)%S1**(-1/poromat(i)%fitting_m)-1)**(1/poromat(i)%fitting_n-1)) * &
                         (poromat(i)%S1**(-(1/poromat(i)%fitting_m+1)))
                    case (2,4,6) ! Brooks & Corey model
                        !here p_b = fitting_n and \lambda_{BC} = fitting_chi
                        !calculate relative permeability
                        poromat(i)%krel1 = poromat(i)%S1eff**((2+3*poromat(i)%fitting_chi)/poromat(i)%fitting_chi)
                        poromat(i)%krel2 = (1-poromat(i)%S1eff)**2 * &
                         (1-poromat(i)%S1eff**((2+poromat(i)%fitting_chi)/poromat(i)%fitting_chi))
                         !calculate capillary pressure influence
                        poromat(i)%dpcdS1 = - poromat(i)%fitting_n/(poromat(i)%fitting_chi*(poromat(i)%Sr2-poromat(i)%Sr1))*&
                         poromat(i)%S1**(-(1+poromat(i)%fitting_chi)/poromat(i)%fitting_chi)
                    case (7,8,9) ! Douglas Jr. et al. model
                        !here A = fitting_n
                        !calculate relative permeability
                        poromat(i)%krel1 = (1-poromat(i)%S1/(1-poromat(i)%Sr2))**2
                        poromat(i)%krel2 = ( (poromat(i)%S1 - poromat(i)%Sr1) / (1-poromat(i)%Sr1) )**2
                        !calculate capillary pressure influence
                        poromat(i)%dpcdS1 = - 2 * poromat(i)%fitting_n * ((poromat(i)%S1-poromat(i)%Sr1)**(-3) + &
                         poromat(i)%Sr2**2 / ((1-poromat(i)%S1)**2 * (1-poromat(i)%Sr1-poromat(i)%Sr2)**2 ) )
                end select
            endif

            !calculate M
            poromat(i)%M = 1/(poromat(i)%invN + poromat(i)%phi*(poromat(i)%S1/poromat(i)%K1+&
             poromat(i)%S2/poromat(i)%K2) - ( poromat(i)%S1*poromat(i)%S2/(poromat(i)%K1*poromat(i)%K2)*poromat(i)%phi&
             + poromat(i)%S1*poromat(i)%S2*poromat(i)%invN*(poromat(i)%S1/poromat(i)%K2+poromat(i)%S2/poromat(i)%K1))*&
             poromat(i)%dpcdS1)

            if (poromat(i)%typ == 3 .or. poromat(i)%typ == 4 .or. poromat(i)%typ == 5 .or. poromat(i)%typ == 6 &
             .or. poromat(i)%typ == 8 .or. poromat(i)%typ == 9) then
                !here: on the rhs K^d = poromat(i)%lambda
                poromat(i)%lambda = poromat(i)%lambda - 2./3. * poromat(i)%my + poromat(i)%M*(poromat(i)%biot**2) * &
                 (1-poromat(i)%S1*poromat(i)%S2*(poromat(i)%S1/poromat(i)%K2+poromat(i)%S2/poromat(i)%K1)*&
                 poromat(i)%dpcdS1)
            endif

            rho    = (1-poromat(i)%phi)*poromat(i)%rhos + poromat(i)%phi*poromat(i)%S1*poromat(i)%rho1 &
                                                        + poromat(i)%phi*poromat(i)%S2*poromat(i)%rho2
            rhoast = rho
            if (poromat(i)%phi > epsilon(poromat(i)%phi)) then
                rhoast = rhoast - poromat(i)%invT*poromat(i)%phi*poromat(i)%S1*poromat(i)%rho1 &
                                - poromat(i)%invT*poromat(i)%phi*poromat(i)%S2*poromat(i)%rho2
            endif

            lambdaast1_inv = (poromat(i)%phi * poromat(i)%S1 * poromat(i)%ny1) / (poromat(i)%kappa * poromat(i)%krel1) !/C1
            lambdaast2_inv = (poromat(i)%phi * poromat(i)%S2 * poromat(i)%ny2) / (poromat(i)%kappa * poromat(i)%krel2) !/C2
            lambdati11_inv = (poromat(i)%phi * poromat(i)%S1 * (poromat(i)%invT**2-poromat(i)%invT) * lambdaast1_inv)/rhoast
            lambdati12_inv = (poromat(i)%phi * poromat(i)%S2 * (poromat(i)%invT**2-poromat(i)%invT) * lambdaast2_inv)/rhoast
            lambdati21_inv = lambdaast1_inv * poromat(i)%invT / poromat(i)%rho1 + lambdati11_inv
            lambdati22_inv = lambdaast2_inv * poromat(i)%invT / poromat(i)%rho2 + lambdati12_inv

            M1   = poromat(i)%M*(1-poromat(i)%S1*poromat(i)%S2/poromat(i)%K1*poromat(i)%dpcdS1)
            M2   = poromat(i)%M*(1-poromat(i)%S1*poromat(i)%S2/poromat(i)%K2*poromat(i)%dpcdS1)
            Mti  = poromat(i)%M*(1+poromat(i)%S1*poromat(i)%S2*poromat(i)%invN/poromat(i)%phi*poromat(i)%dpcdS1)
            Mti1 = poromat(i)%M* &
             (1-(poromat(i)%S1**2 * poromat(i)%invN/poromat(i)%phi + poromat(i)%S1/poromat(i)%K1) * poromat(i)%dpcdS1)
            Mti2 = poromat(i)%M* &
             (1-(poromat(i)%S2**2 * poromat(i)%invN/poromat(i)%phi + poromat(i)%S2/poromat(i)%K2) * poromat(i)%dpcdS1)

            allocate(poromat(i)%A(5+3*par%fluidn,5+3*par%fluidn))
            allocate(poromat(i)%B(5+3*par%fluidn,5+3*par%fluidn))
            allocate(poromat(i)%E(5+3*par%fluidn,5+3*par%fluidn))

            poromat(i)%A( :, :) = 0.
            poromat(i)%A( 1, 4) = -(poromat(i)%lambda + 2*poromat(i)%my)
            poromat(i)%A( 2, 4) = -poromat(i)%lambda
            poromat(i)%A( 3, 5) = -poromat(i)%my
            poromat(i)%A( 4, 1) = -1/rhoast
            poromat(i)%A( 5, 3) = -1/rhoast
            poromat(i)%B( :, :) = 0.
            poromat(i)%B( 1, 5) = -poromat(i)%lambda
            poromat(i)%B( 2, 5) = -(2*poromat(i)%my + poromat(i)%lambda)
            poromat(i)%B( 3, 4) = -poromat(i)%my
            poromat(i)%B( 4, 3) = -1/rhoast
            poromat(i)%B( 5, 2) = -1/rhoast
            poromat(i)%E( :, :) = 0.

            if (poromat(i)%phi > epsilon(poromat(i)%phi)) then
                poromat(i)%A( 1, 4) =  poromat(i)%A( 1, 4) &
                    + poromat(i)%biot*poromat(i)%phi*(poromat(i)%S1*M2 + poromat(i)%S2*M1)
                poromat(i)%A( 1, 7) = -poromat(i)%phi*poromat(i)%S1*M2*poromat(i)%biot
                poromat(i)%A( 2, 4) =  poromat(i)%A( 2, 4) &
                    + poromat(i)%biot*poromat(i)%phi*(poromat(i)%S1*M2 + poromat(i)%S2*M1)
                poromat(i)%A( 2, 7) = -poromat(i)%phi*poromat(i)%S1*M2*poromat(i)%biot
                poromat(i)%A( 4, 6) = -poromat(i)%phi*poromat(i)%S1*poromat(i)%invT/rhoast
                poromat(i)%A( 6, 4) =  M2*poromat(i)%biot - poromat(i)%phi*(poromat(i)%S1*Mti2 + poromat(i)%S2*Mti)
                poromat(i)%A( 6, 7) =  poromat(i)%phi*poromat(i)%S1*Mti2
                poromat(i)%A( 7, 1) =  (poromat(i)%invT-1)/rhoast
                poromat(i)%A( 7, 6) =  poromat(i)%invT/poromat(i)%rho1 &
                    + poromat(i)%phi*poromat(i)%S1*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                poromat(i)%A( 8, 3) =  (poromat(i)%invT-1)/rhoast
                if (par%fluidn == 2) then
                    poromat(i)%A( 1,10) = -poromat(i)%phi*poromat(i)%S2*M1*poromat(i)%biot
                    poromat(i)%A( 2,10) = -poromat(i)%phi*poromat(i)%S2*M1*poromat(i)%biot
                    poromat(i)%A( 4, 9) = -poromat(i)%phi*poromat(i)%S2*poromat(i)%invT/rhoast
                    poromat(i)%A( 6,10) =  poromat(i)%phi*poromat(i)%S2*Mti
                    poromat(i)%A( 7, 9) =  poromat(i)%phi*poromat(i)%S2*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                    poromat(i)%A( 9, 4) =  M1*poromat(i)%biot - poromat(i)%phi*(poromat(i)%S1*Mti + poromat(i)%S2*Mti1)
                    poromat(i)%A( 9, 7) =  poromat(i)%phi*poromat(i)%S1*Mti
                    poromat(i)%A( 9,10) =  poromat(i)%phi*poromat(i)%S2*Mti1
                    poromat(i)%A(10, 1) =  (poromat(i)%invT-1)/rhoast
                    poromat(i)%A(10, 6) =  poromat(i)%phi*poromat(i)%S1*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                    poromat(i)%A(10, 9) =  poromat(i)%invT/poromat(i)%rho2 &
                        + poromat(i)%phi*poromat(i)%S2*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                    poromat(i)%A(11, 3) =  (poromat(i)%invT-1)/rhoast
                endif
                poromat(i)%B( 1, 5) =  poromat(i)%B( 1, 5) &
                    + poromat(i)%biot*poromat(i)%phi*(poromat(i)%S1*M2 + poromat(i)%S2*M1)
                poromat(i)%B( 1, 8) = -poromat(i)%phi*poromat(i)%S1*M2*poromat(i)%biot
                poromat(i)%B( 2, 5) =  poromat(i)%B( 2, 5) &
                    + poromat(i)%biot*poromat(i)%phi*(poromat(i)%S1*M2 + poromat(i)%S2*M1)
                poromat(i)%B( 2, 8) = -poromat(i)%phi*poromat(i)%S1*M2*poromat(i)%biot
                poromat(i)%B( 5, 6) = -poromat(i)%phi*poromat(i)%S1*poromat(i)%invT/rhoast
                poromat(i)%B( 6, 5) =  M2*poromat(i)%biot - poromat(i)%phi*(poromat(i)%S1*Mti2 + poromat(i)%S2*Mti)
                poromat(i)%B( 6, 8) =  poromat(i)%phi*poromat(i)%S1*Mti2
                poromat(i)%B( 7, 3) =  (poromat(i)%invT-1)/rhoast
                poromat(i)%B( 8, 2) =  (poromat(i)%invT-1)/rhoast
                poromat(i)%B( 8, 6) =  poromat(i)%invT/poromat(i)%rho1 &
                    + poromat(i)%phi*poromat(i)%S1*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                if (par%fluidn == 2) then
                    poromat(i)%B( 1,11) = -poromat(i)%phi*poromat(i)%S2*M1*poromat(i)%biot
                    poromat(i)%B( 2,11) = -poromat(i)%phi*poromat(i)%S2*M1*poromat(i)%biot
                    poromat(i)%B( 5, 9) = -poromat(i)%phi*poromat(i)%S2*poromat(i)%invT/rhoast
                    poromat(i)%B( 6,11) =  poromat(i)%phi*poromat(i)%S2*Mti
                    poromat(i)%B( 8, 9) =  poromat(i)%phi*poromat(i)%S2*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                    poromat(i)%B( 9, 5) =  M1*poromat(i)%biot - poromat(i)%phi*(poromat(i)%S1*Mti + poromat(i)%S2*Mti1)
                    poromat(i)%B( 9, 8) =  poromat(i)%phi*poromat(i)%S1*Mti
                    poromat(i)%B( 9,11) =  poromat(i)%phi*poromat(i)%S2*Mti1
                    poromat(i)%B(10, 3) =  (poromat(i)%invT-1)/rhoast
                    poromat(i)%B(11, 2) =  (poromat(i)%invT-1)/rhoast
                    poromat(i)%B(11, 6) =  poromat(i)%phi*poromat(i)%S1*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                    poromat(i)%B(11, 9) =  poromat(i)%invT/poromat(i)%rho2 &
                        + poromat(i)%phi*poromat(i)%S2*(poromat(i)%invT**2-poromat(i)%invT)/rhoast
                endif
                poromat(i)%E( 4, 7) =  (poromat(i)%phi*poromat(i)%S1*poromat(i)%invT*lambdaast1_inv) / rhoast
                poromat(i)%E( 5, 8) =  (poromat(i)%phi*poromat(i)%S1*poromat(i)%invT*lambdaast1_inv) / rhoast
                poromat(i)%E( 7, 7) =  - lambdati21_inv
                poromat(i)%E( 8, 8) =  - lambdati21_inv
                poromat(i)%E( 4, 4) = -poromat(i)%E( 4, 7)
                poromat(i)%E( 5, 5) = -poromat(i)%E( 5, 8)
                poromat(i)%E( 7, 4) = -poromat(i)%E( 7, 7)
                poromat(i)%E( 8, 5) = -poromat(i)%E( 8, 8)
                if (par%fluidn == 2) then
                    poromat(i)%E( 4,10) = (poromat(i)%phi*poromat(i)%S2*poromat(i)%invT*lambdaast2_inv) / rhoast
                    poromat(i)%E( 5,11) = (poromat(i)%phi*poromat(i)%S2*poromat(i)%invT*lambdaast2_inv) / rhoast
                    poromat(i)%E( 7,10) = - lambdati12_inv
                    poromat(i)%E( 8,11) = - lambdati12_inv
                    poromat(i)%E(10, 7) = - lambdati11_inv
                    poromat(i)%E(10,10) = - lambdati22_inv
                    poromat(i)%E(11, 8) = - lambdati11_inv
                    poromat(i)%E(11,11) = - lambdati22_inv
                    poromat(i)%E( 4, 4) = -poromat(i)%E( 4, 7)-poromat(i)%E( 4,10)
                    poromat(i)%E( 5, 5) = -poromat(i)%E( 5, 8)-poromat(i)%E( 5,11)
                    poromat(i)%E( 7, 4) = -poromat(i)%E( 7, 7)-poromat(i)%E( 7,10)
                    poromat(i)%E( 8, 5) = -poromat(i)%E( 8, 8)-poromat(i)%E( 8,11)
                    poromat(i)%E(10, 4) = -poromat(i)%E(10, 7)-poromat(i)%E(10,10)
                    poromat(i)%E(11, 5) = -poromat(i)%E(11, 8)-poromat(i)%E(11,11)
                endif
            else
                poromat(i)%A( 7, 1) =     -1/rhoast
                poromat(i)%A( 8, 3) =     -1/rhoast
                if (par%fluidn == 2) then
                    poromat(i)%A(10, 1) = -1/rhoast
                    poromat(i)%A(11, 3) = -1/rhoast
                endif
                poromat(i)%B( 7, 3) =     -1/rhoast
                poromat(i)%B( 8, 2) =     -1/rhoast
                if (par%fluidn == 2) then
                    poromat(i)%B(10, 3) = -1/rhoast
                    poromat(i)%B(11, 2) = -1/rhoast
                endif
            endif

            !write (*,*) "A",i
            !do iii = 1,5+3*par%fluidn
            !    write (*,"(100g15.5)") (poromat(i)%A(iii,iij),iij=1,5+3*par%fluidn)
            !enddo
            !write (*,*) "B",i
            !do iii = 1,5+3*par%fluidn
            !    write (*,"(100g15.5)") (poromat(i)%B(iii,iij),iij=1,5+3*par%fluidn)
            !enddo
            !write (*,*) "E",i
            !do iii = 1,5+3*par%fluidn
            !    write (*,"(100g15.5)") (poromat(i)%E(iii,iij),iij=1,5+3*par%fluidn)
            !enddo

            allocate(poromat(i)%AP(5+3*par%fluidn,5+3*par%fluidn))
            allocate(poromat(i)%AM(5+3*par%fluidn,5+3*par%fluidn))
            call calcAPAM(poromat(i)%A(:,:),poromat(i)%AP(:,:),poromat(i)%AM(:,:),vmaxA,vminA)
            allocate(poromat(i)%BP(5+3*par%fluidn,5+3*par%fluidn))
            allocate(poromat(i)%BM(5+3*par%fluidn,5+3*par%fluidn))
            call calcAPAM(poromat(i)%B(:,:),poromat(i)%BP(:,:),poromat(i)%BM(:,:),vmaxB,vminB)
            poromat(i)%vmin = min(vminA,vminB)
            poromat(i)%vmax = max(vmaxA,vmaxB)
        enddo
        !mat%types = mattypes
        close(lu)
        deallocate(iregion)
    end subroutine readPorousMaterialProperties

    !Subroutine to read the material index file. Preferable this file should be in materials mod.
    subroutine readMaterialFile(par, matInd, nelem, lu, filename, errmsg)
        !input
        type(parameterVar) :: par
        type(error_message) :: errmsg
        integer :: lu
        integer :: nelem
        character(len=*) :: filename
        !in/out
        type(materialIndizes) :: matInd
        !loacl
        integer :: i
        integer :: dummy
        integer :: ios
        character(len=16) :: myname = "readMaterialFile"

        call addTrace(errmsg, myname)

        if (par%log) write(*,'(a80)') "|                              read material file                              |"
        allocate(matInd%type(nelem))
        allocate(matInd%pml(nelem))
        !filename=trim('cubit/mat')
        open(unit=lu,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif
        do i=1,nelem
           ! materials
           read(lu,*) dummy,matInd%type(i),matInd%pml(i)   !read Material index -> which material is used in that element, and the index identifiing pmls
        end do
        close(lu)
    end subroutine readMaterialFile

    subroutine calcAPAM(A,AP,AM,vmax,vmin)
        !input
        real(kind=custom_real), dimension(:,:), intent(in) :: A
        !output
        real(kind=custom_real), dimension(:,:), intent(out) :: AP, AM
        real(kind=custom_real), intent(out) :: vmax, vmin
        !local variables
        type(error_message) :: errmsg
        character(len=8) :: myname = 'calcAPAM'
        character(len=100) :: errstr
        real(kind=custom_real), dimension(:,:), allocatable :: Awork
        real(kind=custom_real), dimension(:), allocatable :: WR, WI    !contain real and imaginary part of eigenvalues, respectively
        real(kind=custom_real), dimension(:,:), allocatable :: VR      !matrix containing the eigenvectors as columns
        real(kind=custom_real), dimension(:,:), allocatable :: LambdaP, LambdaM  !matrix containing the eigenvalues on the diagonal elements
        real(kind=custom_real), dimension(:,:), allocatable :: VL,invVR
        real(kind=custom_real), dimension(:), allocatable :: work
        integer :: N, lwork, info, i
        integer, dimension(:), allocatable :: ipiv
        integer, dimension(2) :: shapeA

        shapeA=shape(A)
        N=shapeA(1)
        allocate(Awork(N,N))
        allocate(ipiv(N))
        allocate(WR(N))
        allocate(WI(N))
        allocate(VR(N,N))
        allocate(LambdaP(N,N))
        allocate(LambdaM(N,N))
        allocate(invVR(N,N))
        Awork = A

        !calculate eigenvectors and eigenvalues
        if (custom_real==size_double) then
            !do workspace query
            allocate(work(1)); lwork = -1
            !SYNTAX: dgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )   (LAPACK: http://www.netlib.org/lapack/double/dgeev.f)
            call dgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Workspace query failed: LAPACK routine DGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            !compute eigenvalues and eigenvectors
            call dgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Computation of eigenvalues and -vectors failed: LAPACK routine DGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            deallocate(work)
        elseif (custom_real==size_real) then
            !do workspace query
            allocate(work(1)); lwork = -1
            !SYNTAX: sgeev( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )   (LAPACK: http://www.netlib.org/lapack/double/dgeev.f)
            call sgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Workspace query failed: LAPACK routine SGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            !compute eigenvalues and eigenvectors
            call sgeev('N','V', N, Awork, N, WR, WI, VL, N, VR, N, work, lwork, info)
            if (info/=0) then
                write(errstr,*) 'Computation of eigenvalues and -vectors failed: LAPACK routine SGEEV returned INFO = ',info
                call add(errmsg,2,trim(errstr),myname)
                stop
            endif
            deallocate(work)
        else
            call add(errmsg,2,'CUSTOM_REAL is neither SIZE_REAL nor SIZE_DOUBLE',myname)
            stop
        endif

        LambdaP = 0.
        LambdaM = 0.
        do i=1,N
            if (WR(i) > 0) then
                LambdaP(i,i)=WR(i)
            else
                LambdaM(i,i)=WR(i)
            endif
        enddo

        invVR = VR
        if (custom_real==size_double) then
            call dgetrf(N, N, invVR, N, ipiv, info)
            !do workspace query
            allocate(work(1)); lwork = -1
            call dgetri(N, invVR, N, ipiv, work, lwork, info)
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            call dgetri(N, invVR, N, ipiv, work, lwork, info)
            deallocate(work)
        elseif (custom_real==size_real) then
            call sgetrf(N, N, invVR, N, ipiv, info)
            !do workspace query
            allocate(work(1)); lwork = -1
            call sgetri(N, invVR, N, ipiv, work, lwork, info)
            lwork = int(work(1))
            deallocate(work); allocate(work(lwork))
            call sgetri(N, invVR, N, ipiv, work, lwork, info)
            deallocate(work)
        else
            call add(errmsg,2,'CUSTOM_REAL is neither SIZE_REAL nor SIZE_DOUBLE',myname)
            stop
        endif

        AP = 2.*matmul(VR,matmul(LambdaP,invVR))
        AM = 2.*matmul(VR,matmul(LambdaM,invVR))
        vmax = maxval(WR)
        vmin = 300000.
        do i=1,size(WR)
            if (WR(i) > 0) vmin = min(vmin,WR(i))
        enddo

        if (allocated(Awork)) deallocate(Awork)
        if (allocated(ipiv)) deallocate(ipiv)
        if (allocated(WR)) deallocate(WR)
        if (allocated(WI)) deallocate(WI)
        if (allocated(VR)) deallocate(VR)
        if (allocated(LambdaP)) deallocate(LambdaP)
        if (allocated(LambdaM)) deallocate(LambdaM)
        if (allocated(invVR)) deallocate(invVR)
    end subroutine calcAPAM

    subroutine prepareAnelasticCalc(par, mat, errmsg)
        !input
        type(parameterVar)  :: par
        type(materialVar), dimension(:)   :: mat
        type(error_message) :: errmsg
        !local
        integer :: i, mbs
        character(len=20) :: myname = "prepareAnelasticCalc"

        call addTrace(errmsg, myname)

        do i = 1, par%matn
            if (par%log) then
                write(*,'(a80)') "|------------------------------------------------------------------------------|"
                write(*,'(a40, I3.3, a37)')  "|                 Attenuation Material: ",i,"                                    |"
                write(*,'(a40, f12.6, a28)') "|                                   vp: ", mat(i)%vp, "                           |"
                write(*,'(a40, f12.6, a28)') "|                                   vs: ", mat(i)%vs, "                           |"
                call calcAnelasticCoefficients(par, mat(i), i, errmsg)
                mat(i)%vp = sqrt((mat(i)%lambda+2*mat(i)%mu)/mat(i)%rho) !vp_i
                mat(i)%vs = sqrt(mat(i)%mu/mat(i)%rho) !vs
                write(*,'(a40, f12.6, a28)') "|                      vp (attenuated): ", mat(i)%vp, "                           |"
                write(*,'(a40, f12.6, a28)') "|                      vs (attenuated): ", mat(i)%vs, "                           |"
                write(*,'(a80)')  "|                               omegaL:                                        |"
                do mbs = 1, nMB
                    write(*,'(a40, f12.6, a28)')  "|                                       ", mat(i)%omegaL(mbs), "                           |"
                enddo
            endif
        end do
        if (par%log) then
            write(*,'(a80)') "|------------------------------------------------------------------------------|"
            write(*,'(a80)') "|                        end anelastic pre calculation                         |"
            write(*,'(a80)') "|------------------------------------------------------------------------------|"
        end if
    end subroutine prepareAnelasticCalc

    subroutine calcAnelasticCoefficients(par, mat, matnum, errmsg)
        ! For information and equations see Kaeser et al. (2006): An arbitrary high-order Discontinous Galerkin
        ! method for elastic waves on unstructured  meshes - III. Viscoelastic attenuation
        !input
        type(parameterVar) :: par
        type(materialVar) :: mat
        type(error_message) :: errmsg
        integer :: matnum
        !local
        integer :: i,j,m,c
        real(kind=CUSTOM_REAL) :: theta1, theta2, R
        real(kind=CUSTOM_REAL) :: Mup, Mus
        real(kind=CUSTOM_REAL) :: omegaMin
        real(kind=CUSTOM_REAL), dimension(nMB) :: yp,ys
        real(Kind=CUSTOM_REAL), dimension(:), allocatable :: omegaK, invQp, invQs
        real(Kind=CUSTOM_REAL), dimension(:,:), allocatable :: Mp,Ms
        character(len=25) :: myname = "calcAnelasticCoefficients"

        !Variables for testoutput
        real(kind=custom_real), dimension(100) :: qtest, wtest
        real(kind=custom_real) :: temp1, temp2, vptest, vstest
        complex(kind=CUSTOM_REAL) :: itemp1, itemp2,iltemp,imtemp
        integer :: dim=100 ! number of datapoints to plot test functions for Qp and vp
        character(len=255) :: filename

        !Add Module to error message trace
        call addTrace(errmsg, myname)
        !nMB=number of Maxwell bodies, m=number of values for wk (frequencies) to obtain constant quality factor (see Kaeser et al.)
        m = 2*nMB - 1

        ! maximum and minimum frequency of omegaK use f0 and fr from Parfile
        omegaMin = par%f_max_att/par%fac_att
        if (par%log) then
            write(*,'(a40, f12.6, a28)') "|                             omegaMin: ", omegaMin, "                           |"
            write(*,'(a40, f12.6, a28)') "|                             omegaMax: ", par%f_max_att, "                           |"
        end if
        ! prepare frequencies for determination and inverse values for Qp and Qs
        allocate(omegaK(m), invQp(m), invQs(m))

        ! setup equidistant (log) frequencies in the requestet frequency band
        if (nMB > 1) then
            do i = 1, m
                omegaK(i) = exp(log(omegaMin) + (i-1.0)/(2.0*(nMB-1)) * log(par%fac_att) )
            end do
        else
            ! if there is only one Maxwell body use maximum frequency
            omegaK(:)= par%f0_att
        end if
        c = 1
        do i = 1, m
            if(mod(i, 2) == 1) then
                ! get the relaxation frequencies of the nMB mechanisms (see text between eq. 5 and 6 of Kaeser et al.)
                mat%omegaL(c) = omegaK(i)
                c = c + 1
            end if
        end do

        !set up matrix of the system where we wat to optain the anelastic constants for the wavespeeds
        allocate(Mp(m,nMB))
        allocate(Ms(m,nMB))
        do i=1,m
           do j=1,nMB
              Mp(i,j) =  ( mat%omegaL(j)*omegaK(i)+mat%omegaL(j)**2/mat%Qp ) / (mat%omegaL(j)**2+omegaK(i)**2)
              Ms(i,j) =  ( mat%omegaL(j)*omegaK(i)+mat%omegaL(j)**2/mat%Qs ) / (mat%omegaL(j)**2+omegaK(i)**2)
           end do
        end do

        ! define inverse quality factors at every used frequency wk
        invQp(:) = 1.0/mat%Qp
        invQs(:) = 1.0/mat%Qs
        ! solve the system invQp=Mp for yp (same for ys)
        call solveLinearSystemQR(Mp, invQp, yp, m, nMB, errmsg, CUSTOM_REAL)
        call solveLinearSystemQR(Ms, invQs, ys, m, nMB, errmsg, CUSTOM_REAL)

        ! get the unrelaxed Lame parameters mu and lambda
        !P - Waves
        theta1 = 1.0
        theta2 = 0.0
        do j = 1, nMB
           theta1 = theta1 - yp(j)/(1+(par%f0_att/mat%omegaL(j))**2)
           theta2 = theta2 + yp(j)*(par%f0_att/mat%omegaL(j))/(1+(par%f0_att/mat%omegaL(j))**2)
        end do
        R   = sqrt(theta1**2 + theta2**2)
        Mup = mat%rho * mat%vp**2 * (R + theta1)/(2*R**2)

        !S - Waves
        theta1 = 1.0
        theta2 = 0.0
        do j = 1, nMB
           theta1 = theta1 - ys(j)/(1+(par%f0_att/mat%omegaL(j))**2)
           theta2 = theta2 + ys(j)*(par%f0_att/mat%omegaL(j))/(1+(par%f0_att/mat%omegaL(j))**2)
        end do
        R   = sqrt(theta1**2+theta2**2)
        Mus = mat%rho * mat%vs**2 * (R + theta1)/(2*R**2)

        ! Lame parameters:
        mat%mu = Mus
        mat%lambda = Mup - 2*mat%mu

        ! get the anelastic coefficients for the Lame parameters
        do j=1,nMB
           mat%ylambda(j) = (1+ (2 * mat%mu) / mat%lambda) * yp(j) - (2 * mat%mu)/mat%lambda * ys(j)
           mat%ymu(j)     = ys(j)
        end do


        !The remaining part provides two files to print Qp and vp as a function of frequency

        do i = 1, dim
           wtest(i) = (i/2)**2 ! frequencies where Qp and vp will be plotted
        end do

        open(unit=27,file="out/wktest",status='unknown')
        do i = 1, m
            write(27,*) i,omegaK(i)
        end do
        close(27)

        write(filename,"('out/disp_Q_mat',i3.3)") matnum
        open(unit=27,file=trim(filename),status='unknown') !Q(omega)
        do i=1,dim
           temp1=0
           temp2=0
           do j=1,nMB
              temp1=temp1 + ( mat%omegaL(j)**2 * yp(j) ) / ( mat%omegaL(j)**2 + wtest(i)**2 )
              temp2=temp2 + ( mat%omegaL(j)*wtest(i) * yp(j) ) / ( mat%omegaL(j)**2 + wtest(i)**2 )
           end do
           qtest(i) = (1.0-temp1)/temp2
           write(27,*) wtest(i),qtest(i)
        end do
        close(27)

        write(filename,"('out/disp_v_mat',i3.3)") matnum
        open(unit=27,file=trim(filename),status='unknown') ! vp(omega)
        do i=1,dim
           itemp1=0
           itemp2=0
           do j=1,nMB
              itemp1=itemp1 + mat%omegaL(j)*mat%ylambda(j)/(mat%omegaL(j)+cmplx(0.0,1.0)*wtest(i))
              itemp2=itemp2 + mat%omegaL(j)*mat%ymu(j)/(mat%omegaL(j)+cmplx(0.0,1.0)*wtest(i))
           end do
           iltemp = mat%lambda*(1-itemp1)
           imtemp = mat%mu*(1-itemp2)

           vptest = real(sqrt((iltemp+2*imtemp)/mat%rho))
           vstest = real(sqrt(imtemp/mat%rho))
           write(27,*) wtest(i),vptest, vstest
        end do
        close(27)

        deallocate(Mp,Ms)
        deallocate(omegaK, invQp, invQs)

    end subroutine calcAnelasticCoefficients
end module materialsMod
