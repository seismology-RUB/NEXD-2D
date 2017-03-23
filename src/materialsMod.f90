!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of NEXD 2D.
!
!   NEXD 2D is free software: you can redistribute it and/or modify it 
!   under the terms of the GNU General Public License as published by the 
!   Free Software Foundation, either version 3 of the License, or (at your 
!   option) any later version.
!
!   NEXD 2D is distributed in the hope that it will be useful, but WITHOUT
!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!   FITNESS FOR A PARTICULAR PURPOSE. 
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License v3.0
!   along with NEXD 2D. If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------
module materialsMod
    
    use parameterMod
    use constantsMod
    use errorMessage
    use linearSystemMod
    
    implicit none
    
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
    end type materialVar
    
    type materialIndizes
        integer, dimension(:), allocatable :: type           !Index of the material-type that is used in each element 
        integer, dimension(:), allocatable :: pml            !Index on wether an element belongs to a pml or not
    end type materialIndizes
        
    type attenuationVar        
        real(kind=custom_real) :: vp!_u                      !P-wave velocity (unrelaxed) 
        real(kind=custom_real) :: vs!_u                      !S-wave velocity (unrelaxed) 
        real(kind=custom_real) :: lambda!_u                  !1st Lamé parameter (unrelaxed)
        real(kind=custom_real) :: mu!_u                      !2nd Lamé parameter (unrelaxed)
        real(kind=CUSTOM_REAL), dimension(nMB) :: ylambda   !Anelastic coefficient regarding the 1st Lamé parameter
        real(kind=CUSTOM_REAL), dimension(nMB) :: ymu       !Anelastic coefficient regarding the 2st Lamé parameter
        real(kind=CUSTOM_REAL), dimension(nMB) :: omegaL    !Relaxation frequencies of the Maxwellbodys
    end type attenuationVar
    
    contains
        
    subroutine setupMaterials(par, mat, matInd, att, nelem, iregion, errmsg)
        !This subroutine initializes the material variables for the different cases (elastic, attenuation and in the future poroelastic)
        !input
        type(parameterVar) :: par
        type(error_message) :: errmsg
        integer :: nelem
        integer, dimension(:), allocatable :: iregion
        !in/out
        type(materialIndizes) :: matInd
        type(materialVar), dimension(:), allocatable :: mat
        type(attenuationVar) , dimension(:), allocatable :: att     
        !local
        character(len=14) :: myname = "setupMaterials"
        
        call addTrace(errmsg, myname)
        
        ! materials values
        call readMaterialProperties(par, mat, att, 19, iregion, trim('cubit/matprop'), errmsg)
        
        ! materials array
        call readMaterialFile(par, matInd, nelem, 19, trim('cubit/mat'), errmsg)
        
        !Attenuation
        if (par%attenuation) then
            call prepareAnelasticCalc(par, mat, att, errmsg)
        end if
    end subroutine setupMaterials
    
    subroutine readMaterialProperties(par, mat, att, lu, iregion, filename, errmsg)
        !Subroutine to read properties of the materials from the cubit/trelis-generated file
        !input
        type(parameterVar) :: par
        type (error_message) :: errmsg
        integer, intent(in) :: lu       
        character(len=*), intent(in) :: filename
        !in/output
        type(materialVar), dimension(:), allocatable :: mat
        type(attenuationVar), dimension(:), allocatable :: att
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
        allocate(att(par%matn))
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
    end subroutine readMaterialProperties
      
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
        
        if (par%log) write(*,*) "read material file"
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
    
    
    subroutine prepareAnelasticCalc(par, mat, att, errmsg)
        !input
        type(parameterVar)  :: par
        type(materialVar), dimension(:)   :: mat
        type(error_message) :: errmsg
        !in/output
        type(attenuationVar), dimension(:):: att
        !local
        integer :: i
        character(len=20) :: myname = "prepareAnelasticCalc"
        
        call addTrace(errmsg, myname)
        
        do i = 1, par%matn
            if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log) write(*,*) "Attenuation Material ",i,": "
            if (par%log) write(*,*) "vp    : ", mat(i)%vp, "vs    : ", mat(i)%vs
            call calcAnelasticCoefficients(par, mat(i), att(i), errmsg)
            att(i)%vp = sqrt((att(i)%lambda+2*att(i)%mu)/mat(i)%rho) !vp_i
            att(i)%vs = sqrt(att(i)%mu/mat(i)%rho) !vs_u
            if (par%log) write(*,*) "vp (attenuated): ", att(i)%vp, "vs (attenuated): ", att(i)%vs
            if (par%log) write(*,*) "omegaL ", att(i)%omegaL(:)
        end do
        if (par%log) write(*,*) "end anelastic pre calculation"
    end subroutine prepareAnelasticCalc
    
    subroutine calcAnelasticCoefficients(par, mat, att, errmsg)
        ! For information and equations see Kaeser et al. (2006): An arbitrary high-order Discontinous Galerkin
        ! method for elastic waves on unstructured  meshes - III. Viscoelastic attenuation
        !input
        type(parameterVar) :: par      
        type(materialVar) :: mat
        type(error_message) :: errmsg
        !output
        type(attenuationVar) :: att
        !local
        integer :: i,j,m,c
        real(kind=CUSTOM_REAL) :: theta1, theta2, R
        real(kind=CUSTOM_REAL) :: Mup, Mus
        real(kind=CUSTOM_REAL) :: omegaMax, omegaMin        
        real(kind=CUSTOM_REAL), dimension(nMB) :: yp,ys
        real(Kind=CUSTOM_REAL), dimension(:), allocatable :: omegaK, invQp, invQs
        real(Kind=CUSTOM_REAL), dimension(:,:), allocatable :: Mp,Ms
        character(len=25) :: myname = "calcAnelasticCoefficients"
              
        !Variables for testoutput
        real(kind=custom_real), dimension(1000) :: qtest, wtest
        real(kind=custom_real) :: temp1, temp2, vptest
        complex :: itemp1, itemp2,iltemp,imtemp  
        integer :: dim=1000 ! number of datapoints to plot test functions for Qp and vp
        
        !Add Module to error message trace
        call addTrace(errmsg, myname)
        !nMB=number of Maxwell bodies, m=number of values for wk (frequencies) to obtain constant quality factor (see Kaeser et al.)
        m = 2*nMB - 1
        
        ! maximum and minimum frequency of omegaK use f0 and fr from Parfile       
        omegaMax = 2.5*pi*par%f0_att
        omegaMin = par%f0_att/sqrt(par%fr_att)
        write(*,*) "omegaMin: ",omegaMin
        write(*,*) "omegaMax: ", omegaMax
        ! prepare frequencies for determination and inverse values for Qp and Qs
        allocate(omegaK(m), invQp(m), invQs(m))
        
        ! setup equidistant (log) frequencies in the requestet frequency band
        if (nMB > 1) then
            do i = 1, m
                omegaK(i) = exp(log(omegaMin) + (i-1.0)/(2.0*(nMB-1)) * log(par%fr_att) )
            end do
        else
            ! if there is only one Maxwell body use maximum frequency
            omegaK(:)=omegaMax
        end if        
        c = 1 
        do i = 1, m
            if(mod(i, 2) == 1) then
                ! get the relaxation frequencies of the nMB mechanisms (see text between eq. 5 and 6 of Kaeser et al.)
                att%omegaL(c) = omegaK(i) 
                c = c + 1
            end if
        end do

        !set up matrix of the system where we wat to optain the anelastic constants for the wavespeeds
        allocate(Mp(m,nMB))
        allocate(Ms(m,nMB))
        do i=1,m
           do j=1,nMB
              Mp(i,j) =  ( att%omegaL(j)*omegaK(i)+att%omegaL(j)**2/mat%Qp ) / (att%omegaL(j)**2+omegaK(i)**2)
              Ms(i,j) =  ( att%omegaL(j)*omegaK(i)+att%omegaL(j)**2/mat%Qs ) / (att%omegaL(j)**2+omegaK(i)**2)
           end do
        end do

        ! define inverse quality factors at every used frequency wk
        invQp(:) = 1.0/mat%Qp
        invQs(:) = 1.0/mat%Qs
        ! solve the system invQp=Mp for yp (same for ys)
        call solveLinearSystemQR(Mp, invQp, yp, m, nMB, errmsg)
        call solveLinearSystemQR(Ms, invQs, ys, m, nMB, errmsg)

        ! get the unrelaxed Lame parameters mu_u and lambda_u
        !P - Waves
        theta1 = 1.0
        theta2 = 0.0
        do j = 1, nMB
           theta1 = theta1 - yp(j)/(1+(par%f0_att/att%omegaL(j))**2)
           theta2 = theta2 + yp(j)*(par%f0_att/att%omegaL(j))/(1+(par%f0_att/att%omegaL(j))**2)
        end do
        R   = sqrt(theta1**2 + theta2**2)
        Mup = mat%rho * mat%vp**2 * (R + theta1)/(2*R**2)

        !S - Waves
        theta1 = 1.0
        theta2 = 0.0
        do j = 1, nMB
           theta1 = theta1 - ys(j)/(1+(par%f0_att/att%omegaL(j))**2)
           theta2 = theta2 + ys(j)*(par%f0_att/att%omegaL(j))/(1+(par%f0_att/att%omegaL(j))**2)
        end do
        R   = sqrt(theta1**2+theta2**2)
        Mus = mat%rho * mat%vs**2 * (R + theta1)/(2*R**2)
        
        ! Lame parameters:
        att%mu = Mus
        att%lambda = Mup - 2*att%mu
        
        ! get the anelastic coefficients for the Lame parameters
        do j=1,nMB
           att%ylambda(j) = (1+ (2 * att%mu) / att%lambda) * yp(j) - (2 * att%mu)/att%lambda * ys(j)
           att%ymu(j)     = ys(j)
        end do
        
        
        !The remaining part provides two files to print Qp and vp as a function of frequency
        
        do i = 1, dim
           wtest(i) = (i/2)**2 ! frequencies where Qp and vp will be plotted  
        end do
              
        open(unit=27,file="wktest",status='unknown')
        do i = 1, m
            write(27,*) i,omegaK(i)    
        end do
        close(27)
              
        open(unit=27,file="qtest",status='unknown') !Q(omega)
        do i=1,dim
           temp1=0
           temp2=0
           do j=1,nMB
              temp1=temp1 + ( att%omegaL(j)**2 * yp(j) ) / ( att%omegaL(j)**2 + wtest(i)**2 )
              temp2=temp2 + ( att%omegaL(j)*wtest(i) * yp(j) ) / ( att%omegaL(j)**2 + wtest(i)**2 )
           end do
           qtest(i) = (1.0-temp1)/temp2
           write(27,*) wtest(i),qtest(i)
        end do
        close(27)
    
        open(unit=27,file="vptest",status='unknown') ! vp(omega) 
        do i=1,dim
           itemp1=0
           itemp2=0
           do j=1,nMB
              itemp1=itemp1 + att%omegaL(j)*att%ylambda(j)/(att%omegaL(j)+cmplx(0.0,1.0)*wtest(i))
              itemp2=itemp2 + att%omegaL(j)*att%ymu(j)/(att%omegaL(j)+cmplx(0.0,1.0)*wtest(i))
           end do
           iltemp = att%lambda*(1-itemp1)
           imtemp = att%mu*(1-itemp2)
    
           vptest = real(sqrt((iltemp+2*imtemp)/mat%rho))
           write(27,*) wtest(i),vptest
        end do
        close(27)
               
        deallocate(Mp,Ms)
        deallocate(omegaK, invQp, invQs)

    end subroutine calcAnelasticCoefficients
end module materialsMod
