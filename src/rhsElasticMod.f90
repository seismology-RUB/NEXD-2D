!-----------------------------------------------------------------------
!   Copyright 2015-2020 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2020 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2020 Marc S. Boxberg (RWTH Aachen University, GER)
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
module rhsElasticMod

    use constantsMod
    use mpiMod
    use slipInterfaceMod
    use riemannFluxLsiMod
    use waveMod
    use meshMod

    implicit none

    contains

    subroutine rhsElastic(element, iv, ilsi, mesh, elasticfluxvar, lsi, lsipar, etolsi, free, q, qi, rhsQ, f, g, fprime, gprime, pmlcheck, kx, kz)
        !input
        type(meshVar) :: mesh
        type(lsiVar), dimension(:) :: lsi                                   !Container of variables for the interface data
        type(elementToLSI), dimension(:) :: etolsi                          !Mapping: element to lsi. Contains informations on whether an element is an lsi or not
        type(lsi_parameter) :: lsipar                                       !General parameters of the interface
        integer, intent(in) :: ilsi                                         !Current number of the interface to be calculated
        integer, intent(in) :: element                                      !Current element number
        integer, dimension(:), intent(in) :: iv                             !Vector that determines which portion of the whole soltution vector is to be accesed
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: kx, kz          !Some factors for pml claculation
        real(kind=custom_real), dimension(:,:), intent(in) :: free          !Matrix for the selecting boundary conditions
        real(kind=custom_real), dimension(:,:), intent(in) :: q             !Solutionvector
        real(kind=custom_real), dimension(:,:), intent(inout) :: f, g       !Internal flux
        real(kind=custom_real), dimension(:,:), intent(in) :: fprime, gprime!Modified fluxes for pml calculation
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: qi        !solutionvector for mpi connections
        logical, dimension(:), intent(in) :: pmlcheck                       !Is an element part of the pml?
        !output
        real(kind=CUSTOM_REAL), dimension(:,:), intent(out) :: rhsQ         !Right-hand-side of the dg equation
        !local
        integer :: surface                                                  !Surface number (1-3)
        integer :: neighbour                                                !Number of the neighbouring element
        integer :: i                                                        !counter
        integer :: mpi_e, mpi_n
        real(kind=custom_real) :: zin_vp, zout_vp, zin_vs, zout_vs          !local storage for impedances
        real(kind=custom_real) :: nxe,nze                                   !local normals (are newly calculated for each surface)
        real(kind=CUSTOM_REAL), dimension(Np) :: dFdr,dFds,dGdr,dGds        !Derivative matrices for the elastic fluxes
        real(kind=CUSTOM_REAL), dimension(4) :: elasticfluxvar              !variables to increase calculation of elastic fluxes
        real(kind=custom_real), dimension(5,5) :: T,invT                    !Transformation matrix and its inverse
        real(kind=custom_real), dimension(NGLL,5) :: tmp_flux
        real(kind=custom_real), dimension(NGLL,2) :: gamma                  !Expansioncoefficients for P-SV
        real(kind=custom_real), dimension(NGLL,5) :: delta                  !delta q
        real(kind=custom_real), dimension(NGLL,2) :: dV                     !Velocity subset of the solutionvector
        real(kind=custom_real), dimension(NGLL,3) :: dSigma                 !Stress subset of the solutionvector
        real(kind=custom_real), dimension(NGLL,5) :: q_neighbour, q_surface !Solution for one surface in the element and its neighbour
        real(kind=custom_real), dimension(NP,5) :: q_element                !Solution vector for the whole element
        real(kind=custom_real), dimension(3*NGLL,5) :: flux                 !Flux inside the element (left from the element border in the rotated frame).  Dimension(5,5) (Face Points, (Sigma, V))

        !Initialize flux and gamma-arrays as 0
        flux  = 0.
        gamma = 0.

        do surface = 1, nsurface
            neighbour = mesh%neighbor(surface,element)

            !determine impedances for the element
            zin_vp  = mesh%imp_vp(element)
            zin_vs  = mesh%imp_vs(element)

            !Impedances for the neighbouring element are set to the element impedances for absorbing and reflecting BC.
            !These are overwirtten in case a neighbouring element exists.
            zout_vp = zin_vp
            zout_vs = zin_vs

            !Solution in the element and its neighbour (for one surface at a time)
            q_surface   = q(mesh%ibt(:,surface,element),:)

            if (neighbour > 0) then
                q_neighbour = q(mesh%ibn(:,surface,element),:)
                zout_vp = mesh%imp_vp(neighbour)
                zout_vs = mesh%imp_vs(neighbour)
            else if (neighbour == -1) then                                              !surface is an mpi interface
                mpi_e = mesh%mpi_ibool(surface,element)
                mpi_n = mesh%mpi_interface(4,surface,element)
                do i=1,NGLL
                    q_neighbour(i,:) = qi(mesh%mpi_ibt(i,surface,element),:,mpi_e,mpi_n)
                end do
                !Loop over all faces that sit on an mpi-interface for that rank/processor
                do i = 1, size(mesh%mpi_imp_vp(:,1))
                    !Check if the rank and the element in that rank are in the mpi-impedance array (they should)
                    if ((abs(mesh%mpi_imp_vp(i,2) - mesh%mpi_interface(1,surface,element)) < epsilon(mesh%mpi_imp_vp(i,2))) .and. &
                        (abs(mesh%mpi_imp_vp(i,3) - mesh%mpi_interface(2,surface,element)) < epsilon(mesh%mpi_imp_vp(i,3)))) then
                        !associate corresponding impedance value to zout
                        zout_vp = mesh%mpi_imp_vp(i,1)
                        zout_vs = mesh%mpi_imp_vs(i,1)
                    end if
                end do
            else if ((neighbour == 0) .and. mesh%face(surface, element) == -2) then           !Absorbing boundary conditions
                q_neighbour = 0.
            endif

            !local normals
            nxe=mesh%nx(surface*NGLL, element)
            nze=mesh%nz(surface*NGLL, element)

            !Rotation-Matrix and its inverse
            T    = getT(nxe,nze,5)
            invT = getInvT(nxe,nze,5)

            !Get differences of the solutionvectors
            do i = 1, NGLL
                !Relecting boundary (free surface)
                if ((neighbour == 0) .and. (mesh%face(surface, element) == -1)) then
                    delta(i,:) = matmul(free, matmul(invT, q_surface(i,:)))
                else
                    delta(i,:) = matmul(invT, q_neighbour(i,:) - q_surface(i,:))
                endif
            enddo

            !Calculate differences in velocities and stresses
            dV(:,1)      = delta(:,4)
            dV(:,2)      = delta(:,5)
            dSigma(:, 1) = delta(:,1)
            dSigma(:, 2) = delta(:,2)
            dSigma(:, 3) = delta(:,3)

            !calculate expansion coefficients for riemann fluxes
            call riemannExpansionCoefficients(zin_vp, zout_vp, zin_vs, zout_vs, dV, dSigma, gamma)
            !if lsi is to be calculated and the surface is the slip interface add that influence to the flux

            if (lsipar%lsi) then !n > 0
                if (etolsi(element)%isLSI) then
                    if (etolsi(element)%lsisurface(surface)) then
                        call riemannExpansionCoefficientsAddLSI(lsipar, zin_vp, zout_vp, zin_vs, zout_vs,&
                                                                lsi(ilsi)%S_n(surface,:), lsi(ilsi)%S_t(surface,:), gamma)
                    end if
                end if
            end if

            !Calculate Riemann Flux
            tmp_flux(:,1) = - mesh%vp(element) * zin_vp * gamma(:,2)
            tmp_flux(:,2) = - mesh%vp(element) * (zin_vp - 2*mesh%vs(element)/mesh%vp(element) * zin_vs) * gamma(:,2)
            tmp_flux(:,3) = - mesh%vs(element) * zin_vs * gamma(:,1)
            tmp_flux(:,4) =   mesh%vp(element) * gamma(:,2)
            tmp_flux(:,5) =   mesh%vs(element) * gamma(:,1)

            !Transformation has to be done for each GLL-Point individually as only in that case the matmul-operation
            !gets the correct input (a vector with the values for velocity and stress as its entries)
            do i = 1, NGLL
                !tmp_flux(i,:) = - matmul(T, tmp_flux(i,:))
                tmp_flux(i,:) =  matmul(T, tmp_flux(i,:))
            enddo

            !Writing into the (15x5) array has to be done for each collumn seperatly. Otherwise the wrong values are copied!
            do i = 1, 5
                flux(((surface-1)*NGLL+1):(surface*NGLL),i) = tmp_flux(:,i)
            enddo
        end do

        !All calculations below are done in the non-rotated coordinate system!
        q_element = q(iv,:)
        call elasticFluxes(q_element,elasticfluxvar,f,g)     !Calculate f and g

        !PML part
        if (pmlcheck(element)) then
            do i=1,5
                f(:,i) = (fprime(iv,i) + f(:,i))/kx(iv)
                g(:,i) = (gprime(iv,i) + g(:,i))/kz(iv)
            enddo
        end if

        !Putting the two flux parts together to form the surface integral term of the dg-equation.
        do i = 1, 5
            dFdr = matmul(mesh%Dr,f(:,i))
            dFds = matmul(mesh%Ds,f(:,i))
            dGdr = matmul(mesh%Dr,g(:,i))
            dGds = matmul(mesh%Ds,g(:,i))
            rhsQ(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds) + (mesh%rz(iv)*dGdr + mesh%sz(iv)*dGds)
            rhsQ(iv,i) = rhsQ(iv,i) - matmul(mesh%lift, mesh%fscale(:,element)*flux(:,i))
        end do
    end subroutine rhsElastic
end module rhsElasticMod
