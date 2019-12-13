!-----------------------------------------------------------------------
!   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
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
!Module to compute the slip-interface influence for the right-hand-side of the dg equation
module rhsSlipInterface
    use parameterMod
    use constantsMod
    use slipInterfaceMod
    use meshMod
    use errorMessage
    use waveMod
    use mpiMod

    implicit none

    contains

    subroutine initialSlipInterfaceActivity(mesh, lsi, lsi_spec, lsipar, etolsi, q, qi,errmsg)
        !Subroutine to set the initil activity of the slip Interface and allocate some remaining fields.
        !Initialisation of the interface activity has to be done befor the timeloop!
        !input
        type(meshVar) :: mesh                                               !Information on the mesh
        type(lsiVar), dimension(:) :: lsi                                   !Data on all interfaces
        type(interface_spec), dimension(:) :: lsi_spec                      !General specifications for the interfaces
        type(elementToLSI), dimension(:) :: etolsi                          !Mapping: element to slip interface
        type(lsi_parameter) :: lsipar                                       !Type containing switches for certain parts of the interfaces
        type(error_message) :: errmsg                                       !Data on the errormessage
        real(kind=custom_real), dimension(:,:), intent(in) :: q             !Solutionvector
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: qi        !solutionvector for mpi connections
        !local
        character(len=17) :: myname = 'rhsLSI - activity'
        integer :: l, surface, i                                            !Number of the surface and two general counter
        integer :: element, neighbour                                       !Number of the elemnt and its neighbouring element
        integer :: mpi_e, mpi_n                                             !Indices for the mpi-connection
        real(kind=custom_real) :: zin_vp, zout_vp, zin_vs, zout_vs          !local storage for impedances
        real(kind=custom_real), dimension(NGLL,5) :: q_neighbour, q_surface !Solution for one surface in the element and its neighbour

        call addTrace(errmsg,myname)

        !Set Initial activity of the interfaces
        l = 1  !lsicounter
        do element = 1, mesh%nelem
            if (etolsi(element)%isLSI) then
                do surface = 1, nsurface
                    if (etolsi(element)%lsisurface(surface)) then
                        neighbour = mesh%neighbor(surface, element)

                        !Get impedances for the element
                        zin_vp  = mesh%imp_vp(element)
                        zin_vs  = mesh%imp_vs(element)

                        !Set impedances of the neighbouring element to that of the actual element by default for absorbing and reflecting BC
                        zout_vp  = zin_vp
                        zout_vs  = zin_vs

                        !Get the solution vector for the current surface.
                        q_surface   = q(mesh%ibt(:,surface,element),:)

                        !Get information on the neighbour (including boundary conditions and mpi connections)
                        if (neighbour > 0) then
                            q_neighbour = q(mesh%ibn(:,surface,element),:)
                            zout_vp = mesh%imp_vp(neighbour)
                            zout_vs = mesh%imp_vs(neighbour)
                        else if (neighbour == -1) then                                              !surface is an mpi interface
                            mpi_e = mesh%mpi_ibool(surface,element)
                            mpi_n = mesh%mpi_interface(4,surface,element)
                            do i = 1,NGLL
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
                        end if

                        !Calculate the initial activity of the interface jump.
                        select case (trim(lsi_spec(etolsi(element)%lsiprop(surface))%type))
                            case ("elastic")
                                do i = 1, NGLL
                                    if (lsipar%normal) then
                                        lsi(l)%S_n(surface,i) = q_neighbour(i,4) - q_surface(i,4)
                                    end if
                                    if (lsipar%tangential) then
                                        lsi(l)%S_t(surface,i) = q_neighbour(i,5) - q_surface(i,5)
                                    end if
                                end do
                            case default
                                call add(errmsg,2,"No accepted interface-type entered.",myname, "data/interfaces")
                                call print(errmsg); return
                        end select
                    end if
                end do
                lsi(l)%resS_n = 0.
                lsi(l)%rhsS_n = 0.
                lsi(l)%resS_t = 0.
                lsi(l)%rhsS_t = 0.
                l = l+1
            end if
        end do
    end subroutine initialSlipInterfaceActivity

    subroutine activityFluxSlipInterface(element, mesh, lsi, lsipar, lsi_spec, etolsi, rhsQ, rhsQi)
        !This subroutine calcualtes the rhight-hand-side of the evolution equation of the Jump
        ! at the interface for each gll point on the connecting face
        !input
        type(lsiVar) :: lsi                                                                 !Type containing the data on the current interface
        type(lsi_parameter) :: lsipar                                                       !Type containing switches for certain parts of the interface
        type(interface_spec), dimension(:) :: lsi_spec                                      !General specifications of the slip interface
        type(elementToLSI), dimension(:) :: etolsi                                          !Mapping: Element to slip interface
        type(meshVar) :: mesh                                                               !Informations on the mesh
        type(error_message) :: errmsg                                                       !Type containing the error-message
        integer, intent(in) :: element                                                      !Element number
        real(kind=custom_real), dimension(:,:), intent(in) :: rhsQ                          !time derivative of the solution-vector
        real(kind=custom_real), dimension(:,:,:,:), intent(in) :: rhsQi                     !mpi-connection of the time derivative of the solution-vector
        !local
        integer :: surface, i                                                               !Number of surface and general counter
        integer :: neighbour                                                                !neighbouring lsi element
        !integer :: lsi_element                                                             !lsi element
        integer :: mpi_e, mpi_n                                                             !was diese nummern beinhalten muss ich mal noch herausfinden.
        real(kind=custom_real) :: zin_vp, zin_vs, zout_vp, zout_vs                          !local storage for impedances
        real(kind=custom_real) :: nx,nz                                                     !local normals (are newly calculated for each surface)
        real(kind=custom_real), dimension(5,5) :: T,invT                                    !Transformation matrix and its inverse
        real(kind=custom_real), dimension(NGLL,2) :: dV                                     !Velocity subset of the solutionvector
        real(kind=custom_real), dimension(NGLL,5) :: rhsQ_neighbour, rhsQ_surface           !Solution for one surface in the element and its neighbour
        real(kind=custom_real), dimension(NGLL,5) :: rot_rhsQ_neighbour, rot_rhsQ_surface   !Solution for one surface in the element and its neighbour

        do surface = 1, nsurface
            if (etolsi(element)%lsisurface(surface)) then

                !local normals
                nx=mesh%nx(surface*NGLL, element)
                nz=mesh%nz(surface*NGLL, element)

                !Rotation-Matrix and its inverse
                T    = getT(nx,nz,5)
                invT = getInvT(nx,nz,5)

                !Get neighbour
                neighbour = mesh%neighbor(surface, element)

                !Determine Impedances in the Element
                zin_vp  = mesh%imp_vp(element)
                zin_vs  = mesh%imp_vs(element)

                !Set impedances of the neighbouring element to that of the actual element by default for absorbing and reflecting BC
                zout_vp  = zin_vp
                zout_vs  = zin_vs

                rhsQ_surface = rhsQ(mesh%ibt(:,surface,element),:)

                if (neighbour > 0) then
                    rhsQ_neighbour = rhsQ(mesh%ibn(:,surface,element),:)
                    zout_vp = mesh%imp_vp(neighbour)
                    zout_vs = mesh%imp_vs(neighbour)
                else if (neighbour == -1) then                                              !surface is an mpi interface
                    mpi_e = mesh%mpi_ibool(surface,element)
                    mpi_n = mesh%mpi_interface(4,surface,element)
                    do i=1,NGLL
                        rhsQ_neighbour(i,:) = rhsQi(mesh%mpi_ibt(i,surface,element),:,mpi_e,mpi_n)
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
                    rhsQ_neighbour = 0.
                end if

                !Get rotated versions of The solution-Vectors for the element and neighbour
                do i = 1, NGLL
                    rot_rhsQ_surface(i,:) = matmul(invT, rhsQ_surface(i,:))
                    rot_rhsQ_neighbour(i,:) = matmul(invT, rhsQ_neighbour(i,:))
                enddo

                !Get differences in the time derivatives of the solution
                do i = 1, 2
                    !Free surface BC
                    if ((neighbour == 0) .and. (mesh%face(surface, element) == -1)) then
                        dV = 0.
                        !in the matrix "free" all entries regarding velocity are zero.
                        !Since matmul(free, q_surface) is calculating the difference between neighbour and element in that case, dV is zero here.
                    else
                        dV(:,i) = rot_rhsQ_neighbour(:,i+3) - rot_rhsQ_surface(:,i+3)
                    end if
                end do

                !Loop over every gll-point on the connecting face
                do i = 1, NGLL
                    !Jump normal to the interface
                    if (lsipar%normal) then
                        call rhsJump(zin_vp, zout_vp, &
                                     lsi_spec(etolsi(element)%lsiprop(surface))%type, &
                                     lsi_spec(etolsi(element)%lsiprop(surface))%eta_vp, &
                                     lsi%rhsS_n(surface,i), lsi%S_n(surface,i), &
                                     -dV(i,1), rot_rhsQ_surface(i,1), errmsg)
                    endif
                    !Jump Tangential to the interface
                    if (lsipar%tangential) then
                        call rhsJump(zin_vs, zout_vs, &
                                     lsi_spec(etolsi(element)%lsiprop(surface))%type, &
                                     lsi_spec(etolsi(element)%lsiprop(surface))%eta_vs, &
                                     lsi%rhsS_t(surface,i), lsi%S_t(surface,i), &
                                     -dV(i,2), rot_rhsQ_surface(i,3), errmsg)
                    end if
                end do
            end if
        end do
    end subroutine activityFluxSlipInterface

    subroutine rhsJump(zin, zout, type, eta, rhsS, S, dV, Sigma, errmsg)
        !This subourtine contains the equations necessary for the calculation of the rhs of the jump function
        !input
        type(error_message) :: errmsg
        character(len=*), intent(in) :: type
        real(kind=custom_real), intent(in) :: zin, zout, S, Sigma
        real(kind=custom_real), intent(in) :: dV
        real(kind=custom_real), intent(in) :: eta
        !output
        real(kind=custom_real), intent(out) :: rhsS
        !local
        character(len=15) :: myname = 'rhsLSI - Jump'
        real(kind=custom_real) :: nu

        select case (trim(type))
            case ("elastic")
                nu   = (zin + zout)/(eta * zin * zout)
                rhsS = -nu*S + ((zin+zout)/(2*zin*zout))*(Sigma) + 0.5*dV
            case default
                call add(errmsg,2,"No accepted interface-type entered.",myname, "data/interfaces")
                call print(errmsg); return
        end select
    end subroutine rhsJump
end module rhsSlipInterface

