!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2018 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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
module waveMod
    ! the wave equation to solve
    use constantsMod
    use mpiMod

    implicit none

    contains

    subroutine elasticFluxes(q,elasticfluxvar,f,g)
        ! compute elastic fluxes
        !input
        !type(constant_values) :: cv
        real(kind=CUSTOM_REAL), dimension(4) :: elasticfluxvar
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        !output
        real(kind=CUSTOM_REAL), dimension(:,:) :: f,g


        f(:,1)= elasticfluxvar(1)*q(:,4)
        f(:,2)= elasticfluxvar(2)*q(:,4)
        f(:,3)= elasticfluxvar(3)*q(:,5)
        f(:,4)= elasticfluxvar(4)*q(:,1)
        f(:,5)= elasticfluxvar(4)*q(:,3)

        g(:,1)= elasticfluxvar(2)*q(:,5)
        g(:,2)= elasticfluxvar(1)*q(:,5)
        g(:,3)= elasticfluxvar(3)*q(:,4)
        g(:,4)= elasticfluxvar(4)*q(:,3)
        g(:,5)= elasticfluxvar(4)*q(:,2)
    end subroutine elasticFluxes

    subroutine anelasticFluxes(q,wl,f,g)
        ! compute anelastic fluxes
        !input
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        real(kind=CUSTOM_REAL), dimension(nMB) :: wl
        !output
        real(kind=CUSTOM_REAL), dimension(:,:) :: f,g

        real(kind=CUSTOM_REAL), dimension(Np) :: temp1, temp2
        integer :: i

        do i=1,nMB
            temp1=-wl(i)*q(:,4)
            temp2=-wl(i)*q(:,5)
            f(:,(i-1)*3+1) = temp1
            f(:,(i-1)*3+2) = 0.0
            f(:,(i-1)*3+3) = 0.5*temp2

            g(:,(i-1)*3+1) = 0.0
            g(:,(i-1)*3+2) = temp2
            g(:,(i-1)*3+3) = 0.5*temp1
        end do
    end subroutine anelasticFluxes

    function getAPA(vp,vs,rho,lambda,mu)
        !compute A + |A| for exact riemann fluxes
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(5,5) :: getAPA
        getAPA=0.0
        getAPA(1,1) = vp
        getAPA(1,4) = -(lambda+2*mu)
        getAPA(2,1) = (lambda*vp)/(lambda+2*mu)
        getAPA(2,4) = -lambda
        getAPA(3,3) = vs
        getAPA(3,5) = -mu
        getAPA(4,1) = -1/rho
        getAPA(4,4) = vp
        getAPA(5,3) = -1/rho
        getAPA(5,5) = vs
    end function getAPA

    function getAnelasticAPA(vp,vs,rho)
        !compute A + |A| for exact riemann fluxes
        real(kind=CUSTOM_REAL) :: vp,vs,rho
        real(kind=CUSTOM_REAL), dimension(3,5) :: getAnelasticAPA
        getAnelasticAPA=0.0
        getAnelasticAPA(1,1) = 1.0/(vp*rho)
        getAnelasticAPA(1,4) = -1.0
        getAnelasticAPA(3,3) = 1.0/(2*vs*rho)
        getAnelasticAPA(3,5) = -1.0/2.0
    end function getAnelasticAPA

    function getT(nxe,nze)
        !get rotation matrix
        real(kind=CUSTOM_REAL) :: nxe,nze
        real(kind=CUSTOM_REAL), dimension(5,5) :: getT
        getT=0.0
        getT(1,1) = nxe*nxe
        getT(1,2) = nze*nze
        getT(1,3) = -2*nxe*nze
        getT(2,1) = nze*nze
        getT(2,2) = nxe*nxe
        getT(2,3) = 2*nxe*nze
        getT(3,1) = nxe*nze
        getT(3,2) = -nxe*nze
        getT(3,3) = nxe*nxe-nze*nze
        getT(4,4) = nxe
        getT(4,5) = -nze
        getT(5,4) = nze
        getT(5,5) = nxe
    end function getT

    function getInvT(nxe,nze)
        !get inverse rotation matrix
        real(kind=CUSTOM_REAL) :: nxe,nze
        real(kind=CUSTOM_REAL), dimension(5,5) :: getInvT
        getInvT=0.0
        getInvT(1,1) = nxe*nxe
        getInvT(1,2) = nze*nze
        getInvT(1,3) = 2*nxe*nze
        getInvT(2,1) = nze*nze
        getInvT(2,2) = nxe*nxe
        getInvT(2,3) = -2*nxe*nze
        getInvT(3,1) = -nxe*nze
        getInvT(3,2) = nxe*nze
        getInvT(3,3) = nxe*nxe-nze*nze
        getInvT(4,4) = nxe
        getInvT(4,5) = nze
        getInvT(5,4) = -nze
        getInvT(5,5) = nxe
    end function getInvT

   subroutine computeExactRiemannSF(flux,q,qi,neighbor,VT,VTfree,face, mpi_interface,&
                                     mpi_ibool, mpi_ibt,ibt,ibn,num_eq)
        !Riemann fluxes for the Strong form of DG
        integer, intent(in) :: num_eq
        integer, dimension(:) :: neighbor,face
        integer, dimension(:,:) :: ibt,ibn
        real(kind=CUSTOM_REAL), dimension(:,:) :: flux
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
        integer, dimension(:,:) :: mpi_interface
        integer, dimension(:) :: mpi_ibool
        integer, dimension(:,:) :: mpi_ibt
        real(kind=CUSTOM_REAL), dimension(NGLL,5) :: qtemp_n,qtemp_t
        real(kind=CUSTOM_REAL), dimension(NGLL,num_eq) :: temp
        real(kind=CUSTOM_REAL), dimension(:,:,:) :: VT, VTfree
        integer :: is,in,i,j
        integer :: mpi_e, mpi_n
        ! compute exact riemann problem, elementwise

        flux=0.0
        do is=1,3 !surfs
            ! get neighbor element
            in = neighbor(is)
            if (in>0) then ! not boundary of the whole domain
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibn(:,is),:)
                do i=1,NGLL
                    temp(i,:) = matmul(VT(is,:,:),(qtemp_n(i,:)-qtemp_t(i,:)))
                end do
                do j=1,num_eq
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = temp(:,j)
                end do
                !boundaries
                ! "------------------------------------------------------------------"
            else if ((in==0).and.(face(is) == -1)) then ! free
                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    temp(i,:) = matmul(VTfree(is,:,:),qtemp_t(i,:))
                end do
                do j=1,num_eq
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (temp(:,j))
                end do
                ! "------------------------------------------------------------------"
            else if ((in==0).and.face(is) == -2) then !absorb
                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    temp(i,:) = -matmul(VT(is,:,:),qtemp_t(i,:))
                end do
                do j=1,num_eq
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = temp(:,j)
                end do
            else if (in == -1) then
                ! mpi interface
                mpi_e=mpi_ibool(is)
                mpi_n=mpi_interface(4,is)
                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
                end do
                do i=1,NGLL
                    temp(i,:) = matmul(VT(is,:,:),(qtemp_n(i,:)-qtemp_t(i,:)))
                end do
                do j=1,num_eq
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = temp(:,j)
                end do
            end if
        end do
    end subroutine computeExactRiemannSF
end module waveMod
