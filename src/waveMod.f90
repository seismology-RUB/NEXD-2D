!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
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
module waveMod
    ! the wave equation to solve
    use constantsMod
    use matrixMod
    use mpiMod

    implicit none

    contains

    subroutine elasticFluxes(q,lambda,mu,rho,f,g)
        ! compute elastic fluxes
        real(kind=CUSTOM_REAL) :: lambda,mu,rho
        real(kind=CUSTOM_REAL), dimension(:,:) :: f,g
        real(kind=CUSTOM_REAL), dimension(:,:) :: q

        f(:,1)= -(lambda+2.0*mu)*q(:,4)
        f(:,2)= -lambda*q(:,4)
        f(:,3)= -mu*q(:,5)
        f(:,4)= -1.0/rho*q(:,1)
        f(:,5)= -1.0/rho*q(:,3)

        g(:,1)= -lambda*q(:,5)
        g(:,2)= -(lambda+2.0*mu)*q(:,5)
        g(:,3)= -mu*q(:,4)
        g(:,4)= -1.0/rho*q(:,3)
        g(:,5)= -1.0/rho*q(:,2)
    end subroutine elasticFluxes

    subroutine anelasticFluxes(q,wl,f,g)
        ! compute anelastic fluxes
        real(kind=CUSTOM_REAL), dimension(nMB) :: wl
        real(kind=CUSTOM_REAL), dimension(:,:) :: f,g
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        integer :: i

        do i=1,nMB
            f(:,(i-1)*3+1) = -wl(i)*q(:,4)
            f(:,(i-1)*3+2) = 0.0
            f(:,(i-1)*3+3) = -0.5*wl(i)*q(:,5)

            g(:,(i-1)*3+1) = 0.0
            g(:,(i-1)*3+2) = -wl(i)*q(:,5)
            g(:,(i-1)*3+3) = -0.5*wl(i)*q(:,4)
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

    !Only for weak form of DG
    function getAMA(vp,vs,rho,lambda,mu)
        !compute A - |A| for exact riemann fluxes
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(5,5) :: getAMA
        getAMA=0.0
        getAMA(1,1) = -vp
        getAMA(1,4) = -(lambda+2*mu)
        getAMA(2,1) = -(lambda*vp)/(lambda+2*mu)
        getAMA(2,4) = -lambda
        getAMA(3,3) = -vs
        getAMA(3,5) = -mu
        getAMA(4,1) = -1/rho
        getAMA(4,4) = -vp
        getAMA(5,3) = -1/rho
        getAMA(5,5) = -vs
    end function getAMA

    function getAnelasticAPA(vp,vs,rho,lambda,mu)
        !compute A + |A| for exact riemann fluxes
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(3,5) :: getAnelasticAPA
        getAnelasticAPA=0.0
        getAnelasticAPA(1,1) = 1.0/(vp*rho)
        getAnelasticAPA(1,4) = -1.0
        getAnelasticAPA(3,3) = 1.0/(2*vs*rho)
        getAnelasticAPA(3,5) = -1.0/2.0
    end function getAnelasticAPA

    !Only for weak form of DG
    function getAnelasticAMA(vp,vs,rho,lambda,mu)
        !compute A - |A| for exact riemann fluxes
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(3,5) :: getAnelasticAMA
        getAnelasticAMA=0.0
        getAnelasticAMA(1,1) = -1.0/(vp*rho)
        getAnelasticAMA(1,4) = -1.0
        getAnelasticAMA(3,3) = -1.0/(2*vs*rho)
        getAnelasticAMA(3,5) = -1.0/2.0
    end function getAnelasticAMA

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

    subroutine computeExactRiemannWF(flux,q,qi,neighbor,vp,vs,rho,lambda,mu,face, mpi_interface,&
                                     mpi_ibool, mpi_ibt,ibt,ibn,nx,nz,freex,pml,set_pml)
        !Riemann fluxes for the Weak form of DG
        integer, dimension(:) :: neighbor,face
        integer, dimension(:,:) :: ibt,ibn
        logical :: set_pml
        real(kind=CUSTOM_REAL), dimension(:,:) :: freex
        real(kind=CUSTOM_REAL), dimension(3*NGLL,5) :: flux
        real(kind=CUSTOM_REAL), dimension(:) :: nx,nz
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
        integer, dimension(:,:) :: mpi_interface
        integer, dimension(:) :: mpi_ibool
        integer, dimension(:,:) :: mpi_ibt
        real(kind=CUSTOM_REAL), dimension(NGLL,5) :: qtemp_n,qtemp_t
        real(kind=CUSTOM_REAL) :: nxe,nze
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(5,5) :: T,invT
        real(kind=CUSTOM_REAL), dimension(5,5) :: APA,AMA,VT,WT
        integer :: is,in,i,j
        integer :: mpi_e, mpi_n
        logical :: pml
        ! compute exact riemann problem, elementwise

        ! for all equal
        APA = getAPA(vp,vs,rho,lambda,mu)
        AMA = getAMA(vp,vs,rho,lambda,mu)
        flux=0.0

        do is=1,3 !surfs
            ! get neighbor element
            in = neighbor(is)
            if (in>0) then ! not boundary of the whole domain
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibn(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))
                WT= matmul(T,matmul(AMA,invT))
                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                    qtemp_n(i,:) = matmul(WT,qtemp_n(i,:))
                end do
                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (qtemp_t(:,j)+qtemp_n(:,j))
                end do

                ! boundaries
                ! "------------------------------------------------------------------"
            else if ((in==0).and.(face(is) == -1)) then ! free
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibt(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))
                WT= matmul(T,matmul(matmul(AMA,freex),invT))

                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                    qtemp_n(i,:) = matmul(WT,qtemp_n(i,:))
                end do

                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (qtemp_t(:,j)+qtemp_n(:,j))
                end do
            ! "------------------------------------------------------------------"
            else if ((in==0).and.face(is) == -2) then !absorb
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)
                qtemp_t=q(ibt(:,is),:)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))
                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                end do
                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = qtemp_t(:,j)
                end do
            else if (in == -1) then
                !mpi interface
                mpi_e=mpi_ibool(is)
                mpi_n=mpi_interface(4,is)

                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
                end do
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))
                WT= matmul(T,matmul(AMA,invT))
                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                    qtemp_n(i,:) = matmul(WT,qtemp_n(i,:))
                end do
                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (qtemp_t(:,j)+qtemp_n(:,j))
                end do
            end if
        end do
    end subroutine computeExactRiemannWF

    subroutine computeExactRiemannWFAnelastic(flux,q,qi,neighbor,vp,vs,rho,lambda,mu,face, mpi_interface,&
                                              mpi_ibool, mpi_ibt,ibt,ibn,nx,nz,freex,pml,set_pml)
        !Riemann fluxes for the Weak form of DG - anelastic
        integer, dimension(:) :: neighbor,face
        integer, dimension(:,:) :: ibt,ibn
        logical :: set_pml
        real(kind=CUSTOM_REAL), dimension(:,:) :: freex
        real(kind=CUSTOM_REAL), dimension(3*NGLL,3) :: flux
        real(kind=CUSTOM_REAL), dimension(:) :: nx,nz
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
        integer, dimension(:,:) :: mpi_interface
        integer, dimension(:) :: mpi_ibool
        integer, dimension(:,:) :: mpi_ibt
        real(kind=CUSTOM_REAL), dimension(NGLL,5) :: qtemp_n,qtemp_t
        real(kind=CUSTOM_REAL), dimension(NGLL,3) :: atemp_t,atemp_n
        real(kind=CUSTOM_REAL) :: nxe,nze
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(5,5) :: T,invT
        real(kind=CUSTOM_REAL), dimension(3,3) :: aT
        real(kind=CUSTOM_REAL), dimension(3,5) :: VT,WT
        real(kind=CUSTOM_REAL), dimension(3,5) :: aAPA,aAMA
        integer :: is,in,i,j
        integer :: mpi_e, mpi_n
        logical :: pml
        ! compute exact riemann problem, elementwise

        ! for all equal

        aAPA = getAnelasticAPA(vp,vs,rho,lambda,mu)
        aAMA = getAnelasticAMA(vp,vs,rho,lambda,mu)
        flux=0.0
        do is=1,3 !surfs
            ! get neighbor element
            in = neighbor(is)
            if (in>0) then ! not boundary of the whole domain
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibn(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)

                VT= matmul(aT,matmul(aAPA,invT))

                VT= matmul(aT,matmul(aAPA,invT))
                WT= matmul(aT,matmul(aAMA,invT))
                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                    atemp_n(i,:) = matmul(WT,qtemp_n(i,:))
                end do
                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (atemp_t(:,j)+atemp_n(:,j))
                end do

            !boundaries
            !"------------------------------------------------------------------"
            else if ((in==0).and.(face(is) == -1)) then ! free
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibt(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)
                VT= matmul(aT,matmul(aAPA,invT))
                WT= matmul(aT,matmul(matmul(aAMA,freex),invT))

                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                    atemp_n(i,:) = matmul(WT,qtemp_n(i,:))
                end do

                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (atemp_t(:,j)+atemp_n(:,j))
                end do
            ! "------------------------------------------------------------------"
            else if ((in==0).and.face(is) == -2) then !absorb
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)
                qtemp_t=q(ibt(:,is),:)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)
                VT= matmul(aT,matmul(aAPA,invT))
                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                end do

                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = atemp_t(:,j)
                end do
            else if (in == -1) then
                ! mpi interface
                mpi_e=mpi_ibool(is)
                mpi_n=mpi_interface(4,is)

                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
                end do
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)
                VT= matmul(aT,matmul(aAPA,invT))
                WT= matmul(aT,matmul(aAMA,invT))
                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                    atemp_n(i,:) = matmul(WT,qtemp_n(i,:))
                end do
                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (atemp_t(:,j)+atemp_n(:,j))
                end do
            end if
        end do
    end subroutine computeExactRiemannWFAnelastic


    subroutine computeExactRiemannSF(flux,q,qi,neighbor,vp,vs,rho,lambda,mu,face, mpi_interface,&
                                     mpi_ibool, mpi_ibt,ibt,ibn,nx,nz,freex,pml,set_pml)
        !Riemann fluxes for the Strong form of DG
        integer, dimension(:) :: neighbor,face
        integer, dimension(:,:) :: ibt,ibn
        logical :: set_pml
        real(kind=CUSTOM_REAL), dimension(:,:) :: freex
        real(kind=CUSTOM_REAL), dimension(3*NGLL,5) :: flux
        real(kind=CUSTOM_REAL), dimension(:) :: nx,nz
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
        integer, dimension(:,:) :: mpi_interface
        integer, dimension(:) :: mpi_ibool
        integer, dimension(:,:) :: mpi_ibt
        real(kind=CUSTOM_REAL), dimension(NGLL,5) :: qtemp_n,qtemp_t
        real(kind=CUSTOM_REAL) :: nxe,nze
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(5,5) :: T,invT
        real(kind=CUSTOM_REAL), dimension(5,5) :: APA,VT
        integer :: is,in,i,j
        integer :: mpi_e, mpi_n
        logical :: pml
        ! compute exact riemann problem, elementwise

        ! for all equal
        APA = getAPA(vp,vs,rho,lambda,mu)

        flux=0.0

        do is=1,3 !surfs
            ! get neighbor element
            in = neighbor(is)
            if (in>0) then ! not boundary of the whole domain
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibn(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))

                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,(qtemp_n(i,:)-qtemp_t(i,:)))
                end do

                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = qtemp_t(:,j)
                end do
                !boundaries
                ! "------------------------------------------------------------------"
            else if ((in==0).and.(face(is) == -1)) then ! free
                qtemp_t=q(ibt(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(matmul(APA,freex),invT))

                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                end do

                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (qtemp_t(:,j))
                end do
                ! "------------------------------------------------------------------"
            else if ((in==0).and.face(is) == -2) then !absorb
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)
                qtemp_t=q(ibt(:,is),:)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))
                do i=1,NGLL
                    qtemp_t(i,:) = -matmul(VT,qtemp_t(i,:))
                end do

                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = qtemp_t(:,j)
                end do
            else if (in == -1) then
                ! mpi interface
                mpi_e=mpi_ibool(is)
                mpi_n=mpi_interface(4,is)
    
                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
                end do
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)
    
                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)

                VT= matmul(T,matmul(APA,invT))
                do i=1,NGLL
                    qtemp_t(i,:) = matmul(VT,(qtemp_n(i,:)-qtemp_t(i,:)))
                end do
                do j=1,5
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (qtemp_t(:,j))
                end do
            end if
        end do
    end subroutine computeExactRiemannSF

    subroutine computeExactRiemannSFAnelastic(flux,q,qi,neighbor,vp,vs,rho,lambda,mu,face, mpi_interface,&
                                              mpi_ibool, mpi_ibt,ibt,ibn,nx,nz,freex,pml,set_pml)
        !Riemann fluxes for the Strong form of DG - anelastic
        integer, dimension(:) :: neighbor,face
        integer, dimension(:,:) :: ibt,ibn
        logical :: set_pml
        real(kind=CUSTOM_REAL), dimension(:,:) :: freex
        real(kind=CUSTOM_REAL), dimension(3*NGLL,3) :: flux
        real(kind=CUSTOM_REAL), dimension(:) :: nx,nz
        real(kind=CUSTOM_REAL), dimension(:,:) :: q
        real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
        integer, dimension(:,:) :: mpi_interface
        integer, dimension(:) :: mpi_ibool
        integer, dimension(:,:) :: mpi_ibt
        real(kind=CUSTOM_REAL), dimension(NGLL,5) :: qtemp_n,qtemp_t
        real(kind=CUSTOM_REAL), dimension(NGLL,3) :: atemp_t
        real(kind=CUSTOM_REAL) :: nxe,nze
        real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
        real(kind=CUSTOM_REAL), dimension(5,5) :: T,invT
        real(kind=CUSTOM_REAL), dimension(3,3) :: aT
        real(kind=CUSTOM_REAL), dimension(3,5) :: VT
        real(kind=CUSTOM_REAL), dimension(3,5) :: aAPA
        integer :: is,in,i,j
        integer :: mpi_e, mpi_n
        logical :: pml
        ! compute exact riemann problem, elementwise

        ! for all equal

        aAPA = getAnelasticAPA(vp,vs,rho,lambda,mu)
        flux=0.0

        do is=1,3 !surfs
            ! get neighbor element
            in = neighbor(is)
            if (in>0) then ! not boundary of the whole domain
                qtemp_t=q(ibt(:,is),:)
                qtemp_n=q(ibn(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)

                VT= matmul(aT,matmul(aAPA,invT))

                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,(qtemp_n(i,:)-qtemp_t(i,:)))
                end do

                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = atemp_t(:,j)
                end do
                ! boundaries
                ! "------------------------------------------------------------------"
            else if ((in==0).and.(face(is) == -1)) then ! free
                qtemp_t=q(ibt(:,is),:)
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)
                VT= matmul(aT,matmul(matmul(aAPA,freex),invT))

                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,qtemp_t(i,:))
                end do

                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (atemp_t(:,j))
                end do
                ! "------------------------------------------------------------------"
            else if ((in==0).and.face(is) == -2) then !absorb
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)
                qtemp_t=q(ibt(:,is),:)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)
                VT= matmul(aT,matmul(aAPA,invT))
                do i=1,NGLL
                    atemp_t(i,:) = -matmul(VT,qtemp_t(i,:))
                end do

                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = atemp_t(:,j)
                end do
            else if (in == -1) then
                ! mpi interface
                mpi_e=mpi_ibool(is)
                mpi_n=mpi_interface(4,is)

                qtemp_t=q(ibt(:,is),:)
                do i=1,NGLL
                    qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
                end do
                nxe=nx(is*NGLL)
                nze=nz(is*NGLL)

                T = getT(nxe,nze)
                invT = getInvT(nxe,nze)
                aT= T(1:3,1:3)
                VT= matmul(aT,matmul(aAPA,invT))
                do i=1,NGLL
                    atemp_t(i,:) = matmul(VT,(qtemp_n(i,:)-qtemp_t(i,:)))
                end do
                do j=1,3
                    flux(((is-1)*NGLL+1):(is*NGLL),j) = (atemp_t(:,j))
                end do
            end if
        end do
    end subroutine computeExactRiemannSFAnelastic
end module waveMod
