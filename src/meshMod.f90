!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Marc S. Boxberg (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
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
module meshMod
    ! modul to create meshes
    use constantsMod
    use nodesMod
    use triTrafoMod
    use vtkMod
    use dmatricesMod
    use vandermondeMod
    use geometricFactorsMod
    use sourceReceiverMod
    use liftMod
    use normalsMod
    use matrixMod
    use dtMod
    use triangleNeighborMod
    use parameterMod
    use externalModelMod
    !use anelasticMod
    use logoMod
    use errorMessage
    use materialsMod

    implicit none

    type :: meshVar
        sequence
        logical :: has_src = .false.
        logical :: has_rec = .false.

        integer :: nglob
        integer :: nelem
        integer :: nmat
        integer :: ncoord
        integer :: pinterfaces
        integer :: mpi_nn
        integer :: mpi_nnmax
        integer :: mpi_ne
        integer :: mpi_nemax
        integer :: nproc
        integer :: nelx
        integer :: nelz

        integer, dimension(3*NGLL) :: fmaskv
        integer, dimension(NGLL,3) :: fmask

        integer, dimension(:), pointer :: mat => null()
        integer, dimension(:), pointer :: pml
        integer, dimension(:), pointer :: mpi_icon => null()
        integer, dimension(:), pointer :: loc2glob_nodes => null()
        integer, dimension(:), pointer :: mpi_neighbor => null()
        integer, dimension(:), pointer :: mpi_ninterface => null()
        integer, dimension(:,:), pointer :: ibool => null()
        integer, dimension(:,:), pointer :: elem => null()
        integer, dimension(:,:), pointer :: neighbor => null()
        integer, dimension(:,:), pointer :: face => null()
        integer, dimension(:,:), pointer :: pmlface => null()
        integer, dimension(:,:), pointer :: mpi_ibool => null()
        integer, dimension(:,:,:), pointer :: ibn => null()
        integer, dimension(:,:,:), pointer :: ibt => null()
        integer, dimension(:,:,:), pointer :: mpi_interface => null()
        integer, dimension(:,:,:), pointer :: mpi_connection => null()
        integer, dimension(:,:,:), pointer :: mpi_ibt => null()

        real(kind=CUSTOM_REAL) :: dtfactor

        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Dr
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Ds
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Drw
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Dsw
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: invVdm
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Vdm
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Vdmr
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Vdms
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: mass
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL) :: lift

        real(kind=CUSTOM_REAL), dimension(:), pointer :: vx => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vz => null()

        real(kind=CUSTOM_REAL), dimension(:), pointer :: rx => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: rz => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: sx => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: sz => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: jacobian => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vp => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vs => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: rho => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: mu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: lambda => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: qp => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: qs => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vpu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vsu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: muu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: lambdau => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: coord => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: matval => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: ylambda => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: ymu => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: wl => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: nx => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: nz => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: sj => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: fscale => null()


    end type meshVar

    contains

    subroutine createRegularMesh(this,par,errmsg)
        type(error_message) :: errmsg
        type(parameterVar) :: par
        type (srcVar) :: src, tempsrc
        type (recVar) :: rec, temprec
        type(materialIndizes) :: matInd
        type(meshVar), dimension(:), allocatable :: db
        type(materialVar), dimension(:), allocatable :: mat
        type(attenuationVar), dimension(:), allocatable :: att
        type(meshVar), intent(out) :: this


        ! general
        character(len=80) :: filename
        character(len=21) :: myname = 'meshMod'

        integer :: dummy
        integer :: i, j, k, c, d, is, ie, in, n, l, m, iproc, ios
        integer :: ee
        integer :: mpi_ti, mpi_ni
        integer :: temp1, temp2
        integer :: nfree, nabsorb
        integer :: e1, e2, e3
        integer, dimension(Np) :: iv , ivn
        integer, dimension(NGLL) :: ibn2
        integer, dimension(NGLL) :: fmask1, fmask2, fmask3
        real(kind=CUSTOM_REAL) :: minGLL, temp
        real(kind=CUSTOM_REAL), dimension(Np) :: x, z, r, s
        real(kind=CUSTOM_REAL), dimension(Np) :: rx, sx, rz, sz, jacobian
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: tempM1, tempM2
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: sDT, fDT
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: free_nodes
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: absorb_nodes

        !logical :: ext
        integer, dimension(:,:), allocatable :: neighbortemp

        !metis
        integer :: nb_neigh
        integer :: edgecut
        integer :: num_glob, num_part
        integer :: size_glob2loc_nodes
        integer, dimension(3) :: nb_temp
        integer, dimension(3) :: loc_nodes
        integer, dimension(0:4) :: options
        integer, dimension(:), allocatable :: elmnts
        integer, dimension(:), allocatable :: xadj
        integer, dimension(:), allocatable :: adjncy
        integer, dimension(:), allocatable :: part
        integer, dimension(:), allocatable :: glob2loc_elmnts
        integer, dimension(:), allocatable :: num_loc
        integer, dimension(:), allocatable :: nodes_elmnts, nnodes_elmnts
        integer, dimension(:), allocatable :: glob2loc_nodes_nparts
        integer, dimension(:), allocatable :: glob2loc_nodes_parts
        integer, dimension(:), allocatable :: glob2loc_nodes
        integer, dimension(:), allocatable :: part_nodes, num_parts
        integer, dimension(:), allocatable :: inviboolv
        integer, dimension(:), allocatable :: tempv
        integer, dimension(:,:), allocatable :: invibool
        integer, dimension(:,:), allocatable :: icom
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: tempvr

        !source
        integer :: locnsrc
        integer, dimension(:), pointer :: locsrctype => null()
        integer, dimension(:), pointer :: locsrcstf => null()
        integer, dimension(:), allocatable :: locsrcnr
        integer, dimension(:,:), allocatable :: srcnr
        real(kind=CUSTOM_REAL), dimension(:), pointer :: locsrcf0 => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: locsrcfactor => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: locsrcangle_force => null()
        real(kind=custom_real), dimension(:), pointer :: locdelay => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locsrcM => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locsrcxz => null()
        character(len=80), dimension(:), pointer :: locsrcextwavelet => null()

        !receiver
        integer :: locnrec
        integer, dimension(:), pointer :: locrecnr => null()
        integer, dimension(:), pointer :: locrecnum => null()
        integer, dimension(:,:), allocatable :: recnum
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locrecxz => null()

        !external model EXPERIMENTAL
        integer :: em_nx, em_nz
        integer, dimension(:), allocatable :: iregion
        integer, dimension(:,:), allocatable :: em_ipolg
        real(kind=CUSTOM_REAL) :: em_xmin, em_zmin, em_dx, em_dz
        real(kind=CUSTOM_REAL) :: em_alam, em_amu, em_temp
        real(kind=CUSTOM_REAL) :: xs, zs
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: em_rhog
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: em_vpg
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: em_vsg
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: em_xg
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: em_zg

        ! LL LL 23.4.14 was 1e-4
        !very important to choose the right epsilon, to small can cause problems
        real(kind=CUSTOM_REAL) :: epsilon = 1e-3

        ! anelastic part
        logical :: attenuation

        ! debug
        logical :: debug

        ! test
        real(kind=CUSTOM_REAL) :: t1,t2

        call addTrace(errmsg,myname)

        if (par%log) call writeLogo()
        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "Start mesher"
        if (par%log) write(*,*)

        ! if we use MPI we will create databases
        if (par%nproc > 1) then
           allocate(db(par%nproc))
        end if

        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) " use external mesh "

        ! coordinates file
        if (par%log) write(*,*) "read coordinate file"
        filename=trim('cubit/coord')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif
        read(19,*) this%ncoord
        allocate(this%coord(2,this%ncoord))

        do i=1,this%ncoord
            read(19,*) dummy,this%coord(1,i),this%coord(2,i)
        end do
        close(19)

        ! meshfile
        if (par%log) write(*,*) "read meshfile file"
        filename=trim('cubit/mesh')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif

        read(19,*) this%nelem
        allocate(this%elem(3,this%nelem))
        do i=1,this%nelem
            read(19,*) dummy, this%elem(1,i),this%elem(2,i),this%elem(3,i)
        end do
        close(19)
        ! calculate nglob
        this%nglob = this%nelem*Np


        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "number of coordinates ",this%ncoord
        if (par%log) write(*,*) "number of elements ",this%nelem
        if (par%log) write(*,*) "number of points ",this%nglob

        ! allocate arrays
        call allocateMeshArrays(this)
        allocate(neighbortemp(3,this%nelem))
        allocate(sDT(this%nelem),fdt(this%nelem))

        !get neighbors
        if (par%log) write(*,*) "calculate neighbors"
        call triangulation_neighbor_triangles(3,this%nelem,this%elem,neighbortemp)

        ! we have to reorder the neighboring
        this%neighbor(1,:) = neighbortemp(3,:)
        this%neighbor(2,:) = neighbortemp(1,:)
        this%neighbor(3,:) = neighbortemp(2,:)

        ! materials values
        call setupMaterials(par, mat, matInd, att, this%nelem, iregion, errmsg)
        !Add some values to the mesh-type that have been read in at other places.
        this%nmat = par%matn
        this%mat = matInd%type

        ! free_nodes file
        if (par%log) write(*,*) "read free nodes file"
        filename=trim('cubit/free')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif
        read(19,*) nfree
        allocate(free_nodes(nfree))
        do i=1,nfree
           read(19,*) free_nodes(i)
        end do
        close(19)

        ! absorb_nodes file
        if (par%log) write(*,*) "read absorbing file"
        filename=trim('cubit/absorb')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            call print(errmsg)
            stop
        endif
        read(19,*) nabsorb
        allocate(absorb_nodes(nabsorb))
        do i=1,nabsorb
           read(19,*) absorb_nodes(i)
        end do
        close(19)

        !"--------------------------------------------------------------------------------------"
        ! start to create local elements
        !"--------------------------------------------------------------------------------------"

        ! get scaleDT
        call scaleDT(this%coord(1,:),this%coord(2,:),this%elem,this%nelem,sDT)

        ! get local element
        call nodes(balpha(NGLL),x,z)
        call xyToRs(x,z,r,s)
        minGLL=abs(r(1)-r(2))

        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "Minimum gll point distance :", minGLL

        ! precalc matrices
        call vdm2d(this%vdm,r,s)
        call invVdm2D(this%vdm,this%invvdm,0,errmsg)
        this%mass = matmul(transpose(this%invvdm),this%invvdm)
        call dmatrices2d(this%Dr,this%Ds,r,s,this%vdm,errmsg)
        call gradVdm2d(this%vdmr,this%Vdms,r,s)
        ! weak matrices
        tempM1=matmul(this%vdm,transpose(this%vdm))
        tempM2=matmul(this%vdm,transpose(this%vdmr))
        call invert(tempM1, errmsg)
        this%Drw=matmul(tempM2,tempM1)
        tempM2=matmul(this%vdm,transpose(this%vdms))
        this%Dsw=matmul(tempM2,tempM1)

        ! find boundary nodes in the local element
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(s(i)+1.0) < EPS) then
                fmask1(d)=c
                this%fmask(d,1)=c
                d=d+1
            end if
        end do
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(r(i)+s(i)) < EPS) then
                fmask2(d)=c
                this%fmask(d,2)=c
                d=d+1
            end if
        end do
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(r(i)+1.0) < EPS) then
                fmask3(d)=c
                this%fmask(d,3)=c
                d=d+1
            end if
        end do

        call lift2d(this%lift,this%fmask,this%vdm,r,s, errmsg)
        if (.level.errmsg == 2) then; call print(errmsg); stop; endif

        temp = 1e30
        do i=1,this%nelem
            if ( sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(2,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(2,i)))**2) < temp) then
                temp = sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(2,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(2,i)))**2)
            elseif ( sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(3,i)))**2) < temp) then
                temp = sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(3,i)))**2)
            elseif ( sqrt( (this%coord(1,this%elem(2,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(2,i)) - this%coord(2,this%elem(3,i)))**2) < temp) then
                temp = sqrt( (this%coord(1,this%elem(2,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(2,i)) - this%coord(2,this%elem(3,i)))**2)
            endif
        end do

        if (par%log)  write(*,*) "minimal edge length in the mesh: ", temp
        ! LL LL test to determine minimal epsilon for partitioning
        epsilon = temp/1e4

        !"--------------------------------------------------------------------------------------"
        ! tranformation from local to physical coordinates
        !"--------------------------------------------------------------------------------------"
        !if (par%nproc ==1) then
        c=0
        do i=1,this%nelem
            do j=1,NP
                c=c+1
                this%vx(c) = 0.5 * ( -(r(j)+s(j)) * this%coord(1,this%elem(1,i)) + (1+r(j)) * this%coord(1,this%elem(2,i)) +&
                  (1+s(j)) * this%coord(1,this%elem(3,i)))
                this%vz(c) = 0.5 * ( -(r(j)+s(j)) * this%coord(2,this%elem(1,i)) + (1+r(j)) * this%coord(2,this%elem(2,i)) +&
                  (1+s(j)) * this%coord(2,this%elem(3,i)))
                ! set up ibool
                this%ibool(j,i)=c
            end do
        end do

       ! get geometric factors
        do i=1,this%nelem
            iv=this%ibool(:,i)
            call geometricFactors2d(rx,sx,rz,sz,jacobian,this%vx(iv),this%vz(iv),this%Dr,this%Ds)
            if (par%debug) then
                do j=1,size(jacobian)
                    if (jacobian(j) .le. 0.0) then
                        call add(errmsg, 2, "Jacobian negative or zero! Abort...", myname)
                        if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                    end if
                end do
            endif
            this%rx(iv)=rx
            this%sx(iv)=sx
            this%rz(iv)=rz
            this%sz(iv)=sz
            this%jacobian(iv)=jacobian
        end do

        ! LL LL test
        if (par%debug) then
            do i=1,1
                iv=this%ibool(:,i)
                allocate(tempvr(NP))
                do j=1,NP
                    ivn(j)=j
                end do
                tempvr=0.0
                ! write src to vtk
                write(filename,"('out/test.vtk')")
                !    filename=trim(outpath)//trim(filename)
                call writeVtkNodesRealData(filename, this%vx(iv),tempvr,this%vz(iv), real(ivn))
                deallocate(tempvr)
            end do
        endif

        ! get normals
        do i=1,this%nelem
            iv=this%ibool(:,i)
            call normals2d(this%nx(:,i),this%nz(:,i),this%sj(:,i),this%Dr,this%Ds,this%vx(iv),this%vz(iv),this%fmask)
        end do

        ! make scaling vector
        do i=1,this%nelem
            iv=this%ibool(:,i)
            c=1
            do k=1,3
                do j=1,NGLL
                    d=iv(this%fmask(j,k))
                    this%fmaskv(c)=this%fmask(j,k)
                    this%fscale(c,i)=this%sj(c,i)/this%jacobian(d)
                    c=c+1
                end do
            end do
        end do
        ! par%nproc == 1


        ! clean neighbor array
        do i=1,this%nelem
            do j=1,3
                if (this%neighbor(j,i) <= 0) then
                    this%neighbor(j,i)=0
                end if
            end do
        end do

        !"--------------------------------------------------------------------------------------"
        ! build boundaries for external mesh ( very slow at the moment )
        !"--------------------------------------------------------------------------------------"
        this%face(:,:)=0
        ! find connections of faces
        do i=1,this%nelem
            do j=1,3
                c=this%neighbor(j,i)
                if (c>0) then
                    do k=1,3
                        d=this%neighbor(k,c)
                        if (d==i) then
                            this%face(j,i)=k
                        end if
                    end do
                else
                    ! set boundaries -1:free, -2:absorb
                    e1=this%elem(1,i)
                    e2=this%elem(2,i)
                    e3=this%elem(3,i)
                    !write(*,*) e1,e2,e3
                    n=0
                    m=0
                    do l=1,nabsorb
                        if (e1==absorb_nodes(l)) then
                            n=n+1
                        else if (e2==absorb_nodes(l)) then
                            n=n+1
                        else if (e3==absorb_nodes(l)) then
                            n=n+1
                        end if
                    end do
                    do l=1,nfree
                        if (e1==free_nodes(l)) then
                            m=m+1
                        else if (e2==free_nodes(l)) then
                            m=m+1
                        else if (e3==free_nodes(l)) then
                            m=m+1
                        end if
                    end do
                    if (n>1) then
                        this%face(j,i)=-2
                    else if (m>1) then
                        this%face(j,i)=-1
                    else
                        call add(errmsg, 2, "Error in boundaries!", myname)
                        if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                    end if
                end if
            end do
        end do

        this%pmlface=0
        do i=1,this%nelem
            do j=1,3
                c=this%neighbor(j,i)
                if(c>0) then
                    if (this%pml(c)>0) then
                        this%pmlface(i,j)=1
                    end if
                end if
            end do
        end do

        !"--------------------------------------------------------------------------------------"
        ! create boundray nodes association for single mesh
        !"--------------------------------------------------------------------------------------"
        if (iproc == 1) then
            do ie=1,this%nelem
                iv=this%ibool(:,ie)
                do is=1,3
                    in=this%neighbor(is,ie)
                    if (in>0) then
                        ivn=this%ibool(:,in)
                        this%ibt(:,is,ie)=iv(this%fmask(:,is))
                        j=this%face(is,ie)
                        this%ibn(:,is,ie)=ivn(this%fmask(:,j))

                        if ( (abs(this%vx(this%ibt(1,is,ie)) - this%vx(this%ibn(1,is,ie)) ) < epsilon)&
                            .and. (abs(this%vz(this%ibt(1,is,ie)) - this%vz(this%ibn(1,is,ie)) ) < epsilon)) then
                            this%ibn=this%ibn
                        else
                            k=NGLL
                            do j=1,NGLL
                                ibn2(k)=this%ibn(j,is,ie)
                                k=k-1
                            end do
                            this%ibn(:,is,ie)=ibn2
                        end if
                    else
                        this%ibt(:,is,ie)=iv(this%fmask(:,is))
                    end if
                end do
            end do
        end if

        !"--------------------------------------------------------------------------------------"
        ! create materials
        !"--------------------------------------------------------------------------------------"


        !!! LL LL open external model file EXPERIMANTAL
        if(par%extvel) then
            open(1,file = trim(par%externalfilename),status = 'old', iostat = ios)
            if (ios /= 0) then
                call add(errmsg,2,'could not open: '//trim(par%externalfilename),myname)
                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
            endif
            read (1,*) em_nx,em_nz,em_dx,em_dz,em_xmin,em_zmin
            allocate(em_xg(em_nx),em_zg(em_nz),em_ipolg(em_nx,em_nz),em_rhog(em_nx,em_nz),em_vpg(em_nx,em_nz),em_vsg(em_nx,em_nz))
            do j = 1,em_nx
                do k = 1,em_nz
                    read(1,*) em_ipolg(j,k),em_xg(j),em_zg(k),em_temp,em_rhog(j,k),em_vpg(j,k),em_vsg(j,k),em_alam,em_amu
                enddo
            enddo
            close(1)
            em_rhog = em_rhog*1000
            em_vpg = em_vpg*1000
            em_vsg = em_vsg*1000

            do i = 1, this%nelem
                xs=(this%coord(1,this%elem(1,i))+this%coord(1,this%elem(2,i))+this%coord(1,this%elem(3,i)))/3
                zs=(this%coord(2,this%elem(1,i))+this%coord(2,this%elem(2,i))+this%coord(2,this%elem(3,i)))/3
                call define_external_model(xs,zs,iregion(this%mat(i)),this%rho(i),this%vp(i),this%vs(i),em_nx,em_nz,em_dx,em_dz,em_xmin,em_zmin,em_ipolg,em_xg,em_zg,em_rhog,em_vpg,em_vsg )
                this%mu(i) = this%vs(i)**2 * this%rho(i)
                this%lambda(i) = this%vp(i)**2 * this%rho(i)-2*this%mu(i)
            end do
        else
            do i = 1, this%nelem
                this%rho(i) = mat(this%mat(i))%rho
                this%vp(i)  = mat(this%mat(i))%vp
                this%vs(i)  = mat(this%mat(i))%vs
                if(par%attenuation) then
                    this%vpu(i) = att(this%mat(i))%vp
                    this%vsu(i) = att(this%mat(i))%vs
                else
                    this%vpu(i) = 0.0
                    this%vsu(i) = 0.0
                end if
                this%qp(i)        = mat(this%mat(i))%qp
                this%qs(i)        = mat(this%mat(i))%qs
                this%ylambda(:,i) = att(this%mat(i))%ylambda(:)
                this%ymu(:,i)     = att(this%mat(i))%ymu(:)
                this%wl(:,i)      = att(this%mat(i))%omegaL(:)
                this%pml(i)       = matInd%pml(i)
                if (mat(this%mat(i))%vs .gt. 0) then
                    this%mu(i) = mat(this%mat(i))%mu
                    if (par%attenuation) then
                        this%muu(i) = att(this%mat(i))%mu
                    else
                        this%muu(i) = 0.0
                    end if
                else
                    ! acustic test
                    this%mu(i)=0
                end if
                this%lambda(i)=mat(this%mat(i))%lambda
                if (par%attenuation) then
                    this%lambdau(i)=att(this%mat(i))%lambda
                else
                    this%lambdau(i)= 0.0
                end if
            end do
        end if

        write(filename,"('out/vp_model.vtk')")
        call writeVtkTriMeshRealdata2D(filename,this%elem,this%coord,this%vp)
        if (par%set_pml) then
            write(filename,"('out/pml_model.vtk')")
            call writeVtkTriMeshRealdata2D(filename,this%elem,this%coord,real(this%pml))
        endif

        ! compute dtfactor for every element
        temp=1e7
        do i=1,this%nelem
            fDT(i)= (2.0/3.0)*minGLL*(sDT(i)/this%vp(i))
            if (fDT(i) < temp) then
                temp=fDT(i)
            end if
        end do
        this%dtfactor=temp
        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "dtfactor: ",this%dtfactor

        !"--------------------------------------------------------------------------------------"
        !"--------------------------------------------------------------------------------------"
        ! create partions with metis
        !"--------------------------------------------------------------------------------------"
        !"--------------------------------------------------------------------------------------"

        ! create partitioning
        if (par%nproc > 1) then
            ! build elements vector
            allocate(elmnts(this%nelem*3))
            k=1
            do i=1,this%nelem
                do j=1,3
                    elmnts(k)=this%elem(j,i)
                    k=k+1
                end do
            end do
            ! create compact neighbor array
            allocate(xadj(this%nelem+1))
            allocate(adjncy(1:max_neighbor*this%nelem))
            adjncy(:) = 0
            xadj(1)=1
            k=0
            do i=1,this%nelem
                nb_neigh=0
                nb_temp(:)=0
                do j=1,3
                    if (this%neighbor(j,i) > 0) then
                        nb_neigh=nb_neigh+1
                        nb_temp(nb_neigh)=this%neighbor(j,i)
                    end if
                end do
                l=1
                do j=k+1,k+nb_neigh
                    adjncy(j) = nb_temp(l)
                    l=l+1
                end do
                k=k+nb_neigh
                xadj(i+1)=k+1
            end do
            ! create metis partition
            allocate(part(this%nelem))
            options(:)=0

            call METIS_PartGraphRecursive(this%nelem, xadj, adjncy, 0, 0, 0, 1, par%nproc,options, edgecut, part)
            ! create glob2loc_elem
            ! inspired by the specfem partition
            ! be careful, metis gives parts begining from 0 or from 1 depending on the version, so i compile a local version, here 4.0.3
            allocate(glob2loc_elmnts(this%nelem))
            glob2loc_elmnts(:) = 0

            allocate(num_loc(par%nproc))
            num_loc(:) = 1

            ! build a vector with the local numbering of the elements
            do num_glob=1,this%nelem
                num_part=part(num_glob)
                glob2loc_elmnts(num_glob) = num_loc(num_part)
                num_loc(num_part) = num_loc(num_part) + 1
            end do

            !create local node numbering
            ! 2 verctors, nnodes with the positions of the elements in the vektor nodes_elements, similar to the metis numbering
            allocate(nnodes_elmnts(this%ncoord))
            allocate(nodes_elmnts(this%ncoord*nsize))
            nnodes_elmnts(:)=0
            nodes_elmnts(:)=0
            do i=1, 3*this%nelem
                ! nodes is like a matrix with notes as rows and nsize elements as colums
                nodes_elmnts((elmnts(i)-1)*nsize + nnodes_elmnts(elmnts(i))+1) = 1+(i-1)/3
                nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) +1
            end do

            ! now create the local node numbering
            allocate(glob2loc_nodes_nparts(this%ncoord+1))
            allocate(part_nodes(par%nproc), num_parts(par%nproc))

            size_glob2loc_nodes = 1
            part_nodes(:) = 0

            ! go through all coordinates
            do in=1,this%ncoord
                glob2loc_nodes_nparts(in) = size_glob2loc_nodes
                do ie = 1, nnodes_elmnts(in)
                    ! nodes_elmnts((ie)+nsize*(in-1)) gibt eins der maximal _nsize_ elemente zum knoten _in_
                    ! gibt an in welchen partitionen der knoten alles liegt
                    part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
                end do
                do num_part=1,par%nproc ! gehe durch die partitionen
                    if (part_nodes(num_part) == 1) then
                        ! get number of nodes in the global array, there might be some nodes at the interfaces doubble
                        size_glob2loc_nodes = size_glob2loc_nodes +1
                        part_nodes(num_part) = 0
                    end if
                end do
            end do
            glob2loc_nodes_nparts(this%ncoord+1) = size_glob2loc_nodes
            !Old allocation. Will be kept for now to confirm validity of the fix TM
            !allocate(glob2loc_nodes_parts(glob2loc_nodes_nparts(this%ncoord)))
            !allocate(glob2loc_nodes(glob2loc_nodes_nparts(this%ncoord)))
            allocate(glob2loc_nodes_parts(size_glob2loc_nodes))
            allocate(glob2loc_nodes(size_glob2loc_nodes))

            glob2loc_nodes(:) = 1
            part_nodes(:) = 0
            num_parts(:)=1
            size_glob2loc_nodes = 1

            do in=1,this%ncoord
                do ie = 1, nnodes_elmnts(in)
                    ! nodes_elmnts((ie)+nsize*(in-1)) gibt eins der maximal _nsize_ elemente zum knoten _in_
                    ! gibt an in welchen partitionen der knoten alles liegt
                    part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
                end do
                do num_part=1,par%nproc
                    if (part_nodes(num_part) == 1) then
                        ! build arrays with local nodes ordering
                        glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
                        glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
                        size_glob2loc_nodes = size_glob2loc_nodes +1
                        num_parts(num_part) = num_parts(num_part) +1
                        part_nodes(num_part) = 0
                    end if
                end do
            end do

            !"--------------------------------------------------------------------------------------"
            ! build Databases
            !"--------------------------------------------------------------------------------------"

            allocate(invibool(Np,this%nelem))
            allocate(inviboolv(this%nglob))
            invibool(:,:) = 0
            inviboolv(:)=0

            do iproc=1,par%nproc
                db(iproc)%nglob = (num_loc(iproc)-1)*Np
                db(iproc)%nelem = num_loc(iproc)-1
                db(iproc)%nmat = this%nmat
                db(iproc)%ncoord = num_parts(iproc)-1
                db(iproc)%dtfactor = this%dtfactor
                db(iproc)%invVdm = this%invVdm
                db(iproc)%Vdm = this%Vdm
                db(iproc)%Dr = this%Dr
                db(iproc)%Ds = this%Ds
                db(iproc)%Drw = this%Drw
                db(iproc)%Dsw = this%Dsw
                db(iproc)%Vdmr = this%Vdmr
                db(iproc)%Vdms = this%Vdms
                db(iproc)%mass = this%mass
                db(iproc)%lift = this%lift
                db(iproc)%fmask = this%fmask
                db(iproc)%fmaskv = this%fmaskv

                allocate(db(iproc)%vx(db(iproc)%nglob),db(iproc)%vz(db(iproc)%nglob))
                allocate(db(iproc)%matval(db(iproc)%nmat,7))
                allocate(db(iproc)%rx(db(iproc)%nglob), db(iproc)%rz(db(iproc)%nglob), &
                        db(iproc)%sx(db(iproc)%nglob), db(iproc)%sz(db(iproc)%nglob))
                allocate(db(iproc)%jacobian(db(iproc)%nglob))
                allocate(db(iproc)%mat(db(iproc)%nelem))
                allocate(db(iproc)%ibool(Np,db(iproc)%nelem))
                allocate(db(iproc)%nx(3*NGLL,db(iproc)%nelem), db(iproc)%nz(3*NGLL,db(iproc)%nelem), &
                        db(iproc)%sj(3*NGLL,db(iproc)%nelem),db(iproc)%fscale(3*NGLL,db(iproc)%nelem))
                allocate(db(iproc)%coord(2,db(iproc)%ncoord))
                allocate(db(iproc)%elem(3,db(iproc)%nelem))
                allocate(db(iproc)%neighbor(3,db(iproc)%nelem))
                allocate(db(iproc)%mpi_interface(4,3,db(iproc)%nelem))
                allocate(db(iproc)%face(3,db(iproc)%nelem))
                allocate(db(iproc)%pmlface(3,db(iproc)%nelem))
                allocate(db(iproc)%vp(db(iproc)%nelem),db(iproc)%vs(db(iproc)%nelem),db(iproc)%rho(db(iproc)%nelem),&
                        db(iproc)%qp(db(iproc)%nelem),db(iproc)%qs(db(iproc)%nelem),&
                        db(iproc)%mu(db(iproc)%nelem),db(iproc)%lambda(db(iproc)%nelem),db(iproc)%pml(db(iproc)%nelem))

                allocate(db(iproc)%vpu(db(iproc)%nelem),db(iproc)%vsu(db(iproc)%nelem),db(iproc)%muu(db(iproc)%nelem),db(iproc)%lambdau(db(iproc)%nelem))
                allocate(db(iproc)%ylambda(nMB,db(iproc)%nelem),db(iproc)%ymu(nMB,db(iproc)%nelem),db(iproc)%wl(nMB,db(iproc)%nelem))
                allocate(db(iproc)%ibt(NGLL,3,db(iproc)%nelem),db(iproc)%ibn(NGLL,3,db(iproc)%nelem))
                allocate(db(iproc)%mpi_ibt(NGLL,3,db(iproc)%nelem))
                allocate(db(iproc)%loc2glob_nodes(db(iproc)%ncoord))
                ! build local arrays
                ! build local coordinates
                do i=1, this%ncoord
                    do j= glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
                        if (glob2loc_nodes_parts(j) == iproc) then
                            k=glob2loc_nodes(j)
                            db(iproc)%coord(1,k)=this%coord(1,i)
                            db(iproc)%coord(2,k)=this%coord(2,i)
                            db(iproc)%loc2glob_nodes(k)=i
                        end if
                    end do
                end do

                ! build local elements
                do i=1, this%nelem
                    if (part(i) == iproc) then
                        do j=1,3
                            l=elmnts((i-1)*3+j)
                            m=elmnts(((i-1)*3+j))
                            do k=glob2loc_nodes_nparts(l), glob2loc_nodes_nparts(m+1)-1
                                if (glob2loc_nodes_parts(k) == iproc) then
                                    loc_nodes(j) = glob2loc_nodes(k)
                                end if
                            end do
                        end do
                        k = glob2loc_elmnts(i)
                        db(iproc)%elem(1,k) = loc_nodes(1)
                        db(iproc)%elem(2,k) = loc_nodes(2)
                        db(iproc)%elem(3,k) = loc_nodes(3)
                    end if
                end do

                !"--------------------------------------------------------------------------------------"
                ! find mpi neighbors and the neigbor element
                !"--------------------------------------------------------------------------------------"
                db(iproc)%neighbor(:,:) = 0
                db(iproc)%mpi_interface(:,:,:) = 0
                db(iproc)%pinterfaces=0
                do ie=1,this%nelem
                    do j=1,3
                        k=this%neighbor(j,ie)
                        if (k>0) then
                            if (part(ie) == iproc) then
                                if (part(k) == iproc) then
                                    l = glob2loc_elmnts(ie)
                                    m = glob2loc_elmnts(k)
                                    db(iproc)%neighbor(j,l) = m
                                else
                                    l = glob2loc_elmnts(ie)
                                    m = glob2loc_elmnts(k)
                                    db(iproc)%neighbor(j,l) = -1 ! means mpi neighbor
                                    db(iproc)%mpi_interface(1,j,l) = part(k) ! neighbor partition
                                    db(iproc)%mpi_interface(2,j,l) = m ! element in partition
                                    db(iproc)%pinterfaces=db(iproc)%pinterfaces+1
                                end if
                            end if
                        end if
                    end do
                end do

                !"--------------------------------------------------------------------------------------"
                ! recreate all local arrays
                !"--------------------------------------------------------------------------------------"
                db(iproc)%face(:,:) = 0

                c=0
                do i=1,db(iproc)%nelem
                    do j=1,NP
                        c=c+1
                        db(iproc)%vx(c) = 0.5 * ( -(r(j)+s(j)) * db(iproc)%coord(1,db(iproc)%elem(1,i)) + &
                        (1+r(j)) * db(iproc)%coord(1,db(iproc)%elem(2,i)) + (1+s(j)) * db(iproc)%coord(1,db(iproc)%elem(3,i)))
                        db(iproc)%vz(c) = 0.5 * ( -(r(j)+s(j)) * db(iproc)%coord(2,db(iproc)%elem(1,i)) + &
                        (1+r(j)) * db(iproc)%coord(2,db(iproc)%elem(2,i)) + (1+s(j)) * db(iproc)%coord(2,db(iproc)%elem(3,i)))
                        ! set up ibool
                        db(iproc)%ibool(j,i)=c
                    end do
                end do

                ! get geometric factors
                do i=1,db(iproc)%nelem
                    iv=db(iproc)%ibool(:,i)
                    call geometricFactors2d(rx,sx,rz,sz,jacobian,db(iproc)%vx(iv),db(iproc)%vz(iv),db(iproc)%Dr,db(iproc)%Ds)
                    if (par%debug) then
                        do j=1,size(jacobian)
                            if (jacobian(j) .le. 0.0) then
                                call add(errmsg, 2, "Jacobian negative or zero! Abort...", myname)
                                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                            end if
                        end do
                    endif
                    db(iproc)%rx(iv)=rx
                    db(iproc)%sx(iv)=sx
                    db(iproc)%rz(iv)=rz
                    db(iproc)%sz(iv)=sz
                    db(iproc)%jacobian(iv)=jacobian
                end do

                ! get normals
                do i=1,db(iproc)%nelem
                    iv=db(iproc)%ibool(:,i)
                    call normals2d(db(iproc)%nx(:,i),db(iproc)%nz(:,i),db(iproc)%sj(:,i),db(iproc)%Dr,db(iproc)%Ds,db(iproc)%vx(iv),db(iproc)%vz(iv),db(iproc)%fmask)
                end do

                ! make scaling vector
                do i=1,db(iproc)%nelem
                    iv=db(iproc)%ibool(:,i)
                    c=1
                    do k=1,3
                        do j=1,NGLL
                            d=iv(db(iproc)%fmask(j,k))
                            db(iproc)%fmaskv(c)=db(iproc)%fmask(j,k)
                            db(iproc)%fscale(c,i)=db(iproc)%sj(c,i)/db(iproc)%jacobian(d)
                            c=c+1
                        end do
                    end do
                end do

                ! set arrays to local element numbering
                do ie=1, this%nelem
                    if (part(ie)== iproc) then
                        l = glob2loc_elmnts(ie)
                        db(iproc)%face(:,l)=this%face(:,ie)
                        db(iproc)%pmlface(:,l)=this%pmlface(:,ie)
                        db(iproc)%mat(l)=this%mat(ie)
                        db(iproc)%vp(l)=this%vp(ie)
                        db(iproc)%vs(l)=this%vs(ie)
                        db(iproc)%rho(l)=this%rho(ie)
                        db(iproc)%qp(l)=this%qp(ie)
                        db(iproc)%qs(l)=this%qs(ie)

                        db(iproc)%vpu(l)=this%vpu(ie)
                        db(iproc)%vsu(l)=this%vsu(ie)
                        db(iproc)%lambdau(l)=this%lambdau(ie)
                        db(iproc)%muu(l)=this%muu(ie)

                        db(iproc)%ylambda(:,l)=this%ylambda(:,ie)
                        db(iproc)%ymu(:,l)=this%ymu(:,ie)
                        db(iproc)%wl(:,l)=this%wl(:,ie)
                        db(iproc)%lambda(l)=this%lambda(ie)
                        db(iproc)%mu(l)=this%mu(ie)
                        db(iproc)%pml(l)=this%pml(ie)
                    end if
                end do
            end do !iproc

            !"--------------------------------------------------------------------------------------"
            ! recreate boundary nodes association with local numbering
            !"--------------------------------------------------------------------------------------"
            do iproc=1,par%nproc
                db(iproc)%ibn(:,:,:) =0
                do ie=1,db(iproc)%nelem
                    iv=db(iproc)%ibool(:,ie)
                    do is=1,3
                        in=db(iproc)%neighbor(is,ie)
                        if (in>0) then
                            ivn=db(iproc)%ibool(:,in)
                            db(iproc)%ibt(:,is,ie)=iv(db(iproc)%fmask(:,is))
                            j=db(iproc)%face(is,ie)
                            db(iproc)%ibn(:,is,ie)=ivn(db(iproc)%fmask(:,j))
                            if ( .not. ( (abs(db(iproc)%vx(db(iproc)%ibt(1,is,ie)) - db(iproc)%vx(db(iproc)%ibn(1,is,ie)) ) < epsilon)&
                                .and. (abs(db(iproc)%vz(db(iproc)%ibt(1,is,ie)) - db(iproc)%vz(db(iproc)%ibn(1,is,ie)) ) < epsilon)) ) then
                                k=NGLL
                                do j=1,NGLL
                                    ibn2(k)=db(iproc)%ibn(j,is,ie)
                                    k=k-1
                                end do
                                db(iproc)%ibn(:,is,ie)=ibn2
                            end if
                        else if( in<=0) then
                            db(iproc)%ibt(:,is,ie)=iv(db(iproc)%fmask(:,is))
                        end if
                    end do
                end do
            end do
            !"--------------------------------------------------------------------------------------"
            ! create mpi interfaces
            !"--------------------------------------------------------------------------------------"
            allocate(icom(par%nproc,par%nproc))
            icom(:,:) = 0
            do iproc=1,par%nproc
                do ie=1,db(iproc)%nelem
                    do i=1,3
                        if (db(iproc)%mpi_interface(1,i,ie) > 0) then
                            icom(iproc,db(iproc)%mpi_interface(1,i,ie))=1
                        end if
                    end do
                end do
                db(iproc)%mpi_nn=sum(icom(iproc,:))
            end do

            ! set up mpi_neighbor array
            do iproc=1,par%nproc
                allocate(db(iproc)%mpi_neighbor(db(iproc)%mpi_nn))
                c=1
                do i=1,par%nproc
                    if (icom(iproc,i) == 1) then
                        db(iproc)%mpi_neighbor(c)=i
                        c=c+1
                    end if
                end do
            end do

            ! how many mpi interfaces do we have?
            do iproc=1,par%nproc
                allocate(tempv(db(iproc)%mpi_nn))
                allocate(db(iproc)%mpi_ninterface(db(iproc)%mpi_nn))
                tempv(:)=0
                c=1
                do i=1,db(iproc)%mpi_nn
                    in=db(iproc)%mpi_neighbor(i)
                    do ie=1,db(iproc)%nelem
                        do is=1,3
                            if ( (db(iproc)%neighbor(is,ie) == -1) .and. (db(iproc)%mpi_interface(1,is,ie) == in)) then
                                tempv(i)=tempv(i)+1
                            end if
                        end do
                    end do
                    c=c+1
                    db(iproc)%mpi_ninterface(i)=tempv(i)
                end do
                deallocate(tempv)
            end do

            if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log) write(*,*) "found total", sum(db(:)%mpi_nn), "mpi interfaces"

            temp2=0
            do iproc=1,par%nproc
                temp1=db(iproc)%nelem
                if (temp2<temp1) temp2=temp1
            end do

            ! set up mpi_ibool with max mpi neighbors
            do iproc=1,par%nproc
                db(iproc)%mpi_nemax=temp2
                allocate(db(iproc)%mpi_ibool(3,db(iproc)%mpi_nemax))
                db(iproc)%mpi_ibool(:,:)=0
            end do

            temp2=0
            do iproc=1,par%nproc
                temp1=maxval(db(iproc)%mpi_ninterface(:))
                if (temp2<temp1) temp2=temp1
            end do

            ! build interface array
            !"--------------------------------------------------------------------------------------"
            ! create main mpi reconnection arrays.
            !"--------------------------------------------------------------------------------------"
            do iproc=1,par%nproc
                db(iproc)%mpi_nnmax=maxval(db(:)%mpi_nn)
                allocate(db(iproc)%mpi_icon(db(iproc)%mpi_nnmax))
                db(iproc)%mpi_icon(:) = 0
                db(iproc)%mpi_ne=temp2
                ! arbeite mit max values
                allocate(db(iproc)%mpi_connection(db(iproc)%mpi_nn,db(iproc)%mpi_ne,2))
                db(iproc)%mpi_connection(:,:,:)=0
                do i=1,db(iproc)%mpi_nn
                    in=db(iproc)%mpi_neighbor(i)
                    c=1
                    do ie=1,db(iproc)%nelem
                        do is=1,3
                            if ( (db(iproc)%neighbor(is,ie) == -1) .and. (db(iproc)%mpi_interface(1,is,ie) == in)) then
                                db(iproc)%mpi_connection(i,c,1)=ie ! mpi interface to global element
                                db(iproc)%mpi_connection(i,c,2)=is ! which side of the element is the interface
                                db(iproc)%mpi_interface(3,is,ie)=c ! global element to local mpi interface
                                db(iproc)%mpi_interface(4,is,ie)=i ! which local interface?
                                c=c+1
                            end if
                        end do
                    end do
                end do
            end do

            !test plot of mesh silces
            if(par%debug) then
                do iproc=1,par%nproc
                    db(iproc)%nproc=par%nproc
                    write(filename,"('out/mesh',i6.6,'.ps')") iproc
                    call triangulation_order3_plot ( trim(filename), db(iproc)%ncoord, dble(db(iproc)%coord), &
                                                    db(iproc)%nelem, db(iproc)%elem, 2, 2 )
                end do
            endif

            do iproc = 1,par%nproc
                do i = 1,db(iproc)%mpi_nn
                    l = db(iproc)%mpi_neighbor(i)
                    do ie = 1,db(iproc)%mpi_ne
                        if ( db(iproc)%mpi_connection(i,ie,1) >0) then
                            is = db(iproc)%mpi_connection(i,ie,2)
                            ee = db(iproc)%mpi_connection(i,ie,1)
                            k = db(iproc)%mpi_interface(2,is,ee)
                            c = db(l)%mpi_interface(3,db(iproc)%face(is,ee),k)
                            db(iproc)%mpi_icon(i)=db(l)%mpi_interface(4,db(iproc)%face(is,ee),k)
                            !very important for correct mpi boundarys
                            db(iproc)%mpi_ibool(is,ee) = c
                            do j=1,NGLL
                                mpi_ti = db(iproc)%ibt( j,is,ee)
                                mpi_ni = db(l)%ibt(j,db(iproc)%face(is,ee),k)
                                ! LL LL 24.4.14
                                if ((abs(db(iproc)%vx(mpi_ti) - db(l)%vx(mpi_ni) ) < epsilon)&
                                .and. (abs(db(iproc)%vz(mpi_ti) - db(l)%vz(mpi_ni))  < epsilon)) then
                                    db(l)%mpi_ibt(j,db(iproc)%face(is,ee),k)=j
                                else
                                    db(l)%mpi_ibt(j,db(iproc)%face(is,ee),k)=NGLL+1-j
                                end if
                            end do
                            do j=1,NGLL
                                mpi_ti = db(iproc)%ibt( j,is,ee)
                                mpi_ni = db(l)%ibt(db(l)%mpi_ibt(j,db(iproc)%face(is,ee),k),db(iproc)%face(is,ee),k)
                            end do
                        end if
                    end do
                end do
            end do

            !"--------------------------------------------------------------------------------------"
            ! find source and receiver in the global mesh
            !"--------------------------------------------------------------------------------------"
            call initSource(par,src, errmsg)
            call initReceiver(par,rec, errmsg)
            !find sources in global mesh
            call findSource(par,src,this%vx,this%vz,this%nglob,this%nelem,this%ibool,this%coord,this%elem, this%Dr, this%Ds, errmsg)
            ! get global src elements
            allocate(tempv(par%nproc))
            allocate(srcnr(par%nproc,size(src%srcelem)))
            srcnr(:,:)=0
            tempv(:)=0
            j=1
            do i=1,size(src%srcelem)
                db(part(src%srcelem(i)))%has_src=.true.
                tempv(part(src%srcelem(i)))=tempv(part(src%srcelem(i)))+1
                srcnr(part(src%srcelem(i)),tempv(part(src%srcelem(i))))=j
                j=j+1
            end do

            !change element number from global to local numerbing
            do i=1,size(src%srcelem)
                call changeSrcElement(src,i,glob2loc_elmnts(src%srcelem(i)))
            end do
            ! create local source var for proc that carries the source
            j=1
            do iproc=1,par%nproc
                if (db(iproc)%has_src) then
                    !recreate src arrays
                    allocate(locsrcxz(2,tempv(iproc)))
                    allocate(locsrcnr(tempv(iproc)))
                    allocate(locsrctype(tempv(iproc)))
                    allocate(locsrcstf(tempv(iproc)))
                    allocate(locsrcextwavelet(tempv(iproc)))
                    allocate(locsrcf0(tempv(iproc)))
                    allocate(locsrcfactor(tempv(iproc)))
                    allocate(locsrcangle_force(tempv(iproc)))
                    allocate(locsrcm(3,tempv(iproc)))
                    allocate(locdelay(tempv(iproc)))
                    locnsrc=tempv(iproc)
                    do i=1,tempv(iproc)
                        locsrcnr(i) = srcnr(iproc,i)
                        locsrcxz(:,i) = src%srcxz(:,locsrcnr(i))
                        locsrctype(i) = src%srctype(locsrcnr(i))
                        locsrcstf(i) = src%srcstf(locsrcnr(i))
                        locsrcextwavelet(i) = src%extwavelet(locsrcnr(i))
                        locsrcf0(i) = src%srcf0(locsrcnr(i))
                        locsrcfactor(i) = src%srcfactor(locsrcnr(i))
                        locsrcangle_force(i) = src%srcangle_force(locsrcnr(i))
                        locsrcM(:,i) = src%srcm(:,locsrcnr(i))
                        locdelay(i) = src%delay(locsrcnr(i))
                        j=j+1
                    end do
                    call prepareRecalcSrc(tempsrc,locnsrc,locsrcxz,locsrctype,locsrcstf,locsrcextwavelet,locsrcf0,locsrcfactor,locsrcangle_force,locsrcM,locdelay)
                    call findSource(par,tempsrc,db(iproc)%vx,db(iproc)%vz,db(iproc)%nglob,db(iproc)%nelem,db(iproc)%ibool, &
                                    db(iproc)%coord,db(iproc)%elem,db(iproc)%dr,db(iproc)%ds, errmsg)
                    write(filename,"('out/srcVar',i6.6)") iproc
                    !filename=trim(outpath)//trim(filename)
                    call writeSrcVar(tempsrc,filename)
                    deallocate(locsrcnr)
                    deallocate(locsrcxz)
                    deallocate(locsrctype,locsrcstf,locsrcextwavelet,locsrcf0,locsrcfactor,locsrcangle_force,locsrcM)
                    call deallocSrcVar(tempsrc)
                end if
            end do
            deallocate(tempv)
            deallocate(srcnr)

            write(*,*) "start finding receiver"
            call findReceiver(par,rec,this%vx,this%vz,this%nglob,this%nelem,this%ibool,this%coord,this%elem, this%Dr, this%Ds, errmsg)
            write(*,*) "end finding receiver"
            ! get global rec elements
            allocate(tempv(par%nproc))
            allocate(recnum(par%nproc,size(rec%recelem)))
            recnum(:,:)=0
            tempv(:)=0
            j=1
            do i=1,size(rec%recelem)
                db(part(rec%recelem(i)))%has_rec=.true.
                tempv(part(rec%recelem(i)))=tempv(part(rec%recelem(i)))+1
                recnum(part(rec%recelem(i)),tempv(part(rec%recelem(i))))=j
                j=j+1
            end do

            !change element number from global to local numerbing
            do i=1,size(rec%recelem)
                call changeRecElement(rec,i,glob2loc_elmnts(rec%recelem(i)))
            end do
            ! create local source var for proc that carries the source
            j=1
            do iproc=1,par%nproc
                if (db(iproc)%has_rec) then
                    ! recreate rec arrays
                    allocate(locrecxz(2,tempv(iproc)))
                    allocate(locrecnr(tempv(iproc)))
                    allocate(locrecnum(tempv(iproc)))
                    locnrec=tempv(iproc)
                    do i=1,tempv(iproc)
                        locrecnum(i) = recnum(iproc,i)
                        locrecxz(:,i) = rec%recxz(:,locrecnum(i))
                        locrecnr(i) = rec%recnr(locrecnum(i))
                        j=j+1
                    end do
                    call prepareRecalcRec(temprec,locnrec,locrecxz,locrecnr)
                    call findReceiver(par,temprec,db(iproc)%vx,db(iproc)%vz,db(iproc)%nglob,db(iproc)%nelem,db(iproc)%ibool,&
                                      db(iproc)%coord,db(iproc)%elem,db(iproc)%dr,db(iproc)%ds, errmsg)
                    write(filename,"('out/recVar',i6.6)") iproc
                    call writeRecVar(temprec,filename)
                    deallocate(locrecnum)
                    deallocate(locrecxz)
                    deallocate(locrecnr)
                    call deallocRecVar(temprec)
                end if
            end do
            deallocate(tempv)
            allocate(tempvr(size(src%srcxz(1,:))))
            tempvr=0.0
            !write src to vtk
            write(filename,"('out/src.vtk')")
            call writeVtkNodes(filename, src%srcxz(1,:),tempvr,src%srcxz(2,:))
            deallocate(tempvr)

            allocate(tempvr(size(rec%recxz(1,:))))
            tempvr=0.0
            ! write receiver to vtk
            write(filename,"('out/rec.vtk')")
            call writeVtkNodes(filename, rec%recxz(1,:),tempvr,rec%recxz(2,:))
            deallocate(tempvr)

            ! write DATABASE
            do iproc=1,par%nproc
                write(filename,"('out/meshVar',i6.6)") iproc
                call writeMeshVar(db(iproc),filename)
                call deallocMeshVar(db(iproc))
                if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
                if (par%log) write(*,*) "done writing database files"
            end do !nproc create databases
        end if ! nproc >1

        !"--------------------------------------------------------------------------------------"
        ! deallocate arrays
        !"--------------------------------------------------------------------------------------"

        deallocate(sDT,fDT)

        if (par%nproc >1) then
            deallocate(elmnts)
            deallocate(xadj)
            deallocate(adjncy)
            deallocate(part)
            deallocate(glob2loc_elmnts)
            deallocate(num_loc)
            deallocate(nnodes_elmnts)
            deallocate(nodes_elmnts)
            deallocate(glob2loc_nodes_nparts)
            deallocate(part_nodes, num_parts)
            deallocate(glob2loc_nodes_parts)
            deallocate(glob2loc_nodes)
            do iproc=1,par%nproc
                call deallocMeshVar(db(iproc))
            end do
        end if
        write(*,*) "dealloc src"
        call deallocSrcVar(src)
        write(*,*) "dealloc rec"
        call deallocRecVar(rec)
    end subroutine createRegularMesh

    !"--------------------------------------------------------------------------------------"
    ! destructor for meshvar
    !
    subroutine deallocMeshvar(this)
        type(meshVar) :: this
        if (associated(this%matval)) deallocate(this%matval)
        if (associated(this%vx)) deallocate(this%vx)
        if (associated(this%vz)) deallocate(this%vz)
        if (associated(this%rx)) deallocate(this%rx)
        if (associated(this%rz)) deallocate(this%rz)
        if (associated(this%sx)) deallocate(this%sx)
        if (associated(this%sz)) deallocate(this%sz)
        if (associated(this%nx)) deallocate(this%nx)
        if (associated(this%nz)) deallocate(this%nz)
        if (associated(this%sJ)) deallocate(this%sJ)
        if (associated(this%fscale)) deallocate(this%fscale)
        if (associated(this%jacobian)) deallocate(this%jacobian)
        if (associated(this%mat)) deallocate(this%mat)
        if (associated(this%ibool)) deallocate(this%ibool)
        if (associated(this%coord)) deallocate(this%coord)
        if (associated(this%elem)) deallocate(this%elem)
        if (associated(this%neighbor)) deallocate(this%neighbor)
        if (associated(this%face)) deallocate(this%face)
        if (associated(this%pmlface)) deallocate(this%pmlface)
        if (associated(this%vp)) deallocate(this%vp)
        if (associated(this%vs)) deallocate(this%vs)
        if (associated(this%rho)) deallocate(this%rho)
        if (associated(this%qp)) deallocate(this%qp)
        if (associated(this%qs)) deallocate(this%qs)
        if (associated(this%mu)) deallocate(this%mu)
        if (associated(this%vpu)) deallocate(this%vpu)
        if (associated(this%vsu)) deallocate(this%vsu)
        if (associated(this%lambdau)) deallocate(this%lambdau)
        if (associated(this%muu)) deallocate(this%muu)
        if (associated(this%lambda)) deallocate(this%lambda)
        if (associated(this%ymu)) deallocate(this%ymu)
        if (associated(this%wl)) deallocate(this%wl)
        if (associated(this%ylambda)) deallocate(this%ylambda)
        if (associated(this%loc2glob_nodes)) deallocate(this%loc2glob_nodes)
        if (associated(this%pml)) deallocate(this%pml)
        if (associated(this%mpi_interface)) deallocate(this%mpi_interface)
        if (associated(this%mpi_neighbor)) deallocate(this%mpi_neighbor)
        if (associated(this%mpi_connection)) deallocate(this%mpi_connection)
        if (associated(this%mpi_ninterface)) deallocate(this%mpi_ninterface)
        if (associated(this%mpi_ibool)) deallocate(this%mpi_ibool)
        if (associated(this%mpi_ibt)) deallocate(this%mpi_ibt)
        if (associated(this%mpi_icon)) deallocate(this%mpi_icon)
    end subroutine deallocMeshvar

    subroutine allocateMeshArrays(this)
        type(meshVar) :: this
        !"--------------------------------------------------------------------------------------"
        allocate(this%vx(this%nglob),this%vz(this%nglob))
        allocate(this%rx(this%nglob), this%rz(this%nglob), this%sx(this%nglob), this%sz(this%nglob))
        allocate(this%jacobian(this%nglob))
        allocate(this%mat(this%nelem))
        allocate(this%ibool(Np,this%nelem))
        allocate(this%nx(3*NGLL,this%nelem), this%nz(3*NGLL,this%nelem), this%sj(3*NGLL,this%nelem), this%fscale(3*NGLL,this%nelem))
        allocate(this%neighbor(3,this%nelem))
        allocate(this%face(3,this%nelem))
        allocate(this%pmlface(3,this%nelem))
        allocate(this%ibt(NGLL,3,this%nelem),this%ibn(NGLL,3,this%nelem))
        allocate(this%vp(this%nelem),this%vs(this%nelem),this%rho(this%nelem),this%qp(this%nelem),this%qs(this%nelem),this%mu(this%nelem),this%lambda(this%nelem),this%pml(this%nelem))
        allocate(this%vpu(this%nelem),this%vsu(this%nelem),this%muu(this%nelem),this%lambdau(this%nelem))
        allocate(this%ylambda(nMB,this%nelem),this%ymu(nMB,this%nelem), this%wl(nMB,this%nelem))
    end subroutine allocateMeshArrays

    !"--------------------------------------------------------------------------------------"
    ! write meshVar to file
    !
    subroutine writeMeshVar(this,filename)
        implicit none
        type(meshVar) :: this
        character(len=80) filename
        integer :: nglob, nelem, nmat, ncoord
        integer :: pinterfaces
        integer :: mpi_nn
        integer :: mpi_nnmax
        integer :: mpi_ne
        integer :: mpi_nemax
        integer :: nproc
        logical :: has_src
        logical :: has_rec
        real(kind=CUSTOM_REAL) :: dtfactor
        integer :: nelx,nelz

        real(kind=CUSTOM_REAL), dimension(Np,Np) :: invVdm,Vdm,Dr,Ds
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Drw,Dsw,Vdmr,Vdms
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: mass
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL) :: lift
        integer, dimension(NGLL,3) :: fmask
        integer, dimension(3*NGLL) :: fmaskv

        real(kind=CUSTOM_REAL), allocatable, dimension(:) :: vx ,vz
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: matval
        real(kind=CUSTOM_REAL), allocatable, dimension(:) :: rx ,rz ,sx ,sz
        integer, allocatable, dimension(:,:) :: ibool
        real(kind=CUSTOM_REAL), allocatable, dimension(:) :: jacobian
        integer(kind=CUSTOM_REAL), allocatable, dimension(:) :: mat
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: coord
        integer, allocatable, dimension(:,:) :: elem
        integer, allocatable, dimension(:,:) :: neighbor
        integer, allocatable, dimension(:,:) :: face
        integer, allocatable, dimension(:,:) :: pmlface
        real(kind=CUSTOM_REAL), allocatable, dimension(:) :: vp,vs,rho,mu,lambda,qp,qs
        real(kind=CUSTOM_REAL), allocatable, dimension(:) :: vpu,vsu,muu,lambdau
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: ymu,ylambda,wl
        integer, dimension(:), allocatable :: pml
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: nx,nz,sj,fscale
        integer, dimension(:,:,:), allocatable :: ibn ,ibt
        integer, dimension(:), allocatable :: loc2glob_nodes
        integer, allocatable, dimension(:,:,:) :: mpi_interface
        integer, dimension(:), allocatable :: mpi_neighbor
        integer, dimension(:,:,:), allocatable :: mpi_connection
        integer, dimension(:), allocatable :: mpi_ninterface
        integer, dimension(:,:), allocatable :: mpi_ibool
        integer, dimension(:,:,:), allocatable :: mpi_ibt
        integer, dimension(:), allocatable :: mpi_icon

        ! copy meshVar and save to file
        nglob = this%nglob
        nelem = this%nelem
        nmat = this%nmat
        ncoord = this%ncoord
        pinterfaces = this%pinterfaces
        mpi_nn=this%mpi_nn
        mpi_nnmax=this%mpi_nnmax
        mpi_ne=this%mpi_ne
        mpi_nemax=this%mpi_nemax
        nproc=this%nproc
        has_src=this%has_src
        has_rec=this%has_rec
        dtfactor = this%dtfactor
        nelx = this%nelx
        nelz = this%nelz

        invVdm = this%invVdm
        Vdm = this%Vdm
        Dr = this%Dr
        Ds = this%Ds
        Drw = this%Drw
        Dsw = this%Dsw
        Vdmr = this%Vdmr
        Vdms = this%Vdms
        mass = this%mass
        lift = this%lift
        fmask = this%fmask
        fmaskv = this%fmaskv

        allocate(vx(nglob),vz(nglob))
        allocate(matval(nmat,7))
        allocate(rx(nglob), rz(nglob), sx(nglob), sz(nglob))
        allocate(jacobian(nglob))
        allocate(mat(nelem))
        allocate(ibool(Np,nelem))
        allocate(nx(3*NGLL,nelem), nz(3*NGLL,nelem), sj(3*NGLL,nelem), fscale(3*NGLL,nelem))
        allocate(coord(2,ncoord))
        allocate(elem(3,nelem))
        allocate(neighbor(3,nelem))
        allocate(face(3,nelem))
        allocate(pmlface(3,nelem))
        allocate(vp(nelem),vs(nelem),rho(nelem),qp(nelem),qs(nelem),mu(nelem),lambda(nelem),pml(nelem))
        allocate(vpu(nelem),vsu(nelem),muu(nelem),lambdau(nelem))
        allocate(ymu(nMB,nelem),ylambda(nMB,nelem),wl(nMB,nelem))
        allocate(loc2glob_nodes(ncoord))
        allocate(ibt(NGLL,3,nelem),ibn(NGLL,3,nelem))
        allocate(mpi_interface(4,3,nelem))
        allocate(mpi_neighbor(mpi_nn))
        allocate(mpi_connection(mpi_nn,mpi_ne,2))
        allocate(mpi_ninterface(mpi_nn))
        allocate(mpi_ibool(3,mpi_nemax))
        allocate(mpi_ibt(NGLL,3,nelem))
        allocate(mpi_icon(mpi_nnmax))

        vx = this%vx
        vz = this%vz
        matval = this%matval
        rx = this%rx
        rz = this%rz
        sx = this%sx
        sz = this%sz
        jacobian = this%jacobian
        mat = this%mat
        ibool = this%ibool
        nx = this%nx
        nz = this%nz
        sj = this%sj
        fscale = this%fscale
        coord = this%coord
        elem = this%elem
        neighbor = this%neighbor
        face = this%face
        pmlface = this%pmlface
        vp = this%vp
        vs = this%vs
        rho = this%rho
        qp = this%qp
        qs = this%qs
        mu = this%mu
        vpu = this%vpu
        vsu = this%vsu
        lambdau = this%lambdau
        muu = this%muu
        lambda = this%lambda
        ymu = this%ymu
        ylambda = this%ylambda
        wl = this%wl
        pml = this%pml
        ibt = this%ibt
        ibn = this%ibn
        loc2glob_nodes = this%loc2glob_nodes
        mpi_interface = this%mpi_interface
        mpi_neighbor = this%mpi_neighbor
        mpi_connection = this%mpi_connection
        mpi_ninterface = this%mpi_ninterface
        mpi_ibool = this%mpi_ibool
        mpi_ibt = this%mpi_ibt
        mpi_icon = this%mpi_icon

        write(*,*) "write meshVar"

        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        write(27) nglob
        write(27) nelem
        write(27) nmat
        write(27) ncoord
        write(27) pinterfaces
        write(27) mpi_nn
        write(27) mpi_nnmax
        write(27) mpi_ne
        write(27) mpi_nemax
        write(27) nproc
        write(27) has_src
        write(27) has_rec
        write(27) dtfactor
        write(27) nelx
        write(27) nelz

        write(27) invVdm
        write(27) Vdm
        write(27) Dr
        write(27) Ds
        write(27) Drw
        write(27) Dsw
        write(27) Vdmr
        write(27) Vdms
        write(27) mass
        write(27) lift
        write(27) fmask
        write(27) fmaskv

        write(27) vx
        write(27) vz
        write(27) matval
        write(27) rx
        write(27) rz
        write(27) sx
        write(27) sz
        write(27) jacobian
        write(27) mat
        write(27) ibool
        write(27) nx
        write(27) nz
        write(27) sj
        write(27) fscale
        write(27) coord
        write(27) elem
        write(27) neighbor
        write(27) face
        write(27) pmlface
        write(27) vp
        write(27) vs
        write(27) rho
        write(27) qp
        write(27) qs
        write(27) mu
        write(27) lambda
        write(27) vpu
        write(27) vsu
        write(27) muu
        write(27) lambdau
        write(27) ymu
        write(27) ylambda
        write(27) wl
        write(27) pml
        write(27) ibt
        write(27) ibn
        write(27) loc2glob_nodes
        write(27) mpi_interface
        write(27) mpi_neighbor
        write(27) mpi_connection
        write(27) mpi_ninterface
        write(27) mpi_ibool
        write(27) mpi_ibt
        write(27) mpi_icon

        close(27)

        write(*,*) "done writing meshVar"

        deallocate(vx,vz)
        deallocate(matval)
        deallocate(rx, rz, sx, sz)
        deallocate(jacobian)
        deallocate(mat)
        deallocate(ibool)
        deallocate(nx, nz, sj, fscale)
        deallocate(coord)
        deallocate(elem)
        deallocate(neighbor)
        deallocate(face)
        deallocate(pmlface)
        deallocate(vp,vs,rho,qp,qs,mu,lambda,pml,ymu,ylambda,wl)
        deallocate(vpu,vsu,muu,lambdau)
        deallocate(ibt,ibn)
        deallocate(loc2glob_nodes)
        deallocate(mpi_interface)
        deallocate(mpi_neighbor)
        deallocate(mpi_connection)
        deallocate(mpi_ninterface)
        deallocate(mpi_ibool)
        deallocate(mpi_ibt)
        deallocate(mpi_icon)

    end subroutine writeMeshVar
    !"--------------------------------------------------------------------------------------"
    ! read meshVar from file
    !
    subroutine readMeshVar(this,filename, errmsg)
        type(error_message) :: errmsg
        type(meshVar) :: this
        integer :: ios
        character(len=80) :: filename
        character(len=11) :: myname ="readMeshVar"

        call addTrace(errmsg, myname)

        write(*,*) "read meshVar" , filename
        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open: '//trim(filename),myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif
        read(27) this%nglob
        read(27) this%nelem
        read(27) this%nmat
        read(27) this%ncoord
        read(27) this%pinterfaces
        read(27) this%mpi_nn
        read(27) this%mpi_nnmax
        read(27) this%mpi_ne
        read(27) this%mpi_nemax
        read(27) this%nproc
        read(27) this%has_src
        read(27) this%has_rec
        read(27) this%dtfactor
        read(27) this%nelx
        read(27) this%nelz

        read(27) this%invVdm
        read(27) this%Vdm
        read(27) this%Dr
        read(27) this%Ds
        read(27) this%Drw
        read(27) this%Dsw
        read(27) this%Vdmr
        read(27) this%Vdms
        read(27) this%mass
        read(27) this%lift
        read(27) this%fmask
        read(27) this%fmaskv

        allocate(this%vx(this%nglob),this%vz(this%nglob))
        allocate(this%matval(this%nmat,7))
        allocate(this%rx(this%nglob), this%rz(this%nglob), this%sx(this%nglob), this%sz(this%nglob))
        allocate(this%jacobian(this%nglob))
        allocate(this%mat(this%nelem))
        allocate(this%ibool(Np,this%nelem))
        allocate(this%nx(3*NGLL,this%nelem), this%nz(3*NGLL,this%nelem), this%sj(3*NGLL,this%nelem),&
             this%fscale(3*NGLL,this%nelem))
        allocate(this%coord(2,this%ncoord))
        allocate(this%elem(3,this%nelem))
        allocate(this%neighbor(3,this%nelem))
        allocate(this%face(3,this%nelem))
        allocate(this%pmlface(3,this%nelem))
        allocate(this%vp(this%nelem),this%vs(this%nelem),this%rho(this%nelem),this%qp(this%nelem),this%qs(this%nelem),this%mu(this%nelem),this%lambda(this%nelem),this%pml(this%nelem))
        allocate(this%vpu(this%nelem),this%vsu(this%nelem),this%muu(this%nelem),this%lambdau(this%nelem))
        allocate(this%ymu(nMB,this%nelem),this%ylambda(nMB,this%nelem), this%wl(nMB,this%nelem))
        allocate(this%ibt(NGLL,3,this%nelem),this%ibn(NGLL,3,this%nelem))
        allocate(this%loc2glob_nodes(this%ncoord))
        allocate(this%mpi_interface(4,3,this%nelem))
        allocate(this%mpi_neighbor(this%mpi_nn))
        allocate(this%mpi_connection(this%mpi_nn,this%mpi_ne,2))
        allocate(this%mpi_ninterface(this%mpi_nn))
        allocate(this%mpi_ibool(3,this%mpi_nemax))
        allocate(this%mpi_ibt(NGLL,3,this%nelem))
        allocate(this%mpi_icon(this%mpi_nnmax))

        read(27) this%vx
        read(27) this%vz
        read(27) this%matval
        read(27) this%rx
        read(27) this%rz
        read(27) this%sx
        read(27) this%sz
        read(27) this%jacobian
        read(27) this%mat
        read(27) this%ibool
        read(27) this%nx
        read(27) this%nz
        read(27) this%sj
        read(27) this%fscale
        read(27) this%coord
        read(27) this%elem
        read(27) this%neighbor
        read(27) this%face
        read(27) this%pmlface
        read(27) this%vp
        read(27) this%vs
        read(27) this%rho
        read(27) this%qp
        read(27) this%qs
        read(27) this%mu
        read(27) this%lambda
        read(27) this%vpu
        read(27) this%vsu
        read(27) this%muu
        read(27) this%lambdau
        read(27) this%ymu
        read(27) this%ylambda
        read(27) this%wl
        read(27) this%pml
        read(27) this%ibt
        read(27) this%ibn
        read(27) this%loc2glob_nodes
        read(27) this%mpi_interface
        read(27) this%mpi_neighbor
        read(27) this%mpi_connection
        read(27) this%mpi_ninterface
        read(27) this%mpi_ibool
        read(27) this%mpi_ibt
        read(27) this%mpi_icon

        close(27)
    end subroutine readMeshVar
end module meshMod
