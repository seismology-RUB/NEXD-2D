!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
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
    use logoMod
    use errorMessage
    use materialsMod
    use slipInterfaceMod

    implicit none

    type :: meshVar
        sequence
        logical :: has_src = .false.
        logical :: has_rec = .false.
        logical :: poroelastic !poro

        integer :: nglob
        integer :: nelem
        integer :: nmat
        integer :: ncoord
        integer :: pinterfaces
        integer :: mpi_nn                       !Number of neighbouring partitions/mpi-interfaces in Processor
        integer :: mpi_nnmax                    !Maximum number of neighbouring partitions for all Processors
        integer :: mpi_ne
        integer :: mpi_nemax
        integer :: nproc                        !Number of Processors
        integer :: nelx
        integer :: nelz
        integer :: nfluids !poro

        integer, dimension(3*NGLL) :: fmaskv
        integer, dimension(NGLL,3) :: fmask

        integer, dimension(:), pointer :: mat => null()
        integer, dimension(:), pointer :: pml
        integer, dimension(:), pointer :: mpi_icon => null()
        integer, dimension(:), pointer :: loc2glob_nodes => null()
        integer, dimension(:), pointer :: mpi_neighbor => null()
        integer, dimension(:), pointer :: mpi_ninterface => null()
        integer, dimension(:), allocatable :: elemNo
        integer, dimension(:,:), pointer :: ibool => null()
        integer, dimension(:,:), pointer :: elem => null()
        integer, dimension(:,:), pointer :: neighbor => null()
        integer, dimension(:,:), pointer :: face => null()
        integer, dimension(:,:), pointer :: mpi_ibool => null()
        integer, dimension(:,:), allocatable :: smooth_A
        integer, dimension(:,:,:), pointer :: ibn => null()
        integer, dimension(:,:,:), pointer :: ibt => null()
        integer, dimension(:,:,:), pointer :: mpi_interface => null()
        integer, dimension(:,:,:), pointer :: mpi_connection => null()
        integer, dimension(:,:,:), pointer :: mpi_ibt => null()

        real(kind=CUSTOM_REAL) :: dtfactor
        real(kind=CUSTOM_REAL) :: pmlxmin
        real(kind=CUSTOM_REAL) :: pmlxmax
        real(kind=CUSTOM_REAL) :: pmlzmin
        real(kind=CUSTOM_REAL) :: pmlzmax
        real(kind=CUSTOM_REAL) :: minGLL

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
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vpu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: vsu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: rho => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: mu => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: lambda => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: qp => null()
        real(kind=CUSTOM_REAL), dimension(:), pointer :: qs => null()
        real(kind=custom_real), dimension(:), allocatable :: imp_vp, imp_vs
        real(kind=custom_real), dimension(:), allocatable :: vol
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: srcrecmask
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: sDT
        real(kind=custom_real), dimension(:,:), allocatable :: mpi_imp_vp, mpi_imp_vs
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: coord => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: matval => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: ylambda => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: ymu => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: wl => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: nx => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: nz => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: sj => null()
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: fscale => null()
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: A => null() !poro
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: B => null() !poro
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: E => null() !poro
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: AP => null() !poro
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: AM => null() !poro
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: BP => null() !poro
        real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: BM => null() !poro
    end type meshVar

    contains

    subroutine createRegularMesh(this, par, lsipar, errmsg)
        type(error_message) :: errmsg
        type(parameterVar) :: par
        type(lsi_parameter) :: lsipar
        type(lsiVar), dimension(:), allocatable :: lsi
        type(lsi_location), dimension(:), allocatable :: lsi_loc
        type (srcVar) :: src, tempsrc
        type (recVar) :: rec, temprec
        type(materialIndizes) :: matInd
        type(meshVar), dimension(:), allocatable :: db
        type(materialVar), dimension(:), allocatable :: mat
        type(porous_material), dimension(:), allocatable :: poromat
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
        real(kind=CUSTOM_REAL) :: minGLL, temp, temp_min, temp_max
        real(kind=CUSTOM_REAL), dimension(Np) :: x, z, r, s
        real(kind=CUSTOM_REAL), dimension(Np) :: rx, sx, rz, sz, jacobian
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: tempM1, tempM2
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: sDT, fDT
        integer, dimension(:), allocatable :: free_nodes
        integer, dimension(:), allocatable :: absorb_nodes
        real(kind=CUSTOM_REAL), dimension(2) :: center
        real(kind=CUSTOM_REAL) :: maskdist, longdimension

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
        real(kind=CUSTOM_REAL), dimension(:), pointer :: locrec_ang => null()

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
        !very important to choose the right epsi(lon), to small can cause problems
        real(kind=CUSTOM_REAL) :: epsi = 1e-3

        integer, dimension(:,:), allocatable :: lsi_matrix

        character(len=140) :: message

        this%poroelastic = par%poroelastic

        call addTrace(errmsg,myname)

        ! if we use MPI we will create databases
        if (par%nproc > 1) then
           allocate(db(par%nproc))
           db%poroelastic = this%poroelastic
        end if

        ! coordinates file
        if (par%log) write(*,'(a80)') "|                             read coordinate file                             |"
        filename=trim('mesh/coord')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open file!',myname, filename)
            call print(errmsg)
            stop
        endif
        read(19,*) this%ncoord
        allocate(this%coord(2,this%ncoord))

        ! calculate min and max of setup coordinates for use of PML
        this%pmlxmin=1.E10              ! define start values for size search
        this%pmlxmax=-1.E10
        this%pmlzmin=1.E10
        this%pmlzmax=-1.E10

        do i=1,this%ncoord
            read(19,*) dummy,this%coord(1,i),this%coord(2,i)
            if (this%coord(1,i) < this%pmlxmin) this%pmlxmin=this%coord(1,i)
            if (this%coord(1,i) > this%pmlxmax) this%pmlxmax=this%coord(1,i)
            if (this%coord(2,i) < this%pmlzmin) this%pmlzmin=this%coord(2,i)
            if (this%coord(2,i) > this%pmlzmax) this%pmlzmax=this%coord(2,i)
        end do
        close(19)

        !print size of setup
        if (par%log) write(*,'(a21,f8.2,a4,f8.2,a39)') '|  x axis ranges from ', this%pmlxmin, ' to ', this%pmlxmax, '                |'
        if (par%log) write(*,'(a21,f8.2,a4,f8.2,a39)') '|  z axis ranges from ', this%pmlzmin, ' to ', this%pmlzmax, '                |'

        ! meshfile
        if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
        if (par%log) write(*,'(a80)') "|                                read mesh file                                |"
        filename=trim('mesh/mesh')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open file!',myname, filename)
            call print(errmsg)
            stop
        endif

        read(19,*) this%nelem
        allocate(this%elem(3,this%nelem))
        allocate(this%elemNo(this%nelem))
        do i=1,this%nelem
            read(19,*) this%elemNo(i), this%elem(1,i),this%elem(2,i),this%elem(3,i)
        end do
        close(19)
        ! calculate nglob
        this%nglob = this%nelem*Np

        if (par%log) then
            write(*,'(a80)') "|------------------------------------------------------------------------------|"
            write(*,'(a40, i10, a30)') "|                Number of coordinates: ",this%ncoord , "                             |"
            write(*,'(a40, i10, a30)') "|                   Number of elements: ",this%nelem , "                             |"
            write(*,'(a40, i10, a30)') "|                     Number of points: ",this%nglob , "                             |"
            write(*,'(a80)') "|------------------------------------------------------------------------------|"
        end if
        ! allocate arrays
        call allocateMeshArrays(this,par)
        allocate(neighbortemp(3,this%nelem))
        allocate(sDT(this%nelem),fdt(this%nelem))

        !get neighbors
        call triangulation_neighbor_triangles(3,this%nelem,this%elem,neighbortemp)

        ! we have to reorder the neighboring
        this%neighbor(1,:) = neighbortemp(3,:)
        this%neighbor(2,:) = neighbortemp(1,:)
        this%neighbor(3,:) = neighbortemp(2,:)

        deallocate(neighbortemp)

        ! materials values
        if (par%poroelastic) then
            call setupMaterials(par, poromat, matInd, this%nelem, errmsg)
        else
            call setupMaterials(par, mat, matInd, this%nelem, iregion, errmsg)
        endif
        !Add some values to the mesh-type that have been read in at other places.
        this%nmat = par%matn
        this%mat = matInd%type

        ! free_nodes file
        filename=trim('mesh/free')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open file!' ,myname, filename)
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
        filename=trim('mesh/absorb')
        open(unit=19,file=trim(filename), status='old', iostat = ios)
        if (ios /= 0) then
            call add(errmsg,2,'could not open file!' ,myname, filename)
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

        if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
        if (par%log) write(*,'(a40, f12.8, a28)') "|           Minimum gll point distance: ", minGLL, "                           |"

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
        if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif

        temp_min = 1e30
        temp_max = 1e-24
        do i=1,this%nelem
            if ( sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(2,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(2,i)))**2) < temp_min) then
                temp_min = sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(2,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(2,i)))**2)
            elseif ( sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(3,i)))**2) < temp_min) then
                temp_min = sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(3,i)))**2)
            elseif ( sqrt( (this%coord(1,this%elem(2,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(2,i)) - this%coord(2,this%elem(3,i)))**2) < temp_min) then
                temp_min = sqrt( (this%coord(1,this%elem(2,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(2,i)) - this%coord(2,this%elem(3,i)))**2)
            endif
            if ( sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(2,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(2,i)))**2) > temp_max) then
                temp_max = sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(2,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(2,i)))**2)
            elseif ( sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(3,i)))**2) > temp_max) then
                temp_max = sqrt( (this%coord(1,this%elem(1,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(1,i)) - this%coord(2,this%elem(3,i)))**2)
            elseif ( sqrt( (this%coord(1,this%elem(2,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(2,i)) - this%coord(2,this%elem(3,i)))**2) > temp_max) then
                temp_max = sqrt( (this%coord(1,this%elem(2,i)) - this%coord(1,this%elem(3,i)))**2 + (this%coord(2,this%elem(2,i)) - this%coord(2,this%elem(3,i)))**2)
            endif
        end do

        if (par%log)  then
            write(*,'(a40, f12.8, a28)') "|      minimal edge length in the mesh: ", temp_min, "                           |"
            write(*,'(a40, f12.8, a28)') "|      maximal edge length in the mesh: ", temp_max, "                           |"
        end if
        ! LL LL test to determine minimal epsi(lon) for partitioning
        epsi = temp_min/1e4

        do ie=1,this%nelem
            this%vol(ie)=0.5*abs((this%coord(1,this%elem(2,ie))-this%coord(1,this%elem(1,ie)))*(this%coord(2,this%elem(3,ie))-this%coord(2,this%elem(1,ie)))-(this%coord(1,this%elem(3,ie))-this%coord(1,this%elem(1,ie)))*(this%coord(2,this%elem(2,ie))-this%coord(2,this%elem(1,ie))))
        enddo

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
            call geometricFactors2d(rx,sx,rz,sz,jacobian,this%vx(iv),this%vz(iv),this%Dr,this%Ds, errmsg)
            this%rx(iv)=rx
            this%sx(iv)=sx
            this%rz(iv)=rz
            this%sz(iv)=sz
            this%jacobian(iv)=jacobian
        end do

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

                        if ( (abs(this%vx(this%ibt(1,is,ie)) - this%vx(this%ibn(1,is,ie)) ) < epsi)&
                            .and. (abs(this%vz(this%ibt(1,is,ie)) - this%vz(this%ibn(1,is,ie)) ) < epsi)) then
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
        !!! LL LL open external model file EXPERIMENTAL
        if(par%extvel) then
            open(1,file = trim(par%externalfilename),status = 'old', iostat = ios)
            if (ios /= 0) then
                call add(errmsg,2,'could not open file!',myname, par%externalfilename)
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
                !The calculated coordinates are the centroid coordinates of the given triangle
                xs=(this%coord(1,this%elem(1,i))+this%coord(1,this%elem(2,i))+this%coord(1,this%elem(3,i)))/3
                zs=(this%coord(2,this%elem(1,i))+this%coord(2,this%elem(2,i))+this%coord(2,this%elem(3,i)))/3
                call define_external_model(xs,zs,iregion(this%mat(i)),this%rho(i),this%vp(i),this%vs(i),em_nx,em_nz,em_dx,em_dz,em_xmin,em_zmin,em_ipolg,em_xg,em_zg,em_rhog,em_vpg,em_vsg )
                this%mu(i) = this%vs(i)**2 * this%rho(i)
                this%lambda(i) = this%vp(i)**2 * this%rho(i)-2*this%mu(i)
            end do
        else
            this%nfluids = par%fluidn
            do i = 1, this%nelem
                if (this%poroelastic) then
                    this%A(:,:,i) = poromat(this%mat(i))%A
                    this%B(:,:,i) = poromat(this%mat(i))%B
                    this%E(:,:,i) = poromat(this%mat(i))%E
                    this%AP(:,:,i) = poromat(this%mat(i))%AP
                    this%AM(:,:,i) = poromat(this%mat(i))%AM
                    this%BP(:,:,i) = poromat(this%mat(i))%BP
                    this%BM(:,:,i) = poromat(this%mat(i))%BM
                    this%vp(i) = poromat(this%mat(i))%vmax !this is the maximum velocity, not vp
                    this%vs(i) = poromat(this%mat(i))%vmin !this is the minimum velocity, not vs
                    this%vpu(i) = this%vp(i) ! unrelaxed vp and vs is set to elastic values since
                    this%vsu(i) = this%vs(i) ! anelasticity and poroelasticity cannot used at the same time
                else
                    this%rho(i) = mat(this%mat(i))%rho
                    this%vp(i)  = mat(this%mat(i))%vp
                    this%vs(i)  = mat(this%mat(i))%vs
                    this%vpu(i) = mat(this%mat(i))%vpu
                    this%vsu(i) = mat(this%mat(i))%vsu

                    this%imp_vp(i) = mat(this%mat(i))%imp_vp
                    this%imp_vs(i) = mat(this%mat(i))%imp_vs

                    this%qp(i)        = mat(this%mat(i))%qp
                    this%qs(i)        = mat(this%mat(i))%qs
                    this%ylambda(:,i) = mat(this%mat(i))%ylambda(:)
                    this%ymu(:,i)     = mat(this%mat(i))%ymu(:)
                    this%wl(:,i)      = mat(this%mat(i))%omegaL(:)
                    this%pml(i)       = matInd%pml(i)
                    if (mat(this%mat(i))%vs .gt. 0) then
                        this%mu(i) = mat(this%mat(i))%mu
                    else
                        ! acoustic test
                        this%mu(i)=0
                    end if
                    this%lambda(i)=mat(this%mat(i))%lambda
                endif
            end do
        end if

        if (this%poroelastic) then
            write(filename,"('vmax_model.vtk')")
            call writeVtkTriMeshRealdata2D(trim(outpath)//filename,this%elem,this%coord,this%vp, 'vmax') !this is the maximum velocity, not vp
            write(filename,"('vmin_model.vtk')")
            call writeVtkTriMeshRealdata2D(trim(outpath)//filename,this%elem,this%coord,this%vs, 'vmin') !this is the minimum velocity, not vs
        else
            write(filename,"('vp_model.vtk')")
            call writeVtkTriMeshRealdata2D(trim(outpath)//filename,this%elem,this%coord,this%vp, 'vp')
            write(filename,"('vs_model.vtk')")
            call writeVtkTriMeshRealdata2D(trim(outpath)//filename,this%elem,this%coord,this%vs, 'vs')
            write(filename,"('rho_model.vtk')")
            call writeVtkTriMeshRealdata2D(trim(outpath)//filename,this%elem,this%coord,this%rho, 'rho')
        endif

        if (par%set_pml) then
            write(filename,"('pml_model.vtk')")
            call writeVtkTriMeshRealdata2D(trim(outpath)//filename,this%elem,this%coord,real(this%pml, CUSTOM_REAL), 'pml')
        endif

        ! compute dtfactor for every element
        temp=1e7
        do i=1,this%nelem
            fDT(i)= (2.0/3.0)*minGLL*(sDT(i)/this%vpu(i))
            if (fDT(i) < temp) then
                temp=fDT(i)
            end if
        end do
        this%dtfactor=temp
        if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
        if (par%log) write(*,'(a40, f8.2, a32)') "|                                 vmax: ",maxval(this%vp), "                               |"
        if (par%log) write(*,'(a40, f8.2, a32)') "|                                 vmin: ",minval(this%vs), "                               |"
        if (par%attenuation) then
            if (par%log) write(*,'(a40, f8.2, a32)') "|                       unrelaxed vmax: ",maxval(this%vp), "                               |"
            if (par%log) write(*,'(a40, f8.2, a32)') "|                       unrelaxed vmin: ",minval(this%vs), "                               |"
        endif
        if (par%log) write(*,'(a40, es12.5, a28)') "|recommended time step with CFL factor: ",this%dtfactor*par%cfl, "                           |"

        if (lsipar%lsi) then
            if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
            if (par%log) write(*,'(a80)') "|                         creating fracture(s)...                              |"
            call readLSIlocation(lsi_loc, this%coord, par%set_pml, par%pml_delta, 19, "data/fracs", errmsg)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
            call createFracture(this, lsi, lsi_loc, lsipar)
            call lsiMatrix(lsi, lsipar, this, lsi_matrix)
        else
            lsipar%lsin = 0
        end if

        if (.not. par%autodt .and. (par%dt > this%dtfactor*par%cfl)) then
            write(message,'(a18, es13.6, a44, es13.6, a52)') "The input for dt (", par%dt,") is larger then the recommended time step (", this%dtfactor*par%cfl,") for this simulation. This may cause instabilities."
            call add(errmsg, 1, message, myname)
        end if

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

            db(:)%pmlxmin=this%pmlxmin
            db(:)%pmlxmax=this%pmlxmax
            db(:)%pmlzmin=this%pmlzmin
            db(:)%pmlzmax=this%pmlzmax
            db(:)%minGLL = minGLL

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
                allocate(db(iproc)%sDT(db(iproc)%nelem))
                if (this%poroelastic) then
                    allocate(db(iproc)%A(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                             db(iproc)%B(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                             db(iproc)%E(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                             db(iproc)%AM(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                             db(iproc)%AP(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                             db(iproc)%BM(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                             db(iproc)%BP(5+3*par%fluidn,5+3*par%fluidn,this%nelem))
                    allocate(db(iproc)%vp(this%nelem),db(iproc)%vs(db(iproc)%nelem),db(iproc)%vpu(this%nelem),db(iproc)%vsu(db(iproc)%nelem))
                else
                    allocate(db(iproc)%vp(db(iproc)%nelem),db(iproc)%vs(db(iproc)%nelem),db(iproc)%rho(db(iproc)%nelem),&
                             db(iproc)%qp(db(iproc)%nelem),db(iproc)%qs(db(iproc)%nelem),&
                             db(iproc)%mu(db(iproc)%nelem),db(iproc)%lambda(db(iproc)%nelem))
                    allocate(db(iproc)%vpu(db(iproc)%nelem),db(iproc)%vsu(db(iproc)%nelem))
                endif
                allocate(db(iproc)%pml(db(iproc)%nelem))
                allocate(db(iproc)%ylambda(nMB,db(iproc)%nelem),db(iproc)%ymu(nMB,db(iproc)%nelem),db(iproc)%wl(nMB,db(iproc)%nelem))
                allocate(db(iproc)%ibt(NGLL,3,db(iproc)%nelem),db(iproc)%ibn(NGLL,3,db(iproc)%nelem))
                allocate(db(iproc)%mpi_ibt(NGLL,3,db(iproc)%nelem))
                allocate(db(iproc)%loc2glob_nodes(db(iproc)%ncoord))

                allocate(db(iproc)%elemNo(db(iproc)%nelem))
                allocate(db(iproc)%imp_vp(db(iproc)%nelem))
                allocate(db(iproc)%imp_vs(db(iproc)%nelem))
                allocate(db(iproc)%vol(db(iproc)%nelem))
                allocate(db(iproc)%srcrecmask(db(iproc)%nelem))
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

                ! build local elements / bulid local node indizes
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
                    call geometricFactors2d(rx,sx,rz,sz,jacobian,db(iproc)%vx(iv),db(iproc)%vz(iv),db(iproc)%Dr,db(iproc)%Ds, errmsg)
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
                        db(iproc)%sDT(l)=sDT(ie)
                        db(iproc)%mat(l)=this%mat(ie)
                        if (this%poroelastic) then
                            db(iproc)%nfluids      = this%nfluids
                            db(iproc)%A(:,:,l)     = this%A(:,:,ie)
                            db(iproc)%B(:,:,l)     = this%B(:,:,ie)
                            db(iproc)%E(:,:,l)     = this%E(:,:,ie)
                            db(iproc)%AM(:,:,l)    = this%AM(:,:,ie)
                            db(iproc)%AP(:,:,l)    = this%AP(:,:,ie)
                            db(iproc)%BM(:,:,l)    = this%BM(:,:,ie)
                            db(iproc)%BP(:,:,l)    = this%BP(:,:,ie)
                            db(iproc)%vp(l)        = this%vp(ie)
                            db(iproc)%vs(l)        = this%vs(ie)
                        else
                            db(iproc)%face(:,l)    = this%face(:,ie)
                            db(iproc)%mat(l)       = this%mat(ie)
                            db(iproc)%vp(l)        = this%vp(ie)
                            db(iproc)%vs(l)        = this%vs(ie)
                            db(iproc)%vpu(l)       = this%vpu(ie)
                            db(iproc)%vsu(l)       = this%vsu(ie)
                            db(iproc)%rho(l)       = this%rho(ie)
                            db(iproc)%qp(l)        = this%qp(ie)
                            db(iproc)%qs(l)        = this%qs(ie)

                            db(iproc)%imp_vp(l)    = this%imp_vp(ie)
                            db(iproc)%imp_vs(l)    = this%imp_vs(ie)

                            db(iproc)%ylambda(:,l) = this%ylambda(:,ie)
                            db(iproc)%ymu(:,l)     = this%ymu(:,ie)
                            db(iproc)%wl(:,l)      = this%wl(:,ie)
                            db(iproc)%lambda(l)    = this%lambda(ie)
                            db(iproc)%mu(l)        = this%mu(ie)

                            db(iproc)%elemNo(l)    = this%elemNo(ie)
                        endif
                        db(iproc)%pml(l)       = this%pml(ie)
                        db(iproc)%vol(l)       = this%vol(ie)
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
                            if ( .not. ( (abs(db(iproc)%vx(db(iproc)%ibt(1,is,ie)) - db(iproc)%vx(db(iproc)%ibn(1,is,ie)) ) < epsi)&
                                .and. (abs(db(iproc)%vz(db(iproc)%ibt(1,is,ie)) - db(iproc)%vz(db(iproc)%ibn(1,is,ie)) ) < epsi)) ) then
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

            !create Impedance mpi connection array
            do iproc = 1, par%nproc
                allocate(db(iproc)%mpi_imp_vp(db(iproc)%pinterfaces, 3))
                allocate(db(iproc)%mpi_imp_vs(db(iproc)%pinterfaces, 3))
                i = 1
                do ie = 1, db(iproc)%nelem
                    do is = 1, nsurface
                        if (db(iproc)%neighbor(is,ie) == -1) then
                            db(iproc)%mpi_imp_vp(i, 1) = db(db(iproc)%mpi_interface(1,is,ie))%imp_vp(db(iproc)%mpi_interface(2,is,ie))
                            db(iproc)%mpi_imp_vp(i, 2) = db(iproc)%mpi_interface(1,is,ie)
                            db(iproc)%mpi_imp_vp(i, 3) = db(iproc)%mpi_interface(2,is,ie)
                            db(iproc)%mpi_imp_vs(i, 1) = db(db(iproc)%mpi_interface(1,is,ie))%imp_vs(db(iproc)%mpi_interface(2,is,ie))
                            db(iproc)%mpi_imp_vs(i, 2) = db(iproc)%mpi_interface(1,is,ie)
                            db(iproc)%mpi_imp_vs(i, 3) = db(iproc)%mpi_interface(2,is,ie)
                            i = i+1
                        end if
                    end do
                end do
            end do

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
                                if ((abs(db(iproc)%vx(mpi_ti) - db(l)%vx(mpi_ni) ) < epsi)&
                                .and. (abs(db(iproc)%vz(mpi_ti) - db(l)%vz(mpi_ni))  < epsi)) then
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
            ! find source in the global mesh
            !"--------------------------------------------------------------------------------------"
            call initSource(par,src,this%coord, errmsg)

            !find sources in global mesh
            call findSource(par,src,this%vx,this%vz,this%nelem,this%ibool,this%coord,this%elem, this%Dr, this%Ds, this%pml, 'glob', errmsg)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
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
                    call prepareRecalcSrc(tempsrc,locnsrc,locsrcxz,locsrctype,locsrcstf,locsrcextwavelet,locsrcf0,locsrcfactor,locsrcangle_force,locsrcM,locdelay,locsrcnr)
                    call findSource(par,tempsrc,db(iproc)%vx,db(iproc)%vz,db(iproc)%nelem,db(iproc)%ibool, &
                                    db(iproc)%coord,db(iproc)%elem,db(iproc)%dr,db(iproc)%ds, db(iproc)%pml, 'locl', errmsg)
                    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                    write(filename,"('srcVar',i6.6)") iproc
                    call writeSrcVar(tempsrc,trim(outpath)//filename)
                    deallocate(locsrcnr)
                    deallocate(locsrcxz)
                    deallocate(locsrctype,locsrcstf,locsrcextwavelet,locsrcf0,locsrcfactor,locsrcangle_force,locsrcM)
                    call deallocSrcVar(tempsrc)
                end if
            end do

            deallocate(tempv)
            deallocate(srcnr)

            !"--------------------------------------------------------------------------------------"
            ! find receiver in the global mesh
            !"--------------------------------------------------------------------------------------"
            call initReceiver(par,rec, this%coord, errmsg)

            call findReceiver(par,rec,this%vx,this%vz,this%nelem,this%ibool,this%coord,this%elem, this%Dr, this%Ds, this%pml, 'glob', errmsg)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
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
                    allocate(locrec_ang(tempv(iproc)))
                    allocate(locrecnr(tempv(iproc)))
                    allocate(locrecnum(tempv(iproc)))
                    locnrec=tempv(iproc)
                    do i=1,tempv(iproc)
                        locrecnum(i) = recnum(iproc,i)
                        locrecxz(:,i) = rec%recxz(:,locrecnum(i))
                        locrecnr(i) = rec%recnr(locrecnum(i))
                        locrec_ang(i) = rec%rec_ang(locrecnum(i))
                        j=j+1
                    end do
                    call prepareRecalcRec(temprec,locnrec,locrecxz,locrecnr,locrec_ang)
                    call findReceiver(par,temprec,db(iproc)%vx,db(iproc)%vz,db(iproc)%nelem,db(iproc)%ibool,&
                                      db(iproc)%coord,db(iproc)%elem,db(iproc)%dr,db(iproc)%ds, db(iproc)%pml, 'locl',errmsg)
                    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
                    write(filename,"('recVar',i6.6)") iproc
                    call writeRecVar(temprec,trim(outpath)//filename)
                    deallocate(locrecnum)
                    deallocate(locrecxz)
                    deallocate(locrec_ang)
                    deallocate(locrecnr)
                    call deallocRecVar(temprec)
                end if
            end do

            deallocate(tempv)
            allocate(tempvr(size(src%srcxz(1,:))))
            tempvr=0.0
            !write src to vtk
            write(filename,"('src.vtk')")
            call writeVtkNodes(trim(outpath)//filename, src%srcxz(1,:),tempvr,src%srcxz(2,:))
            deallocate(tempvr)

            allocate(tempvr(size(rec%recxz(1,:))))
            tempvr=0.0
            ! write receiver to vtk
            write(filename,"('rec.vtk')")
            call writeVtkNodes(trim(outpath)//filename, rec%recxz(1,:),tempvr,rec%recxz(2,:))
            deallocate(tempvr)

            !LSI only
            if (lsipar%lsi) then
                if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
                if (par%log) write(*,'(a80)') "|                          writing lsi database file                           |"
                if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
                call writeLSIDatabase(lsi_matrix, this, db, par%nproc, part, glob2loc_elmnts, errmsg)
            end if

            ! Mask for adjoint inversion
            if (this%pmlxmax-this%pmlxmin >= this%pmlzmax-this%pmlzmin) then
                longdimension=this%pmlxmax-this%pmlxmin
            else
                longdimension=this%pmlzmax-this%pmlzmin
            endif

            maskdist=par%maskradius*longdimension
            this%srcrecmask=1.

            if (par%inversion) then
                do i=1,size(src%srcelem)
                    do ie=1,this%nelem
                        center=(this%coord(:,this%elem(1,ie))+this%coord(:,this%elem(2,ie))+this%coord(:,this%elem(3,ie)))/3.
                        if ((src%srcxz(1,i)-center(1))**2+(src%srcxz(2,i)-center(2))**2 < (maskdist/2.)**2) then
                            this%srcrecmask(ie)=0.
                        else if ((src%srcxz(1,i)-center(1))**2+(src%srcxz(2,i)-center(2))**2 < maskdist**2) then
                            this%srcrecmask(ie)=this%srcrecmask(ie)*0.5*(1.+(cos(2*pi/maskdist*sqrt((src%srcxz(1,i)-center(1))**2+(src%srcxz(2,i)-center(2))**2))))
                        endif
                    enddo
                enddo
                do i=1,size(rec%recelem)
                    do ie=1,this%nelem
                        center=(this%coord(:,this%elem(1,ie))+this%coord(:,this%elem(2,ie))+this%coord(:,this%elem(3,ie)))/3.
                        if ((rec%recxz(1,i)-center(1))**2+(rec%recxz(2,i)-center(2))**2 < (maskdist/2.)**2) then
                            this%srcrecmask(ie)=0.
                        else if ((rec%recxz(1,i)-center(1))**2+(rec%recxz(2,i)-center(2))**2 < maskdist**2) then
                            this%srcrecmask(ie)=this%srcrecmask(ie)*0.5*(1.+(cos(2*pi/maskdist*sqrt((rec%recxz(1,i)-center(1))**2+(rec%recxz(2,i)-center(2))**2))))
                        endif
                    enddo
                enddo

                do ie=1,this%nelem
                    if (this%pml(ie)>0.and.par%set_pml) then
                        this%srcrecmask(ie)=0.
                    endif
                enddo

                do ie=1,this%nelem
                    center=(this%coord(:,this%elem(1,ie))+this%coord(:,this%elem(2,ie))+this%coord(:,this%elem(3,ie)))/3.
                    do l=1,nfree
                        if ((this%coord(1,int(free_nodes(l)))-center(1))**2+(this%coord(2,int(free_nodes(l)))-center(2))**2 < (maskdist/2.)**2) then
                            this%srcrecmask(ie)=0.
                        endif
                    enddo
                enddo

                do iproc=1,par%nproc
                    do ie=1, this%nelem
                        if (part(ie)== iproc) then
                            l = glob2loc_elmnts(ie)
                            db(iproc)%srcrecmask(l) = this%srcrecmask(ie)
                        endif
                    enddo
                enddo

                write(filename,"('/mask.vtk')")
                filename=trim(outpath)//trim(filename)
                call writeVtkTriMeshRealdata2D(filename, this%elem, this%coord, this%srcrecmask,'Mask')
            endif
            ! build smoothing matrix

            do iproc=1,par%nproc
                allocate(db(iproc)%smooth_A(db(iproc)%nelem,4))
                if (par%inversion) then
                    do ie=1,db(iproc)%nelem
                        i=0
                        do j=1,3
                            if (db(iproc)%neighbor(j,ie) /= 0) i= i + 1
                        enddo
                        db(iproc)%smooth_A(ie,1) = - i
                        do j=1,3
                            if (db(iproc)%neighbor(j,ie)==0) then
                                db(iproc)%smooth_A(ie,j+1)=0
                            else if (db(iproc)%neighbor(j,ie)>0) then
                                db(iproc)%smooth_A(ie,j+1)=db(iproc)%neighbor(j,ie)
                            else if (db(iproc)%neighbor(j,ie)==-1) then
                                db(iproc)%smooth_A(ie,j+1)=-db(iproc)%mpi_ibool(j,ie) ! negative sign to mark MPI boundary
                            endif
                        enddo
                    enddo
                endif
            enddo


             ! write DATABASE
            do iproc=1,par%nproc
                write(filename,"('meshVar',i6.6)") iproc
                call writeMeshVar(db(iproc),trim(outpath)//filename)
                call deallocMeshVar(db(iproc))
            end do !nproc create databases
            deallocate(absorb_nodes)
            deallocate(free_nodes)
        end if ! par%nproc >1

        !"--------------------------------------------------------------------------------------"
        ! deallocate arrays
        !"--------------------------------------------------------------------------------------"
        if (par%nproc > 1) then
            do iproc =1, par%nproc
                call deallocMeshVar(db(iproc))
            enddo
        end if
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
        call deallocSrcVar(src)
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
        if (associated(this%A)) deallocate(this%A)
        if (associated(this%B)) deallocate(this%B)
        if (associated(this%E)) deallocate(this%E)
        if (associated(this%AP)) deallocate(this%AP)
        if (associated(this%AM)) deallocate(this%AM)
        if (associated(this%BP)) deallocate(this%BP)
        if (associated(this%BM)) deallocate(this%BM)
        if (associated(this%vp)) deallocate(this%vp)
        if (associated(this%vs)) deallocate(this%vs)
        if (associated(this%vp)) deallocate(this%vpu)
        if (associated(this%vs)) deallocate(this%vsu)
        if (associated(this%rho)) deallocate(this%rho)
        if (associated(this%qp)) deallocate(this%qp)
        if (associated(this%qs)) deallocate(this%qs)
        if (associated(this%mu)) deallocate(this%mu)
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
        if (allocated(this%elemNo)) deallocate(this%elemNo)
        if (allocated(this%imp_vp)) deallocate(this%imp_vp)
        if (allocated(this%imp_vs)) deallocate(this%imp_vs)
        if (allocated(this%mpi_imp_vp)) deallocate(this%mpi_imp_vp)
        if (allocated(this%mpi_imp_vs)) deallocate(this%mpi_imp_vs)
        if (allocated(this%vol)) deallocate(this%vol)
        if (allocated(this%srcrecmask)) deallocate(this%srcrecmask)
        if (allocated(this%smooth_A)) deallocate(this%smooth_A)
        if (allocated(this%sDT)) deallocate(this%sDT)
    end subroutine deallocMeshvar

    subroutine allocateMeshArrays(this,par)
        type(meshVar) :: this
        type(parameterVar) :: par
        !"--------------------------------------------------------------------------------------"
        allocate(this%vx(this%nglob),this%vz(this%nglob))
        allocate(this%rx(this%nglob), this%rz(this%nglob), this%sx(this%nglob), this%sz(this%nglob))
        allocate(this%jacobian(this%nglob))
        allocate(this%mat(this%nelem))
        allocate(this%ibool(Np,this%nelem))
        allocate(this%nx(3*NGLL,this%nelem), this%nz(3*NGLL,this%nelem), this%sj(3*NGLL,this%nelem), this%fscale(3*NGLL,this%nelem))
        allocate(this%neighbor(3,this%nelem))
        allocate(this%face(3,this%nelem))
        allocate(this%ibt(NGLL,3,this%nelem),this%ibn(NGLL,3,this%nelem))
        if (this%poroelastic) then
            allocate(this%A(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                     this%B(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                     this%E(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                     this%AM(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                     this%AP(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                     this%BM(5+3*par%fluidn,5+3*par%fluidn,this%nelem),&
                     this%BP(5+3*par%fluidn,5+3*par%fluidn,this%nelem))
            allocate(this%vp(this%nelem),this%vs(this%nelem), this%vpu(this%nelem),this%vsu(this%nelem))
        else
            allocate(this%vp(this%nelem),this%vs(this%nelem),this%rho(this%nelem),this%qp(this%nelem),this%qs(this%nelem),this%mu(this%nelem),this%lambda(this%nelem))
            allocate(this%vpu(this%nelem),this%vsu(this%nelem))
        endif
        allocate(this%pml(this%nelem))
        allocate(this%ylambda(nMB,this%nelem),this%ymu(nMB,this%nelem), this%wl(nMB,this%nelem))
        allocate(this%imp_vp(this%nelem), this%imp_vs(this%nelem))
        allocate(this%vol(this%nelem))
        allocate(this%srcrecmask(this%nelem))
        allocate(this%smooth_A(this%nelem,4))
        allocate(this%sDT(this%nelem))
    end subroutine allocateMeshArrays

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
        real(kind=CUSTOM_REAL) :: dtfactor, minGLL
        integer :: nelx,nelz
        integer :: nfluids

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
        real(kind=CUSTOM_REAL), allocatable, dimension(:) :: vp,vs,vpu,vsu,rho,mu,lambda,qp,qs
        logical :: poroelastic
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:) :: A,B,E,AP,AM,BP,BM
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
        minGLL = this%minGLL

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
        if (this%poroelastic) then
            allocate(A(size(this%A(:,1,1)),size(this%A(1,:,1)),nelem),&
                     B(size(this%B(:,1,1)),size(this%B(1,:,1)),nelem),&
                     E(size(this%E(:,1,1)),size(this%E(1,:,1)),nelem),&
                     AM(size(this%AM(:,1,1)),size(this%AM(1,:,1)),nelem),&
                     AP(size(this%AP(:,1,1)),size(this%AP(1,:,1)),nelem),&
                     BM(size(this%BM(:,1,1)),size(this%BM(1,:,1)),nelem),&
                     BP(size(this%BP(:,1,1)),size(this%BP(1,:,1)),nelem))
            allocate(vp(nelem),vs(nelem), vpu(nelem),vsu(nelem))
        else
            allocate(vp(nelem),vs(nelem),rho(nelem),qp(nelem),qs(nelem),mu(nelem),lambda(nelem))
            allocate(vpu(nelem),vsu(nelem))
        endif
        allocate(pml(nelem))
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
        poroelastic = this%poroelastic
        if (this%poroelastic) then
            nfluids = this%nfluids
            A = this%A
            B = this%B
            E = this%E
            AM = this%AM
            AP = this%AP
            BM = this%BM
            BP = this%BP
            vp = this%vp
            vs = this%vs
            vpu= this%vpu
            vsu= this%vsu
        else
            vp = this%vp
            vs = this%vs
            vpu= this%vpu
            vsu= this%vsu
            rho = this%rho
            qp = this%qp
            qs = this%qs
            mu = this%mu
            lambda = this%lambda
        endif
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

        !write(*,'(a31, a18, a31)') "|                         write", trim(filename), "                              |"

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
        write(27) minGLL

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
        write(27) poroelastic
        if (this%poroelastic) then
            write(27) nfluids
            write(27) A
            write(27) B
            write(27) E
            write(27) AM
            write(27) AP
            write(27) BM
            write(27) BP
            write(27) vp
            write(27) vs
            write(27) vpu
            write(27) vsu
        else
            write(27) vp
            write(27) vs
            write(27) vpu
            write(27) vsu
            write(27) rho
            write(27) qp
            write(27) qs
            write(27) mu
            write(27) lambda
        endif
        write(27) ymu
        write(27) ylambda
        write(27) wl
        write(27) pml
        write(27) this%pmlxmin
        write(27) this%pmlxmax
        write(27) this%pmlzmin
        write(27) this%pmlzmax
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

        write(27) this%elemNo
        write(27) this%imp_vp
        write(27) this%imp_vs
        write(27) this%mpi_imp_vp
        write(27) this%mpi_imp_vs
        write(27) this%vol
        write(27) this%srcrecmask
        write(27) this%smooth_A
        write(27) this%sDT

        close(27)

        !write(*,*) "done writing meshVar"

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
        if (this%poroelastic) then
            deallocate(A,B,E,AM,AP,BM,BP,vp,vs,vpu,vsu)
        else
            deallocate(vp,vs,vpu,vsu,rho,qp,qs,mu,lambda)
        endif
        deallocate(pml,ymu,ylambda,wl)
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

        !write(*,'(a31, a18, a31)') "|                       reading", trim(filename), "                              |"
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
        read(27) this%minGLL

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
        allocate(this%vp(this%nelem),this%vs(this%nelem),this%rho(this%nelem),this%qp(this%nelem),this%qs(this%nelem),this%mu(this%nelem),this%lambda(this%nelem),this%pml(this%nelem))
        allocate(this%vpu(this%nelem),this%vsu(this%nelem))
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

        allocate(this%elemNo(this%nelem))
        allocate(this%imp_vp(this%nelem))
        allocate(this%imp_vs(this%nelem))
        allocate(this%mpi_imp_vp(this%pinterfaces, 3))
        allocate(this%mpi_imp_vs(this%pinterfaces, 3))
        allocate(this%vol(this%nelem))
        allocate(this%srcrecmask(this%nelem))
        allocate(this%smooth_A(this%nelem,4))
        allocate(this%sDT(this%nelem))

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
        read(27) this%poroelastic
        if (this%poroelastic) then
            read(27) this%nfluids
            allocate(this%A(5+3*this%nfluids,5+3*this%nfluids,this%nelem),&
                     this%B(5+3*this%nfluids,5+3*this%nfluids,this%nelem),&
                     this%E(5+3*this%nfluids,5+3*this%nfluids,this%nelem),&
                     this%AM(5+3*this%nfluids,5+3*this%nfluids,this%nelem),&
                     this%AP(5+3*this%nfluids,5+3*this%nfluids,this%nelem),&
                     this%BM(5+3*this%nfluids,5+3*this%nfluids,this%nelem),&
                     this%BP(5+3*this%nfluids,5+3*this%nfluids,this%nelem))
            read(27) this%A
            read(27) this%B
            read(27) this%E
            read(27) this%AM
            read(27) this%AP
            read(27) this%BM
            read(27) this%BP
            read(27) this%vp
            read(27) this%vs
            read(27) this%vpu
            read(27) this%vsu
        else
            allocate(this%rho(this%nelem),this%qp(this%nelem),this%qs(this%nelem),this%mu(this%nelem),this%lambda(this%nelem))
            read(27) this%vp
            read(27) this%vs
            read(27) this%vpu
            read(27) this%vsu
            read(27) this%rho
            read(27) this%qp
            read(27) this%qs
            read(27) this%mu
            read(27) this%lambda
        endif
        read(27) this%ymu
        read(27) this%ylambda
        read(27) this%wl
        read(27) this%pml
        read(27) this%pmlxmin
        read(27) this%pmlxmax
        read(27) this%pmlzmin
        read(27) this%pmlzmax
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

        read(27) this%elemNo
        read(27) this%imp_vp
        read(27) this%imp_vs
        read(27) this%mpi_imp_vp
        read(27) this%mpi_imp_vs
        read(27) this%vol
        read(27) this%srcrecmask
        read(27) this%smooth_A
        read(27) this%sDT

        close(27)
    end subroutine readMeshVar

    subroutine writeLSIDatabase(lsi_matrix, mesh, db, nproc, part, glob2loc_elements, errmsg)
        !input
        type(meshVar) :: mesh
        type(meshVar), dimension(:) :: db
        type(error_message) :: errmsg
        integer, intent(in) :: nproc
        integer, dimension(:), intent(in) :: part               !Vector that maps the processor number to the element number size(part) = number of elements
        integer, dimension(:), intent(in) :: glob2loc_elements
        integer, dimension(:,:), intent(in) :: lsi_matrix
        !local
        character(len=16) :: myname = 'writeLSIDatabase'
        character(len=22) :: outfile
        integer :: iproc, ilsi, element, surface, lsi_in_db !counter
        integer :: lu
        integer :: lsi_proc
        integer :: lsi_local
        integer, dimension(3) :: lsiprop
        logical :: isLSI
        logical, dimension(3) :: lsisurface


        call addTrace(errmsg, myname)

        lu = 2

        if (nproc > 1) then
            do iproc = 1, nproc
                lsi_in_db = 0
                !count how many lsi-elements are in the specific partition. This is later needed to allocate the
                !correct amount of interfaces in the solver.
                do ilsi = 1, size(lsi_matrix(:,1))
                    lsi_proc = part(lsi_matrix(ilsi, 1))
                    if (lsi_proc == iproc) then
                        lsi_in_db = lsi_in_db + 1
                    end if
                end do

                write(outfile, "('elementToLSI', i6.6)") iproc
                open(lu, file = trim(outpath)//outfile, status = "unknown")
                write(*,'(a29, a22, a29)') "|                      write ", trim(outfile), "                            |"
                write(2,*) lsi_in_db

                do element = 1, db(iproc)%nelem
                    isLSI = .false.
                    lsisurface = .false.
                    lsiprop = 0
                    do ilsi = 1, size(lsi_matrix(:,1))
                        lsi_proc = part(lsi_matrix(ilsi, 1))
                        if (lsi_proc == iproc) then
                            lsi_local = glob2loc_elements(lsi_matrix(ilsi,1))   !Local numbering of the lsi-elements
                            if (element == lsi_local) then
                                isLSI = .true.
                                do surface = 1, nsurface
                                    if (lsi_matrix(ilsi, surface + 1) > 0) then
                                        lsisurface(surface) = .true.
                                        lsiprop(surface) = lsi_matrix(ilsi, surface + 4)
                                    end if
                                end do
                            end if
                        end if
                    end do
                    write(lu,*) isLSI, lsisurface(1), lsisurface(2), lsisurface(3), lsiprop(1), lsiprop(2), lsiprop(3)
                end do
                close(lu)
            enddo
        else
            open (lu, file = trim(outpath)//"elementToLSI", status = "unknown")
            do element = 1, mesh%nelem
                isLSI = .false.
                lsisurface = .false.
                lsiprop = 0
                do ilsi = 1, size(lsi_matrix(:,1))
                    if (element == lsi_matrix(ilsi,1)) then
                        isLSI = .true.
                        do surface = 1, nsurface
                            if (lsi_matrix(ilsi, surface + 1) > 0) then
                                lsisurface(surface) = .true.
                                lsiprop(surface) = lsi_matrix(ilsi, surface + 4)
                            end if
                        end do
                    end if
                end do
                write(lu,*) isLSI, lsisurface(1), lsisurface(2), lsisurface(3), lsiprop(1), lsiprop(2), lsiprop(3)
            end do
            close(lu)
        end if

    end subroutine writeLSIDatabase

    subroutine lsiMatrix(lsi, lsipar, mesh, lsi_matrix)
        !input
        type(meshVar) :: mesh
        type(lsiVar), dimension(:) :: lsi
        type(lsi_parameter) :: lsipar
        !output
        integer, dimension(:,:), allocatable, intent(out) :: lsi_matrix !the first entry in every row is the element.
                                                                        !the next three are the neighbours that connect over an lsi.
                                                                        !if no such neighbour is found the entry is 0
                                                                        !the last three contain the propertyindex for the surface
        !local
        integer :: i, element, surface, ilsi, cnt
        integer :: neighbour
        integer :: lsinumber
        integer, dimension(7) :: tmpline = 0
        integer, dimension(:,:), allocatable :: tmpmatrix

        allocate(tmpmatrix(2*lsipar%lsin,7)) !Allocate temporary matrix with the maximum amount of lsi-elements possible.
                                             !(this is always 2 times the amount of the interfaces as each interface has two sides)

        tmpmatrix = 0
        lsinumber = 0
        cnt = 1
        do element = 1, mesh%nelem
            do surface = 1, nsurface
                neighbour = mesh%neighbor(surface, element)
                do ilsi = 1, lsipar%lsin
                    if ((element == lsi(ilsi)%elements(1) .and. neighbour == lsi(ilsi)%elements(2)) &
                    .or. (element == lsi(ilsi)%elements(2) .and. neighbour == lsi(ilsi)%elements(1))) then
                        if (tmpline(1) == element) then
                            tmpline(surface+1) = neighbour
                            tmpline(surface+4) = lsi(ilsi)%prop
                        else
                            tmpline(1) = element
                            tmpline(surface+1) = neighbour
                            tmpline(surface+4) = lsi(ilsi)%prop
                        endif
                    end if
                end do
            end do
            if (tmpline(1) > 0) then
                lsinumber = lsinumber + 1
                do i = 1, 7
                    tmpmatrix(cnt, i) = tmpline(i)
                enddo
                cnt = cnt +1
                tmpline = 0
            end if
        end do

        allocate(lsi_matrix(lsinumber, 7))   !Allocate the lsi element-matrix to the total number of elements that have lsi!
        lsi_matrix = tmpmatrix(:lsinumber,:)
        deallocate(tmpmatrix)
    end subroutine lsiMatrix

    subroutine createFracture(mesh, lsi, lsi_loc, lsipar)
        !This suborutine created the fracture along the element faces from a defined starting point to an end point.
        !The node (Vertex) clostest to the input-value will be used as a starting end ending location.
        !input
        type(meshVar) :: mesh
        type(lsi_parameter) :: lsipar
        type(lsiVar), dimension(:), allocatable :: lsi
        type(lsi_location), dimension(:) :: lsi_loc
        !output
        !local
        integer :: i, frac, element, node
        integer :: lsi_coord_number, frac_number, lsi_number
        integer :: start_node
        integer :: start_element
        integer :: start_point, end_point
        integer, dimension(:,:), allocatable :: frac_points
        integer, dimension(:,:), allocatable :: tmp_lsi_elements
        integer, dimension(:), allocatable :: tmp_property
        real(kind=custom_real) :: node_x, node_z
        real(kind=custom_real) :: mindiff_start, diff_start
        real(kind=custom_real) :: mindiff_end, diff_end
        real(kind=custom_real), dimension(:,:), allocatable :: lsi_coords

        !Allocate enough memeory to save the lsi coordinates. The absolute max is the total ammount of nodes (coordinates)

        allocate(lsi_coords(mesh%ncoord, 2))
        allocate(tmp_lsi_elements(mesh%ncoord, 2))
        allocate(tmp_property(mesh%ncoord))
        allocate(frac_points(size(lsi_loc), 3))

        !Placing thes variables above the fracture loop ensures that these are counted continuously over multiple fractures
        lsi_number = 0
        lsi_coord_number = 0

        do frac = 1, size(lsi_loc)
            !frac_number counts the number of slip interfaces that make up the current farcture
            frac_number = 1
            !For now this value is hard coded. This is just a starting point which is to be minimized.
            !One option as a start value might be the maximum length of an element face.
            mindiff_start = 20.
            mindiff_end = 20.
            !Loop over the complete mesh in order to find the start and end point of the fracture in the mesh.
            do element = 1, mesh%nelem
                do node = 1, 3 !The trinangle has three nodes... ;)
                    node_x = mesh%coord(1,mesh%elem(node, element))
                    node_z = mesh%coord(2,mesh%elem(node, element))
                    diff_start = sqrt((node_x - lsi_loc(frac)%start_x)**2 + (node_z - lsi_loc(frac)%start_z)**2)
                    diff_end   = sqrt((node_x - lsi_loc(frac)%end_x)**2 + (node_z - lsi_loc(frac)%end_z)**2)
                    if (diff_start < mindiff_start) then
                        mindiff_start = diff_start
                        start_point = mesh%elem(node, element)
                        start_element = element
                        start_node = node
                    endif
                    if (diff_end < mindiff_end) then
                        mindiff_end = diff_end
                        end_point = mesh%elem(node, element)
                    endif
                enddo
            enddo
            !Save the start point as first entry to the coordinate list
            lsi_coord_number = lsi_coord_number + 1
            lsi_coords(lsi_coord_number, 1) = mesh%coord(1, start_point)
            lsi_coords(lsi_coord_number, 2) = mesh%coord(2, start_point)

            call fracture(mesh, lsi_loc(frac), start_point, start_element, end_point, lsi_coords, lsi_number, lsi_coord_number, frac_number, tmp_lsi_elements, tmp_property)

            !Save necessary information to plot the lines. The first entry is the number of the slip interface from which the current,
            !fracture starts. The second entry is the total amount of coordinates up to the current fracture and equally the number
            !final coordinate in the current fracture. The last entry has the number of slip interfaces each fracture is made from.
            frac_points(frac,1) = (lsi_coord_number - frac_number)+1
            frac_points(frac,2) = lsi_coord_number
            frac_points(frac,3) = frac_number
        enddo

        call writeVTKfractures(trim(outpath)//'fractures.vtk', lsi_coords(:lsi_coord_number,:), frac_points, size(lsi_loc))

        allocate(lsi(lsi_number))

        lsipar%lsin = lsi_number
        do i = 1, lsi_number
            lsi(i)%elements(1) = tmp_lsi_elements(i, 1)
            lsi(i)%elements(2) = tmp_lsi_elements(i, 2)
            lsi(i)%prop = tmp_property(i)
        end do

        deallocate(lsi_coords)
        deallocate(tmp_lsi_elements)
    end subroutine

    function point_distance(mesh, point1, point2)
        !input
        type(meshVar) :: mesh
        integer, intent(in) :: point1, point2
        !output
        real(kind=custom_real) :: point_distance
        !local
        real(kind=custom_real) :: x1, z1, x2, z2

        x1 = mesh%coord(1,point1)
        z1 = mesh%coord(2,point1)
        x2 = mesh%coord(1,point2)
        z2 = mesh%coord(2,point2)

        point_distance = sqrt((x1 - x2)**2 + (z1 - z2)**2)
    end function

    function dtl(mesh, lsi_loc, point)! dtl = distance to line
        !input
        type(meshVar) :: mesh
        type(lsi_location) :: lsi_loc
        integer, intent(in) :: point
        !output
        real(kind = custom_real) ::  dtl
        !local
        real(kind = custom_real) ::  s, qx, qz, lx, lz


        qx = lsi_loc%end_x - lsi_loc%start_x
        qz = lsi_loc%end_z - lsi_loc%start_z
        s = (qx * (mesh%coord(1,point) - lsi_loc%start_x) + qz * (mesh%coord(2,point) - lsi_loc%start_z))/(qx**2 + qz**2)
        lx = lsi_loc%start_x + s*qx
        lz = lsi_loc%start_z + s*qz
        dtl = sqrt((mesh%coord(1,point) - lx)**2 + (mesh%coord(2,point) - lz)**2)
    end function


    subroutine fracture(mesh, lsi_loc, start_point, start_element, end_point, lsi_coords, lsi_number, lsi_coord_number, frac_number, tmp_lsi_elements, tmp_property)
        !input
        type(meshVar) :: mesh
        type(lsi_location) :: lsi_loc
        integer, intent(in) :: start_point, end_point, start_element
        !output
        integer :: lsi_number, lsi_coord_number, frac_number
        integer, dimension(:,:) :: tmp_lsi_elements
        integer, dimension(:) :: tmp_property
        real(kind=custom_real), dimension(:,:) :: lsi_coords
        !local
        integer :: i,j, elem, node, node2                !iterators
        integer :: nelem, node_number, pt                !counters
        integer :: tmpelem, tmppoint, neighbour, element !(temporary) storage variables
        integer :: next_point, point
        integer, dimension(:), allocatable :: tmp, elements, nodes, pospoints
        real(kind=custom_real) :: dist_end, dist_point, tmp_dist
        logical :: in, found
        logical :: newneighbour

        allocate(tmp(20))
        element = start_element
        point = start_point
        next_point = 0
        !below here will be done until prev_point == end_point
        do while (point /= end_point)
            in = .false.
            newneighbour = .true.
            tmp = 0
            !Find all elements that contain the point to be tested.
            tmp(1) = element
            nelem = 1
            do while (newneighbour)
                tmpelem = tmp(nelem)
                found = .false.
                do i = 1, 3
                    neighbour = mesh%neighbor(i, tmpelem)
                    if (neighbour == 0) cycle
                    do node = 1, 3
                        if (mesh%elem(node, neighbour) == point) then
                            do j = 1, nelem
                                if (neighbour == tmp(j)) then
                                    in = .true.
                                    exit
                                else
                                    in = .false.
                                end if
                            end do
                            if (.not. in) then
                                found = .true.
                                nelem = nelem + 1
                                tmp(nelem) = neighbour
                            end if
                        end if
                    end do
                end do
                if (.not. found) newneighbour = .false.
            end do

            allocate(elements(nelem))
            elements = tmp(:nelem)
            tmp = 0
            in  = .false.

            !Find all points in the above elements execpt the point to be tested
            node_number = 0
            do elem = 1, nelem
                tmpelem = elements(elem)
                do node = 1, 3
                    if (mesh%elem(node, tmpelem) == point) then
                        do node2 = 1, 3
                            if (point /= mesh%elem(node2, tmpelem)) then
                                do i = 1, node_number
                                    if (tmp(i) == mesh%elem(node2, tmpelem)) then
                                        in = .true.
                                    exit
                                    else
                                        in = .false.
                                    end if
                                end do
                                if (.not. in) then
                                    node_number = node_number + 1
                                    tmp(node_number) = mesh%elem(node2, tmpelem)
                                end if
                            end if
                        end do
                    end if
                end do
            end do
            allocate(nodes(node_number))
            nodes(:) = tmp(:node_number)
            tmp = 0
            pt = 0
            dist_end = point_distance(mesh, point, end_point)
            do node = 1, node_number
                tmppoint = nodes(node)
                dist_point = point_distance(mesh, tmppoint, end_point)
                if (dist_point < dist_end) then
                    pt = pt + 1
                    tmp(pt) = tmppoint
                end if
            end do
            allocate(pospoints(pt))
            pospoints = tmp(:pt)
            dist_point = 100.
            do node = 1, pt
                tmp_dist = dtl(mesh, lsi_loc, pospoints(node))
                if (tmp_dist < dist_point) then
                    dist_point = tmp_dist
                    next_point = pospoints(node)
                end if
            end do

            !Find the elements that contain both point -> those are the elements that have the lsi.
            j = 0
            lsi_number = lsi_number + 1
            lsi_coord_number = lsi_coord_number + 1
            frac_number = frac_number + 1
            do elem = 1, nelem
                tmpelem = elements(elem)
                do node = 1, 3
                    if (next_point == mesh%elem(node, tmpelem)) then
                        j = j + 1
                        tmp_lsi_elements(lsi_number, j) = tmpelem
                    end if
                end do
            end do
            tmp_property(lsi_number) = lsi_loc%property
            lsi_coords(lsi_coord_number, 1) = mesh%coord(1, next_point)
            lsi_coords(lsi_coord_number, 2) = mesh%coord(2, next_point)

            point = next_point
            element = tmp_lsi_elements(lsi_number, 1)
            deallocate(elements, nodes, pospoints)
        end do
    end subroutine

    subroutine copyMesh(oldmesh, newmesh)

        type(meshVar) :: oldmesh, newmesh

        newmesh%nglob = oldmesh%nglob
        newmesh%nelem = oldmesh%nelem
        newmesh%nmat = oldmesh%nmat
        newmesh%ncoord = oldmesh%ncoord
        newmesh%pinterfaces = oldmesh%pinterfaces
        newmesh%mpi_nn = oldmesh%mpi_nn
        newmesh%mpi_nnmax = oldmesh%mpi_nnmax
        newmesh%mpi_ne = oldmesh%mpi_ne
        newmesh%mpi_nemax = oldmesh%mpi_nemax
        newmesh%nproc = oldmesh%nproc
        newmesh%has_src = oldmesh%has_src
        newmesh%has_rec = oldmesh%has_rec
        newmesh%dtfactor = oldmesh%dtfactor
        newmesh%nelx = oldmesh%nelx
        newmesh%nelz = oldmesh%nelz
        newmesh%minGLL = oldmesh%minGLL

        newmesh%invVdm = oldmesh%invVdm
        newmesh%Vdm = oldmesh%Vdm
        newmesh%Dr = oldmesh%Dr
        newmesh%Ds = oldmesh%Ds
        newmesh%Drw = oldmesh%Drw
        newmesh%Dsw = oldmesh%Dsw
        newmesh%Vdmr = oldmesh%Vdmr
        newmesh%Vdms = oldmesh%Vdms
        newmesh%mass = oldmesh%mass
        newmesh%lift = oldmesh%lift
        newmesh%fmask = oldmesh%fmask
        newmesh%fmaskv = oldmesh%fmaskv

        allocate(newmesh%vx(newmesh%nglob),newmesh%vz(newmesh%nglob))
        allocate(newmesh%matval(newmesh%nmat,7))
        allocate(newmesh%rx(newmesh%nglob), newmesh%rz(newmesh%nglob), newmesh%sx(newmesh%nglob), newmesh%sz(newmesh%nglob))
        allocate(newmesh%jacobian(newmesh%nglob))
        allocate(newmesh%mat(newmesh%nelem))
        allocate(newmesh%ibool(Np,newmesh%nelem))
        allocate(newmesh%nx(3*NGLL,newmesh%nelem), newmesh%nz(3*NGLL,newmesh%nelem), newmesh%sj(3*NGLL,newmesh%nelem),&
             newmesh%fscale(3*NGLL,newmesh%nelem))
        allocate(newmesh%coord(2,newmesh%ncoord))
        allocate(newmesh%elem(3,newmesh%nelem))
        allocate(newmesh%neighbor(3,newmesh%nelem))
        allocate(newmesh%face(3,newmesh%nelem))
        allocate(newmesh%vp(newmesh%nelem),newmesh%vs(newmesh%nelem),newmesh%rho(newmesh%nelem),newmesh%qp(newmesh%nelem),newmesh%qs(newmesh%nelem),newmesh%mu(newmesh%nelem),newmesh%lambda(newmesh%nelem),newmesh%pml(newmesh%nelem))
        allocate(newmesh%vpu(newmesh%nelem),newmesh%vsu(newmesh%nelem))
        allocate(newmesh%ymu(nMB,newmesh%nelem),newmesh%ylambda(nMB,newmesh%nelem), newmesh%wl(nMB,newmesh%nelem))
        allocate(newmesh%ibt(NGLL,3,newmesh%nelem),newmesh%ibn(NGLL,3,newmesh%nelem))
        allocate(newmesh%loc2glob_nodes(newmesh%ncoord))
        allocate(newmesh%mpi_interface(4,3,newmesh%nelem))
        allocate(newmesh%mpi_neighbor(newmesh%mpi_nn))
        allocate(newmesh%mpi_connection(newmesh%mpi_nn,newmesh%mpi_ne,2))
        allocate(newmesh%mpi_ninterface(newmesh%mpi_nn))
        allocate(newmesh%mpi_ibool(3,newmesh%mpi_nemax))
        allocate(newmesh%mpi_ibt(NGLL,3,newmesh%nelem))
        allocate(newmesh%mpi_icon(newmesh%mpi_nnmax))

        allocate(newmesh%elemNo(newmesh%nelem))
        allocate(newmesh%imp_vp(newmesh%nelem))
        allocate(newmesh%imp_vs(newmesh%nelem))
        allocate(newmesh%mpi_imp_vp(newmesh%pinterfaces, 3))
        allocate(newmesh%mpi_imp_vs(newmesh%pinterfaces, 3))
        allocate(newmesh%vol(newmesh%nelem))
        allocate(newmesh%srcrecmask(newmesh%nelem))
        allocate(newmesh%smooth_A(newmesh%nelem,4))

        newmesh%vx = oldmesh%vx
        newmesh%vz = oldmesh%vz
        newmesh%matval = oldmesh%matval
        newmesh%rx = oldmesh%rx
        newmesh%rz = oldmesh%rz
        newmesh%sx = oldmesh%sx
        newmesh%sz = oldmesh%sz
        newmesh%jacobian = oldmesh%jacobian
        newmesh%mat = oldmesh%mat
        newmesh%ibool = oldmesh%ibool
        newmesh%nx = oldmesh%nx
        newmesh%nz = oldmesh%nz
        newmesh%sj = oldmesh%sj
        newmesh%fscale = oldmesh%fscale
        newmesh%coord = oldmesh%coord
        newmesh%elem = oldmesh%elem
        newmesh%neighbor = oldmesh%neighbor
        newmesh%face = oldmesh%face
        newmesh%poroelastic = oldmesh%poroelastic
        if (newmesh%poroelastic) then
            newmesh%nfluids = oldmesh%nfluids
            allocate(newmesh%A(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem),&
                     newmesh%B(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem),&
                     newmesh%E(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem),&
                     newmesh%AM(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem),&
                     newmesh%AP(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem),&
                     newmesh%BM(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem),&
                     newmesh%BP(5+3*newmesh%nfluids,5+3*newmesh%nfluids,newmesh%nelem))
            newmesh%A = oldmesh%A
            newmesh%B = oldmesh%B
            newmesh%E = oldmesh%E
            newmesh%AM = oldmesh%AM
            newmesh%AP = oldmesh%AP
            newmesh%BM = oldmesh%BM
            newmesh%BP = oldmesh%BP
            newmesh%vp = oldmesh%vp
            newmesh%vs = oldmesh%vs
            newmesh%vpu= oldmesh%vpu
            newmesh%vsu= oldmesh%vsu
        else
            allocate(newmesh%rho(newmesh%nelem),newmesh%qp(newmesh%nelem),newmesh%qs(newmesh%nelem),newmesh%mu(newmesh%nelem),newmesh%lambda(newmesh%nelem))
            newmesh%vp = oldmesh%vp
            newmesh%vs = oldmesh%vs
            newmesh%vpu= oldmesh%vpu
            newmesh%vsu= oldmesh%vsu
            newmesh%rho = oldmesh%rho
            newmesh%qp = oldmesh%qp
            newmesh%qs = oldmesh%qs
            newmesh%mu = oldmesh%mu
            newmesh%lambda = oldmesh%lambda
        endif
        newmesh%ymu = oldmesh%ymu
        newmesh%lambda = oldmesh%lambda
        newmesh%wl = oldmesh%wl
        newmesh%pml = oldmesh%pml
        newmesh%pmlxmin = oldmesh%pmlxmin
        newmesh%pmlxmax = oldmesh%pmlxmax
        newmesh%pmlzmin = oldmesh%pmlzmin
        newmesh%pmlzmax = oldmesh%pmlzmax
        newmesh%ibt = oldmesh%ibt
        newmesh%ibn = oldmesh%ibn
        newmesh%loc2glob_nodes = oldmesh%loc2glob_nodes
        newmesh%mpi_interface = oldmesh%mpi_interface
        newmesh%mpi_neighbor = oldmesh%mpi_neighbor
        newmesh%mpi_connection = oldmesh%mpi_connection
        newmesh%mpi_ninterface = oldmesh%mpi_ninterface
        newmesh%mpi_ibool = oldmesh%mpi_ibool
        newmesh%mpi_ibt = oldmesh%mpi_ibt
        newmesh%mpi_icon = oldmesh%mpi_icon

        newmesh%elemNo = oldmesh%elemNo
        newmesh%imp_vp = oldmesh%imp_vp
        newmesh%imp_vs = oldmesh%imp_vs
        newmesh%mpi_imp_vp = oldmesh%mpi_imp_vp
        newmesh%mpi_imp_vs = oldmesh%mpi_imp_vs
        newmesh%vol = oldmesh%vol
        newmesh%srcrecmask = oldmesh%srcrecmask
        newmesh%smooth_A = oldmesh%smooth_A
        newmesh%sDT = oldmesh%sDT

        close(27)

    end subroutine

end module meshMod
