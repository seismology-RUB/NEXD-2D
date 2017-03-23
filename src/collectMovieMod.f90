!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2016 Elena Busch (Rice University, United States)
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
module collectMovieMod
    use parameterMod
    use meshMod
    use vtkMod
    use errorMessage
    use progressbarMod
    use plotMod
    
    implicit none

    contains

    subroutine collectMovie(par, movie, errmsg)
        !input
        type(parameterVar) :: par
        type(movie_parameter) :: movie
        type(error_message) :: errmsg
        !local
        type(meshVar),dimension(:), allocatable :: db
        integer, pointer :: nproc
        integer :: iproc
        integer, pointer :: nt
        integer :: nstep
        integer :: it,nframe,c,d,e,i
        integer, dimension(:), allocatable :: ndata,ndata2,ndata3
        ! file
        logical :: file_exists
        character(len=80) ::filename, outpath
        ! data
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: uplot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ux
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: uz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vplot
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vel_x
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: vel_z
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigmaxx
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigmazz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sigmaxz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyzplot
        !real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tempdata
        integer :: glob_nelem, glob_nglob, glob_ncoord
        real(kind=CUSTOM_REAL) , dimension(:,:), allocatable :: glob_coord, glob_coord2
        integer , dimension(:,:), allocatable :: glob_elem

        !For the errormessage:
        character(len=12) :: myname = "collectMovie"

        !Add module to the trace of the errormessage
        call addTrace(errmsg, myname)

        outpath = "out"

        if (par%log) then
            write (*, "(a80)") "|                               generating movie                               |"
            write (*, "(a80)") "|------------------------------------------------------------------------------|"
        end if
        ! load database

        nframe=par%nt/movie%frame

        allocate(db(par%nproc))
        allocate(ndata(par%nproc))
        allocate(ndata2(par%nproc))
        allocate(ndata3(par%nproc))

        do iproc=1,par%nproc
            write(filename,"('/meshVar',i6.6)") iproc
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readMeshVar(db(iproc),filename, errmsg)

            else
                call add(errmsg, 2, "Error in the database. File "//trim(filename)//" do not exist. Abort...", myname)
                if (.level.errmsg == 2) then; call print(errmsg); stop; endif
           end if
        end do
        if (par%log) write (*, "(a80)") "|------------------------------------------------------------------------------|"
        
        ! get pointers    
        glob_nelem=0
        glob_nglob=0
        glob_ncoord=0
        do iproc=1,par%nproc
            glob_nelem=glob_nelem+db(iproc)%nelem
            glob_nglob=glob_nglob+db(iproc)%nglob
            glob_ncoord=glob_ncoord+db(iproc)%ncoord ! only guess for number of coords
            ndata(iproc) = db(iproc)%nglob
            ndata2(iproc) = db(iproc)%ncoord
            ndata3(iproc) = db(iproc)%nelem
        end do
        
        allocate(xyzplot(3,glob_nglob))

        ! read coords
        d=1
        e=ndata(1)
        do iproc=1,par%nproc
            xyzplot(1,d:e) = db(iproc)%vx(:)
            xyzplot(2,d:e) = 0.0
            xyzplot(3,d:e) = db(iproc)%vz(:)
            if (iproc < par%nproc) then
                d=d+ndata(iproc)
                e=e+ndata(iproc+1)
            end if
        end do

        !build global mesh
        allocate(glob_elem(3,glob_nelem))
        allocate(glob_coord(3,glob_ncoord))

        !build global elements
        c=1
        do iproc=1,par%nproc
            do i=1,db(iproc)%nelem
                glob_elem(1,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(1,i))
                glob_elem(2,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(2,i))
                glob_elem(3,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(3,i))
                c=c+1
            end do
        end do

        !build global coordinates
        glob_coord=-1e17
        do iproc=1,par%nproc
            do i=1,db(iproc)%ncoord
                glob_coord(1,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(1,i)
                glob_coord(2,db(iproc)%loc2glob_nodes(i)) = 0.0
                glob_coord(3,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(2,i)
            end do
        end do

        c=0
        do i=1,glob_ncoord
            if (glob_coord(1,i)>-1e17) then
                c=c+1
            end if
        end do
        allocate(glob_coord2(3,c))
        c=0
        do i=1,glob_ncoord
            if (glob_coord(1,i)>-1e17) then
                c=c+1
                glob_coord2(:,c)=glob_coord(:,i)
            end if
        end do

        !do the things that involve allocating memory
        if (movie%displacement) then
            if (par%log) write (*,*) "Creating Displacement files"
            allocate(uplot(glob_nglob,nframe))
            allocate(ux(glob_nglob,nframe))
            allocate(uz(glob_nglob,nframe))
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,uplot,"normU",glob_nglob,nframe,ndata)
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,ux,"ux",glob_nglob,nframe,ndata)
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,uz,"uz",glob_nglob,nframe,ndata)
            if (movie%trimesh) then
                if (par%log) write(*,*) "Creating trimesh displacement files"
                call generateTriMeshFiles(db, uplot, "norm_u", nframe, par%nproc, glob_elem, glob_coord2)
                call generateTriMeshFiles(db, ux, "ux", nframe, par%nproc, glob_elem, glob_coord2)
                call generateTriMeshFiles(db, uz, "uz", nframe, par%nproc, glob_elem, glob_coord2)
            end if
            if (movie%points) then
                if (par%log) write(*,*) "Creating points displacement files"
                call generatePointFiles(uplot, "normu", xyzplot, nframe)
                call generatePointFiles(ux, "ux", xyzplot, nframe)
                call generatePointFiles(uz, "uz", xyzplot, nframe)
            end if
            deallocate(uplot)
            deallocate(ux)
            deallocate(uz)
        end if
        
        if (movie%velocity) then
            if (par%log) write (*,*) "Creating Velocity files"
            allocate(vplot(glob_nglob,nframe))
            allocate(vel_x(glob_nglob,nframe))
            allocate(vel_z(glob_nglob,nframe))
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,vplot,"normV",glob_nglob,nframe,ndata)
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,vel_x,"vx",glob_nglob,nframe,ndata)
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,vel_z,"vz",glob_nglob,nframe,ndata)
            if (movie%trimesh) then
                if (par%log) write(*,*) "Creating trimesh Velocity files"
                call generateTriMeshFiles(db, vplot, "norm_v", nframe, par%nproc, glob_elem, glob_coord2)
                call generateTriMeshFiles(db, vel_x, "vx", nframe, par%nproc, glob_elem, glob_coord2)
                call generateTriMeshFiles(db, vel_z, "vz", nframe, par%nproc, glob_elem, glob_coord2)
            end if
            if (movie%points) then
                if (par%log) write(*,*) "Creating points Velocity files"
                call generatePointFiles(vplot, "normv", xyzplot, nframe)
                call generatePointFiles(vel_x, "vx", xyzplot, nframe)
                call generatePointFiles(vel_z, "vz", xyzplot, nframe)
            end if
            deallocate(vplot)
            deallocate(vel_x)
            deallocate(vel_z)
        end if
        
        if (movie%stress) then
            if (par%log) write (*,*) "Creating Stress files"
            allocate(sigmaxx(glob_nglob,nframe))
            allocate(sigmazz(glob_nglob,nframe))
            allocate(sigmaxz(glob_nglob,nframe))
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,sigmaxx,"sigmaXX",glob_nglob,nframe,ndata)
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,sigmazz,"sigmaZZ",glob_nglob,nframe,ndata)
            call readplotbin(par%nproc,movie%frame,par%nt,errmsg,sigmaxz,"sigmaXZ",glob_nglob,nframe,ndata)
            if (movie%trimesh) then
                if (par%log) write(*,*) "Creating trimesh stress files"
                call generateTriMeshFiles(db, sigmaxx, "sigmaxx", nframe, par%nproc, glob_elem, glob_coord2)
                call generateTriMeshFiles(db, sigmazz, "sigmazz", nframe, par%nproc, glob_elem, glob_coord2)
                call generateTriMeshFiles(db, sigmaxz, "sigmaxz", nframe, par%nproc, glob_elem, glob_coord2)
            end if
            if (movie%points) then
                if (par%log) write(*,*) "Creating points stress files"
                call generatePointFiles(sigmaxx, "sigmaxx", xyzplot, nframe)
                call generatePointFiles(sigmazz, "sigmazz", xyzplot, nframe)
                call generatePointFiles(sigmaxz, "sigmaxz", xyzplot, nframe)
            end if
            deallocate(sigmaxx)
            deallocate(sigmazz)
            deallocate(sigmaxz)
        end if

        if (allocated(db)) deallocate(db)
        if (allocated(ndata)) deallocate(ndata)
        if (allocated(ndata2)) deallocate(ndata2)
        if (allocated(ndata3)) deallocate(ndata3)
        if (allocated(xyzplot)) deallocate(xyzplot)
        if (allocated(glob_elem)) deallocate(glob_elem)
        if (allocated(glob_coord)) deallocate(glob_coord)
    end subroutine collectMovie
    
    subroutine readplotbin(nproc,nstep,nt,errmsg,dat,name1,glob_nglob,nframe,ndata)
        !input
        integer :: nproc
        integer :: nt
        integer :: glob_nglob
        integer :: nframe
        integer :: nstep
        type(error_message) :: errmsg
        character(len=*) :: name1
        
        real(kind=CUSTOM_REAL), dimension(glob_nglob,nframe):: dat
        integer, dimension(:):: ndata
        !local
        integer :: iproc, c, it,e,d
        
        d=1
        e=ndata(1)
        
        do iproc=1,nproc
            c=1
            do it=nstep,nt,nstep
                if (name1 == 'uz') then
                    call readPlotBinaries(d, e, c, dat, name1, iproc, it, errmsg)
                else
                    call readPlotBinaries(c, d, e, dat, name1, iproc, it, errmsg)
                end if
            c=c+1
            end do
            if (iproc < nproc) then
                d=d+ndata(iproc)
                e=e+ndata(iproc+1)
            end if
        end do
    end subroutine readplotbin
end module collectMovieMod
