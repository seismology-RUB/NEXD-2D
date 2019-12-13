!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2019 Andre Lamert (Ruhr-Universität Bochum, GER)
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
program solver
    ! Program to solve the 2D (poro)elastic wave equation with the discontinous galerkin methode in 2D on triangular mesh
    use parameterMod
    use meshMod
    use timeloopMod
    use sourceReceiverMod
    use errorMessage
    use slipInterfaceMod

    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type (meshVar) :: mesh
    type (srcVar) :: src
    type (recVar) :: rec
    type(interface_spec), dimension(:), allocatable :: lsi_spec
    type(lsiVar), dimension(:), allocatable :: lsi
    type(elementToLSI), dimension(:), allocatable :: etolsi
    type(lsi_parameter) :: lsipar
    type(movie_parameter) :: movie
    integer :: world_size
    integer :: myrank
    logical :: file_exists
    character(len=80) :: filename
    character(len=6) :: myname = "solver"

    !Create new errormessagetrace
    call new(errmsg,myname)

    ! start MPI
    call init_mpi()
    call comm_size(world_size)
    call comm_rank(myrank)

    if (myrank==0) call writeLogo()

    ! read Parfile
    call readParfile(par, lsipar, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    call sync_mpi()

    if (par%log .and. myrank == 0) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                             start dg2d mpi solver                            |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if

    call sync_mpi()
    if (par%log) write(*,'(a40, i5, a35)')  "|                        Start MPI job: " , myrank, "                                  |"
    call sync_mpi()
    if (par%log .and. myrank == 0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
    call sync_mpi()

    call deallocMeshVar(mesh)

    write(filename,"('out/meshVar',i6.6)") myrank+1
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
        call readMeshVar(mesh,filename, errmsg)
    else
        call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
        call print(errmsg)
        call stop_mpi()
    end if

    call sync_mpi()
    if (par%log .and. myrank == 0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
    call sync_mpi()

    !LSI only
    if (lsipar%lsi) then
        write(filename,"('out/elementToLSI',i6.6)") myrank+1
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            call setup_lsi(mesh%nelem, etolsi, lsi, lsi_spec, lsipar, filename, myrank, errmsg)
        else
            call add(errmsg, 2, "Error in databases, file "//trim(filename)//" does not exist!", myname, filename)
            call print(errmsg)
            call stop_mpi()
        end if
        call sync_mpi()
        if (par%log .and. myrank == 0) write(*,'(a80)') "|------------------------------------------------------------------------------|"
        call sync_mpi()
    end if

    ! start timeloop
    call timeloop2d(par,mesh,src,rec, lsi, lsi_spec, lsipar, etolsi, movie, myrank, errmsg)

    if (par%log .and. myrank == 0) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                              end dg2d mpi solver                             |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if

    !LSI Only - Deallocate LSI Arrays
    if (lsipar%lsi) then
        if (allocated(lsi)) deallocate(lsi)
        if (allocated(lsi_spec)) deallocate(lsi_spec)
    endif
end program solver
