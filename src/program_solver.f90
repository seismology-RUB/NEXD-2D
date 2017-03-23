!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
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
program solver
    ! Program to solve the 2D elastic wave equation with the discontinous galerkin methode in 2D on triangular mesh
    use parameterMod
    use meshMod
    use timeloopMod
    use sourceReceiverMod
    use errorMessage

    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type (meshVar) :: mesh
    type (srcVar) :: src
    type (recVar) :: rec
    type(movie_parameter) :: movie
    integer :: world_size
    integer :: myrank
    logical :: log
    character(len=6) :: myname = "solver"

    !Create new errormessagetrace
    call new(errmsg,myname)

    ! start MPI
    call init_mpi()
    call comm_size(world_size)
    call comm_rank(myrank)

    if (myrank==0) call writeLogo()

    ! read Parfile
    call readParfile(par, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
  
    ! log?
    if (par%log .and. myrank == 0) then
        write(*,*) "--------------------------------------------------------------------------------------"
        write(*,*) "start dg2d mpi solver"
        write(*,*) "--------------------------------------------------------------------------------------"
    end if

    ! start timeloop
    call timeloop2d(par,mesh,src,rec, movie, myrank, errmsg)

    if (par%log .and. myrank == 0) then
        write(*,*) "--------------------------------------------------------------------------------------"
        write(*,*) "end dg2d mpi solver"
        write(*,*) "--------------------------------------------------------------------------------------"
    end if
end program solver
