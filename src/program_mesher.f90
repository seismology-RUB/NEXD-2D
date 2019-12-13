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
program mesher
    !program to mesh the given cubit files to build a MPI version of the code
    use parameterMod
    use meshMod
    use plotMod
    use errorMessage
    use logoMod
    use slipInterfaceMod

    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type (meshVar) :: mesh
    type(lsi_parameter) :: lsipar
    type(movie_parameter) :: movie
    character(len=80) :: filename
    integer :: myrank = 0
    character(len=6) :: myname = "mesher"

    !Create new errormessagetrace
    call new(errmsg,myname)

    if (par%log) call writeLogo()

    ! read Parfile
    call readParfile(par, lsipar, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    ! log?
    if (par%log) then
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
        write(*,'(a80)') "|                     start mesher for mpi version of dg2d                     |"
        write(*,'(a80)') "|------------------------------------------------------------------------------|"
    end if

    ! create the mesh and write it to database files
    call createRegularMesh(mesh, par, lsipar, errmsg)

    ! need meshVar to plot the mesh ( not good code but to debug)
    ! plot mesh
    if (par%log) write(*,'(a80)') "|                                  plot mesh                                   |"
    call triangulation_order3_plot ( "out/mesh.ps", mesh%ncoord, dble(mesh%coord), mesh%nelem, mesh%elem, 2, 2 )
    filename = "out/pointsmesh.txt"
    call plotPoints2d(mesh%vx,mesh%vz,filename)
    if (par%log) write(*,'(a80)') "|------------------------------------------------------------------------------|"
    if (.level.errmsg == 1) call print(errmsg)
    call deallocMeshvar(mesh)
end program mesher
