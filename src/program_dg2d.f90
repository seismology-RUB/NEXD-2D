!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
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
program dg2d
    ! Program to solve the 2D elastic wave equation with the discontinous galerkin methode in 2D on triangular mesh
    !Is this program still in use? TM TM
    use constantsMod
    use parameterMod
    use gllMod
    use warpfactorMod
    use nodesMod
    use triTrafoMod
    use simplexMod
    use vandermondeMod
    use dmatricesMod
    use meshMod
    use plotMod
    use timeloopMod
    use sourceReceiverMod
    use errorMessage


    implicit none

    type(error_message) :: errmsg
    type (parameterVar) :: par
    type (meshVar) :: mesh
    type (sourceVar) :: src
    type (recVar) :: rec
    logical :: log
    character(len=80) :: filename
    character(len=4) :: myname = "dg2d"

    !Create new errormessagetrace
    call new(errmsg,myname)

    ! read Parfile
    call readParfile(par, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    ! log?
    log=ParLogS(par)
    if (log) then
        write(*,*) "--------------------------------------------------------------------------------------"
        write(*,*) "start dg2d"
        write(*,*) "--------------------------------------------------------------------------------------"
    end if

    ! create the mesh
    call createRegularMesh(mesh,par, errmsg)

    ! plot mesh
    write(*,*) "plot mesh"
    call triangulation_order3_plot ( "out/mesh.ps", mesh%ncoord, dble(mesh%coord), mesh%nelem, mesh%elem, 2, 2 )
    filename = "out/pointsmesh.txt"
    call plotPoints2d(mesh%vx,mesh%vz,filename)

    ! start timeloop
    call timeloop2d(par,mesh,src,rec)

    ! dealloc
    call deallocMeshVar(mesh)

    if (log) then
        write(*,*) "--------------------------------------------------------------------------------------"
        write(*,*) "end dg2d"
        write(*,*) "--------------------------------------------------------------------------------------"
    end if
end program dg2d
