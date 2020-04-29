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
module plotMod
    use meshMod
    use vtkMod
    use errorMessage
    use constantsMod
    use progressbarMod

    implicit none
    ! module to save some files to disk for plotting purpose
    contains

    subroutine plotPoints2d(x,z,filename)
        implicit none
        real(kind=CUSTOM_REAL), dimension(:) :: x,z
        character(len=80) :: filename
        integer :: i

        open(unit=19,file=trim(filename))
        do i=1,size(x)
            write(19,*) x(i),z(i)
        end do
        close(19)
    end subroutine plotPoints2d

    subroutine writePlotBinaries(data, path, name, myrank, it)
        !input
        character(len=*) :: name                        !Name of the variable
        integer :: myrank                               !Number of the core
        integer :: it                                   !current timestep
        real(kind=custom_real), dimension(:) :: data    !Vector to store the data to be plotted
        !Local
        character(len=6) :: number_of_rank
        character(len=7) :: timestep
        character(len=50) :: filename, path

        write(number_of_rank, "(i6.6)") myrank+1
        write(timestep, "(i7.7)") it

        filename = trim(path)//trim(name)//"_"//trim(number_of_rank)//"_it"//trim(timestep)//".bin"
        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        write(27) data
        close(27)
    end subroutine

    subroutine writePlotBinariesInv(data, path, name, myrank, iter)
        !input
        character(len=*) :: name                        !Name of the variable
        integer :: myrank                               !Number of the core
        integer :: iter                                 !current iterationstep
        real(kind=custom_real), dimension(:) :: data    !Vector to store the data to be plotted
        !Local
        character(len=6) :: number_of_rank
        character(len=7) :: iterstep
        character(len=50) :: filename, path

        write(number_of_rank, "(i6.6)") myrank+1
        write(iterstep, "(i3.3)") iter

        filename = trim(path)//trim(name)//"_"//trim(number_of_rank)//"_iter"//trim(iterstep)//".bin"
        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        write(27) data
        close(27)
    end subroutine

    subroutine readPlotBinaries(frame, leln, ueln, data, path, name, proc, it, errmsg)
        !input
        type(error_message) :: errmsg
        character(len=*) :: name                        !Name of the variable
        integer :: proc                                 !Number of the core
        integer :: it                                   !current time step
        !These parameters have to do with the parallelisation of the code. All the information has to be merged into 1 array,
        !thus these counters are necessary
        integer :: frame                                !c in collectMovie
        integer :: leln                                 !lower_element_number (d in collectMovie)
        integer :: ueln                                 !upper_element_number (e in collectMovie)
        real(kind=custom_real), dimension(:,:) :: data  !Array to store the data to be plotted
        character(len=16) :: myname = "readPlotBinaries"

        !Local
        character(len=6) :: number_of_rank
        character(len=7) :: timestep
        character(len=50) :: filename, path
        logical :: file_exists

        write(number_of_rank, "(i6.6)") proc
        write(timestep, "(i7.7)") it

        filename = trim(path)//trim(name)//"_"//trim(number_of_rank)//"_it"//trim(timestep)//".bin"
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
            read(27) data(leln:ueln,frame )
            close(27)
        else
            call add(errmsg, 2, "Error in reading plot files! File "//trim(filename)//" does not exist!", myname, filename)
            call print(errmsg)
            stop
        end if
    end subroutine

    subroutine readPlotBinariesInv(leln, ueln, data, path, name, proc, iter, errmsg)
        !input
        type(error_message) :: errmsg
        character(len=*) :: name                        !Name of the variable
        integer :: proc                                 !Number of the core
        integer :: iter                                 !current iteration step
        !These parameters have to do with the parallelisation of the code. All the information has to be merged into 1 array,
        !thus these counters are necessary
        integer :: leln                                 !lower_element_number (d in collectMovie)
        integer :: ueln                                 !upper_element_number (e in collectMovie)
        real(kind=custom_real), dimension(:) :: data  !Array to store the data to be plotted
        character(len=19) :: myname = "readPlotBinariesInv"

        !Local
        character(len=6) :: number_of_rank
        character(len=7) :: iterstep
        character(len=50) :: filename, path
        logical :: file_exists

        write(number_of_rank, "(i6.6)") proc
        write(iterstep, "(i3.3)") iter

        filename = trim(path)//trim(name)//"_"//trim(number_of_rank)//"_iter"//trim(iterstep)//".bin"
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
            read(27) data(leln:ueln)
            close(27)
        else
            call add(errmsg, 2, "Error in reading plot files! File "//trim(filename)//" does not exist!", myname, filename)
            call print(errmsg)
            stop
        end if
    end subroutine


    subroutine generateTriMeshFiles(db, data, path, name, nframe, nproc, glob_elem, glob_coord2, make_progress)
        !input
        type(meshVar), dimension(:) :: db               !Databases containing each the information of on core
        integer :: nframe                               !Total number of movieframes
        integer :: nproc                                !Total number of cores used
        integer, dimension(:,:) :: glob_elem
        real(kind=custom_real), dimension(:,:) :: data  !Dataarray
        real(kind=custom_real), dimension(:,:) :: glob_coord2
        character(len=*) :: name                        !Addition to the filename
        logical :: make_progress
        !local
        integer :: element
        integer :: point
        integer :: it
        integer :: iproc
        !These counters are used create a single file from all database files.
        integer :: total_element_number                 !The local element number is added up for each processor to ensure continuity
        integer :: total_point_number                   !The total point number is created here. nglob = nelement*Np
        real(kind=custom_real), dimension(size(data(:,1)), size(data(1,:))) :: tempdata
        character(len=50) :: filename, path
        character(len=7) :: timestep


        !get data and average over element
        total_element_number = 0
        total_point_number = 0
        tempdata(:,:) = 0

        !This loop goes for each core through all movie binaries to create a single array containing all the data for plotting.
        do iproc=1, nproc
            do it=1, nframe
                do element=1, db(iproc)%nelem
                    do point=1, Np
                        tempdata(element+total_element_number, it) = tempdata(element+total_element_number, it) &
                                                                    + data(db(iproc)%ibool(point, element)+total_point_number, it)
                    end do
                end do
            end do
            if (iproc < nproc) then
                total_element_number = total_element_number + db(iproc)%nelem
                total_point_number   = total_point_number + db(iproc)%nglob
            end if
        end do
        tempdata=tempdata/Np

        do it = 1,nframe
            write(timestep,"(i7.7)") it
            filename = "movie_element_"//trim(name)//"_it"//trim(timestep)//".vtk"
            call writeVtkTriMeshRealdata(trim(path)//filename, glob_elem, glob_coord2, tempdata(:,it), name)
            if (make_progress) call progress(it, nframe)
        end do
    end subroutine

    subroutine generateTriMeshFilesAdj(db, data, path, name, nframe, nproc, glob_elem, glob_coord2, make_progress)
        !input
        type(meshVar), dimension(:) :: db               !Databases containing each the information of on core
        integer :: nframe                               !Total number of movieframes
        integer :: nproc                                !Total number of cores used
        integer, dimension(:,:) :: glob_elem
        real(kind=custom_real), dimension(:,:) :: data  !Dataarray
        real(kind=custom_real), dimension(:,:) :: glob_coord2
        character(len=*) :: name                        !Addition to the filename
        logical :: make_progress
        !local
        integer :: element
        integer :: point
        integer :: it
        integer :: iproc
        !These counters are used create a single file from all database files.
        integer :: total_element_number                 !The local element number is added up for each processor to ensure continuity
        integer :: total_point_number                   !The total point number is created here. nglob = nelement*Np
        real(kind=custom_real), dimension(size(data(:,1)), size(data(1,:))) :: tempdata
        character(len=50) :: filename, path
        character(len=7) :: timestep


        !get data and average over element
        total_element_number = 0
        total_point_number = 0
        tempdata(:,:) = 0

        !This loop goes for each core through all movie binaries to create a single array containing all the data for plotting.
        do iproc=1, nproc
            do it=1, nframe
                do element=1, db(iproc)%nelem
                    do point=1, Np
                        tempdata(element+total_element_number, it) = tempdata(element+total_element_number, it) &
                                                                    + data(db(iproc)%ibool(point, element)+total_point_number, it)
                    end do
                end do
            end do
            if (iproc < nproc) then
                total_element_number = total_element_number + db(iproc)%nelem
                total_point_number   = total_point_number + db(iproc)%nglob
            end if
        end do
        tempdata=tempdata/Np

        do it = 1,nframe
            write(timestep,"(i7.7)") it
            filename = "movie_element_"//trim(name)//"_it"//trim(timestep)//".vtk"
            call writeVtkTriMeshRealdata(trim(path)//filename, glob_elem, glob_coord2, tempdata(:,it), name)
            if (make_progress) call progress(it, nframe)
        end do
    end subroutine

    subroutine generatePointFiles(data, path, name, xyzplot, nframe, make_progress)
        !input
        integer :: nframe                                 !Total number of movieframes
        real(kind=custom_real), dimension(:,:) :: data    !Dataarray
        real(kind=custom_real), dimension(:,:) :: xyzplot !Array containing the coordinate information
        character(len=*) :: name                          !Addition to the filename
        logical :: make_progress
        !local
        integer :: it
        character(len=50) :: filename, path
        character(len=7) :: timestep

        do it = 1, nframe
            write(timestep,"(i7.7)") it
            filename = "movie_points_"//trim(name)//"_it"//trim(timestep)//".vtk"
            call writeVtkNodesRealData(trim(path)//filename, xyzplot(1,:),xyzplot(2,:),xyzplot(3,:), data(:,it), name)
            if (make_progress) call progress(it, nframe)
        end do
    end subroutine

    subroutine generatePointFilesAdj(data, path, name, xyzplot, nframe, make_progress)
        !input
        integer :: nframe                                 !Total number of movieframes
        real(kind=custom_real), dimension(:,:) :: data    !Dataarray
        real(kind=custom_real), dimension(:,:) :: xyzplot !Array containing the coordinate information
        character(len=*) :: name                          !Addition to the filename
        logical :: make_progress
        !local
        integer :: it
        character(len=50) :: filename, path
        character(len=7) :: timestep

        do it = 1, nframe
            write(timestep,"(i7.7)") nframe-it+1
            filename = "movie_points_"//trim(name)//"_it"//trim(timestep)//".vtk"
            call writeVtkNodesRealData(trim(path)//filename, xyzplot(1,:),xyzplot(2,:),xyzplot(3,:), data(:,it), name)
            if (make_progress) call progress(it, nframe)
        end do
    end subroutine

    subroutine generateTriMeshFilesElem(data, path, name, nframe, glob_elem, glob_coord2, make_progress)
        !input
        integer :: nframe                               !Total number of movieframes
        integer, dimension(:,:) :: glob_elem
        real(kind=custom_real), dimension(:,:) :: data  !Dataarray
        real(kind=custom_real), dimension(:,:) :: glob_coord2
        character(len=*) :: name                        !Addition to the filename
        logical :: make_progress
        !local
        !These counters are used create a single file from all database files.
        character(len=50) :: filename, path
        character(len=7) :: timestep

        write(timestep,"(i7.7)") nframe
        filename = "movie_element_"//trim(name)//"_it"//trim(timestep)//".vtk"
        call writeVtkTriMeshRealdata(trim(path)//filename, glob_elem, glob_coord2, data(:,1), name)
        if (make_progress) call progress(nframe, nframe)

    end subroutine

end module plotMod
