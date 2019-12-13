!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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
module vtkMod
    ! module to create files in VTK format which can be read by paraview
    use constantsMod

    implicit none

    contains

    subroutine writeVtkTetraMeshIntdata(filename, elem, coord, data, fieldname)
        character(len=*) :: filename, fieldname
        real(kind=CUSTOM_REAL), dimension(:,:) :: coord
        integer, dimension(:) :: data
        integer, dimension(:,:) :: elem
        integer :: nelem, ncoord
        integer :: i

        nelem=size(elem(1,:))
        ncoord=size(coord(1,:))

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Mesh with Data'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
        do i=1,ncoord
            write(19,'(3e18.6)') coord(1,i),coord(2,i),coord(3,i)
        enddo
        write(19,*) ""

        write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*5
        do i=1,nelem
            write(19,'(5i12)') 4 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1,elem(4,i)-1
        end do
        write(19,*) ""

        write(19,'(a,i12)') "CELL_TYPES ",nelem
        write(19,*) (10,i=1,nelem)
        write(19,*) ""

        write(19,'(a,i12)') "CELL_DATA ",nelem
        write(19,'(a,a,a)') "SCALARS ", fieldname, " integer"
        write(19,'(a)') "LOOKUP_TABLE default"
        do i = 1,nelem
            write(19,*) data(i)
        enddo
        write(19,*) ""
        close(19)
    end subroutine writeVtkTetraMeshIntdata

    subroutine writeVtkTetraMeshRealdata(filename, elem, coord, data, fieldname)
        character(len=*) :: filename, fieldname
        real(kind=CUSTOM_REAL), dimension(:,:) :: coord
        real(kind=CUSTOM_REAL), dimension(:) :: data
        integer, dimension(:,:) :: elem
        integer :: nelem, ncoord
        integer :: i

        nelem=size(elem(1,:))
        ncoord=size(coord(1,:))

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Mesh with Data'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
        do i=1,ncoord
            write(19,'(3e18.6)') coord(1,i),coord(2,i),coord(3,i)
        enddo
        write(19,*) ""

        write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*5
        do i=1,nelem
            write(19,'(5i12)') 4 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1,elem(4,i)-1
        end do
        write(19,*) ""

        write(19,'(a,i12)') "CELL_TYPES ",nelem
        write(19,*) (10,i=1,nelem)
        write(19,*) ""

        write(19,'(a,i12)') "CELL_DATA ",nelem
        write(19,'(a,a,a)') "SCALARS ", fieldname, " float"
        write(19,'(a)') "LOOKUP_TABLE default"
        do i = 1,nelem
            write(19,*) data(i)
        enddo
        write(19,*) ""
        close(19)
    end subroutine writeVtkTetraMeshRealdata

    subroutine writeVtkTriMeshRealdata(filename, elem, coord, data, fieldname)
        character(len=*) :: filename, fieldname
        real(kind=CUSTOM_REAL), dimension(:,:) :: coord
        real(kind=CUSTOM_REAL), dimension(:) :: data
        integer, dimension(:,:) :: elem
        integer :: nelem, ncoord
        integer :: i

        nelem=size(elem(1,:))
        ncoord=size(coord(1,:))

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Mesh with Data'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
        do i=1,ncoord
            write(19,'(3e18.6)') coord(1,i),coord(2,i),coord(3,i)
        enddo
        write(19,*) ""

        write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*4
        do i=1,nelem
            write(19,'(5i12)') 3 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1
        end do
        write(19,*) ""

        write(19,'(a,i12)') "CELL_TYPES ",nelem
        do i=1,nelem
            write(19,*) 5
        end do
        write(19,*) ""

        write(19,'(a,i12)') "CELL_DATA ",nelem
        write(19,'(a,a,a)') "SCALARS ", fieldname, " float"
        write(19,'(a)') "LOOKUP_TABLE default"
        do i = 1,nelem
            write(19,*) data(i)
        enddo
        write(19,*) ""
        close(19)
    end subroutine writeVtkTriMeshRealdata

    subroutine writeVtkTriMeshRealdata2D(filename, elem, coord, data, fieldname)
        character(len=*) :: filename, fieldname
        real(kind=CUSTOM_REAL), dimension(:,:) :: coord
        real(kind=CUSTOM_REAL), dimension(:) :: data
        integer, dimension(:,:) :: elem
        integer :: nelem, ncoord
        integer :: i

        nelem=size(elem(1,:))
        ncoord=size(coord(1,:))

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Mesh with Data'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(19, '(a,i12,a)') 'POINTS ', ncoord, ' float'
        do i=1,ncoord
            write(19,'(3e18.6)') coord(1,i),0.0,coord(2,i)
        enddo
        write(19,*) ""

        write(19,'(a,i12,i12)') "CELLS ",nelem,nelem*4
        do i=1,nelem
            write(19,'(5i12)') 3 ,elem(1,i)-1,elem(2,i)-1,elem(3,i)-1
        end do
        write(19,*) ""

        write(19,'(a,i12)') "CELL_TYPES ",nelem
        do i=1,nelem
            write(19,*) 5
        end do
        write(19,*) ""

        write(19,'(a,i12)') "CELL_DATA ",nelem
        write(19,'(a,a,a)') "SCALARS ", fieldname, " float"
        write(19,'(a)') "LOOKUP_TABLE default"
        do i = 1,nelem
            write(19,*) data(i)
        enddo
        write(19,*) ""
        close(19)
    end subroutine writeVtkTriMeshRealdata2D

    subroutine writeVtkNodes(filename, x,y,z)
        character(len=*) :: filename
        real(kind=CUSTOM_REAL), dimension(:) ::x,y,z
        integer :: n
        integer :: i

        n=size(x)

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Mesh with Data'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET POLYDATA'
        write(19, '(a,i12,a)') 'POINTS ', n, ' float'
        do i=1,n
            write(19,'(3e18.6)') x(i),y(i),z(i)
        enddo
        write(19,*) ""
        write(19,'(a,i12,i12)') "VERTICES ",n,n*4
        do i=1,n
            write(19,'(5i12)') 3 ,i-1,i-1,i-1
        end do
        write(19,*) ""
        close(19)
    end subroutine writeVtkNodes

    subroutine writeVtkNodesRealData(filename, x,y,z,data, fieldname)
        character(len=*) :: filename, fieldname
        real(kind=CUSTOM_REAL), dimension(:) ::x,y,z,data
        integer :: n
        integer :: i

        n=size(x)

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Mesh with Data'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET POLYDATA'
        write(19, '(a,i12,a)') 'POINTS ', n, ' float'
        do i=1,n
            write(19,'(3e18.6)') x(i),y(i),z(i)
        enddo
        write(19,*) ""
        write(19,'(a,i12,i12)') "VERTICES ",n,n*4
        do i=1,n
            write(19,'(5i12)') 3 ,i-1,i-1,i-1
        end do
        write(19,'(a,i12)') "POINT_DATA ",n
        write(19,'(a,a,a)') "SCALARS ", fieldname, " float"
        write(19,'(a)') "LOOKUP_TABLE default"
        do i = 1,n
            write(19,*) data(i)
        enddo
        write(19,*) ""
        close(19)
    end subroutine writeVtkNodesRealData

    subroutine writeVtkNodesRealDataBin(filename, x,y,z,data, fieldname)
        character(len=*) :: filename, fieldname
        real(kind=CUSTOM_REAL), dimension(:) ::x,y,z,data
        integer :: n
        integer :: i

        n=size(x)

        open(19,file=trim(filename),status='unknown',form='UNFORMATTED')
        write(19) '# vtk DataFile Version 3.1'
        write(19) 'Mesh with Data'
        write(19) 'ASCII'
        write(19) 'DATASET POLYDATA'
        write(19) 'POINTS ', n, ' float'
        do i=1,n
            write(19) x(i),y(i),z(i)
        enddo
        write(19) ""
        write(19) "VERTICES ",n,n*4
        do i=1,n
           write(19) 3 ,i-1,i-1,i-1
        end do
        write(19) "POINT_DATA ",n
        write(19) "SCALARS ", fieldname, " float"
        write(19) "LOOKUP_TABLE default"
        do i = 1,n
            write(19) data(i)
        enddo
        write(19) ""
        close(19)
    end subroutine writeVtkNodesRealDataBin

    subroutine writeVTKfractures(filename, points, fracture_size, nfracs)
        !input
        character(len=*), intent(in) :: filename
        integer :: nfracs
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: points
        integer, dimension(:,:), intent(in) :: fracture_size
        !local
        integer :: i, j, n
        character(len = 50) :: myfmt

        n = size(points(:,1))

        open(19,file=trim(filename),status='unknown')
        write(19,'(a)') '# vtk DataFile Version 3.1'
        write(19,'(a)') 'Line representation of the fracture'
        write(19,'(a)') 'ASCII'
        write(19,'(a)') 'DATASET POLYDATA'
        write(19,'(a,i12,a)') 'POINTS ', n, ' float'
        do i = 1, n
            write(19, '(3e18.6)') points(i,1),0.,points(i,2)
        enddo
        write(19, *) ""

        write(19,'(a, 2i12)') 'LINES ', nfracs, n+nfracs

        do i = 1, nfracs
            write(19, '(i12)', advance='no') fracture_size(i, 3)
            write(myfmt, '("(",I0,"(I12))")') fracture_size(i, 3)
            write(19, myfmt) (j-1, j=fracture_size(i,1), fracture_size(i,2))
        enddo
        write(19,*) ""
        close(19)
    end subroutine


end module vtkMod
