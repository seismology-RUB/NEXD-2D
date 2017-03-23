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
module partMod
    use constantsMod
    ! module to deal with the partitioning of the mesh
    ! we will use metis at the moment

    !This module is currently not in use TM TM
    !essentially this module covers the code used in MeshMod after if (nproc > 1) then (currently line 754 - 882
    contains

    ! built partition of the mesh
    subroutine partMesh(nproc, nelem, elem, neighbor)
        integer :: nelem
        integer, dimension(:,:) :: elem
        integer, dimension(:,:) :: neighbor
        integer, pointer :: nproc
        integer, dimension(:), allocatable :: elmnts
        integer, dimension(:), allocatable  :: xadj
        integer, dimension(:), allocatable  :: adjncy
        integer :: nb_neigh
        integer , dimension(3) :: nb_temp
        integer, dimension(0:4) :: options
        integer :: edgecut
        integer, dimension(:),allocatable :: part
        integer, dimension(:), allocatable :: glob2loc_elmnts
        integer :: num_glob, num_part
        integer, dimension(:), allocatable :: num_loc
        integer, dimension(:), allocatable :: nodes_elmnts, nnodes_elmnts

        integer, dimension(:), allocatable :: glob2loc_nodes_nparts
        integer, dimension(:), allocatable :: glob2loc_nodes_parts
        integer, dimension(:), allocatable :: glob2loc_nodes
        integer :: size_glob2loc_nodes
        integer, dimension(:), allocatable :: part_nodes, num_parts

    
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
        !part(:) => null()
        options(:)=0
        write(*,*) "METIS"
        call METIS_PartGraphRecursive(this%nelem, xadj, adjncy, 0, 0, 0, 1, nproc,options, edgecut, part)
        write(*,*) "END METIS"
        write(*,*) part

        ! create glob2loc_elem
        ! inspired by the specfem partition
        ! be careful, metis gives parts begining from 0 or from 1 depending on the version, so i compile a local version, here 4.0.3
        allocate(glob2loc_elmnts(this%nelem))
        glob2loc_elmnts(:) = 0

        allocate(num_loc(nproc))
        num_loc(:) = 1

        ! build a vector with the local numbering of the elements
        do num_glob=1,this%nelem
            num_part=part(num_glob)
            glob2loc_elmnts(num_glob) = num_loc(num_part)
            num_loc(num_part) = num_loc(num_part) + 1
        end do

        ! create local node numbering
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
        write(*,*)  "nnodes_elmnts", nnodes_elmnts
        write(*,*)  "nodes_elmnts", nodes_elmnts

        ! now create the local node numbering
        allocate(glob2loc_nodes_nparts(this%ncoord))
        allocate(part_nodes(nproc), num_parts(nproc))

        size_glob2loc_nodes = 0
        part_nodes(:) = 0

        ! go through all coordinates
        do in=1,this%ncoord
            glob2loc_nodes_nparts(in) = size_glob2loc_nodes
            do ie = 1, nnodes_elmnts(in)
                ! nodes_elmnts((ie)+nsize*(in-1)) gibt eins der maximal _nsize_ elemente zum knoten _in_
                ! gibt an in welchen partitionen der knoten alles liegt
                part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
            end do
            do num_part=1,nproc ! gehe durch die partitionen
                if (part_nodes(num_part) == 1) then
                    ! get number of nodes in the global array, there might be some nodes at the interfaces doubble
                    size_glob2loc_nodes = size_glob2loc_nodes +1
                    part_nodes(num_part) = 0
                end if
            end do
            glob2loc_nodes_nparts(this%ncoord) = size_glob2loc_nodes
        end do

        allocate(glob2loc_nodes_parts(glob2loc_nodes_nparts(this%ncoord)))
        allocate(glob2loc_nodes(glob2loc_nodes_nparts(this%ncoord)))

        glob2loc_nodes(:) = 0
        part_nodes(:) = 0
        num_parts(:)=1
        size_glob2loc_nodes = 1

        do in=1,this%ncoord
            do ie = 1, nnodes_elmnts(in)
                ! nodes_elmnts((ie)+nsize*(in-1)) gibt eins der maximal _nsize_ elemente zum knoten _in_
                ! gibt an in welchen partitionen der knoten alles liegt
                part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
            end do
            do num_part=1,nproc
                if (part_nodes(num_part) == 1) then
                    glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
                    glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
                    size_glob2loc_nodes = size_glob2loc_nodes +1
                    num_parts(num_part) = num_parts(num_part) +1
                    part_nodes(num_part) = 0
                end if
            end do
        end do
        write(*,*) "glob2loc_nodes_parts",glob2loc_nodes_parts
        write(*,*) "glob2loc_nodes",glob2loc_nodes
    end subroutine partMesh
end module partMod
