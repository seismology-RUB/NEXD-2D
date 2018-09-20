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
module derMod
    ! compute grad, div curl
    use constantsMod
    implicit none

    contains

    subroutine grad2d(ux,uz,u,Dr,Ds,rx,sx,rz,sz)
        ! compute gradient in 2d
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: rx,sx,rz,sz
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: u
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds
        real(kind=CUSTOM_REAL), dimension(size(u)), intent(out) :: ux,uz
        real(kind=CUSTOM_REAL), dimension(size(u)) :: ur,us

        ur = matmul(Dr,u)
        us = matmul(Ds,u)
        ux = rx*ur + sx*us
        uz = rz*ur + sz*us
    end subroutine grad2d

    subroutine div2d(divu,u,v,Dr,Ds,rx,sx,rz,sz)
        !compute divergence for vertorfield (u,v)
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: rx,sx,rz,sz
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: u,v
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds
        real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(out) :: divu
        real(kind=CUSTOM_REAL), dimension(size(u)) :: ur,us,vr,vs

        allocate(divu(size(u)))

        ur = matmul(Dr,u)
        us = matmul(Ds,u)
        vr = matmul(Dr,v)
        vs = matmul(Ds,v)

        divu = rx*ur + sx*us + rz*vr + sz*vs
    end subroutine div2d

    subroutine curl2d(ux,uz,rx,sx,rz,sz,Dr,Ds,curlu)
        !compute 2D curl operator in (x,z)-plane
        !Input
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: rx,sx,rz,sz
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: ux, uz
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds
        !output
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable, intent(out) :: curlu
        !local
        real(kind=CUSTOM_REAL), dimension(size(ux)) :: uxr, uxs, uzr, uzs

        allocate(curlu(size(ux), 3))

        uxr = matmul(Dr,ux)
        uxs = matmul(Ds,ux)
        uzr = matmul(Dr,uz)
        uzs = matmul(Ds,uz)

        curlu = 0
        curlu(:,2) = rz*uxr + sz*uxs - rx*uzr - sx*uzs

    end subroutine curl2d


end module derMod
