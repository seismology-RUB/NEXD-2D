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
        real(kind=CUSTOM_REAL), dimension(size(u)), intent(out) :: divu
        real(kind=CUSTOM_REAL), dimension(size(u)) :: ur,us,vr,vs

        ur = matmul(Dr,u)
        us = matmul(Ds,u)
        vr = matmul(Dr,v)
        vs = matmul(Ds,v)

        divu = rx*ur + sx*us + rz*vr + sz*vs
    end subroutine div2d
end module derMod
