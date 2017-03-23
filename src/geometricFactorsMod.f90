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
module geometricFactorsMod
    ! module to compute the metric elements for local mappings of the element
    use constantsMod
    implicit none

    contains

    subroutine geometricFactors2d(rx,sx,ry,sy,J,x,y,Dr,Ds)
        !Input
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds
        !Output
        real(kind=CUSTOM_REAL), dimension(size(x)) :: rx,sx,ry,sy,J
        !Local
        character(len=18) :: myname = "geometricFactors2D"
        real(kind=CUSTOM_REAL), dimension(size(x)) :: xr,xs,yr,ys

        xr = matmul(Dr,x)
        xs = matmul(Ds,x)
        yr = matmul(Dr,y)
        ys = matmul(Ds,y)

        J = -xs*yr + xr*ys

        rx = ys/J
        sx = -yr/J
        ry = -xs/J
        sy = xr/J
    end subroutine geometricFactors2d
end module geometricFactorsMod
