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
module dmatricesMod
    ! module to calculate the differentiation matrices on the simplex at (r,s)
    use constantsMod
    use vandermondeMod
    use errorMessage

    implicit none

    contains

    subroutine dmatrices2d(dr,ds,r,s,v, errmsg)
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
        real(kind=CUSTOM_REAL), dimension(size(v(:,1)),Np) :: dr,ds,w,vr,vs
        character(len=11) :: myname = "dmatrices2d"

        call addTrace(errmsg, myname)

        call gradVdm2D(vr,vs,r,s)
        call invVdm2D(v,w,0, errmsg)
    
        dr=matmul(vr,w)
        ds=matmul(vs,w)
    end subroutine dmatrices2d
end module dmatricesMod
