!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2018 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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
module vandermondeMod
    ! module to create the vandermonde matrix
    use constantsMod
    use simplexMod
    use triTrafoMod
    use jacobiMod
    use errorMessage
    use matrixMod

    implicit none

    contains

    subroutine vdm1d(v,r)
        ! initialize th 1D vandermonde matrix
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r
        real(kind=CUSTOM_REAL), dimension(size(r),NGLL), intent(out) :: v
        integer :: i

        do i=1,NGLL
            call jacobiP(v(:,i),r(:),zero,zero,i-1)
        end do
    end subroutine vdm1d


    subroutine vdm2d(v,r,s)
        ! initialize the 2D vandermonde matrix
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(size(r),Np), intent(out) :: v
        real(kind=CUSTOM_REAL), dimension(size(r)) :: a,b
        integer :: i,j,sk
        !
        call rstoab(r,s,a,b)

        sk=1
        do i=0,NGLL-1
            do j=0,NGLL-1-i
                call simplex2DP(v(:,sk),a,b,i,j)
                sk=sk+1
            end do
        end do
    end subroutine vdm2d

    subroutine invVdm2D(v,w,trans, errmsg)
        ! calulates v^-1' if trans==0 -> inverse if trans == 1 => inverse transponierte
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
        real(kind=CUSTOM_REAL), dimension(size(v(:,1)),Np) :: w
        integer, intent(in) :: trans

        character(len=8) :: myname = "invVdm2D"

        call addTrace(errmsg, myname)
        w=v
        call invert(w, errmsg)

        if (trans==1) then
           w=transpose(w)
        end if
    end subroutine invVdm2D

    subroutine gradVdm2D(v2dr,v2ds,r,s)
        ! build the gradient of the modal basis (i,j) at (r,s) in referenze tri
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(size(r),Np), intent(out) :: v2dr,v2ds
        real(kind=CUSTOM_REAL), dimension(size(r)) :: a,b
        integer :: i,j,sk
        !
        call rstoab(r,s,a,b)

        sk=1
        do i=0,NGLL-1
            do j=0,NGLL-1-i
                call gradSimplex2DP(v2dr(:,sk),v2ds(:,sk),a,b,i,j)
                sk=sk+1
            end do
        end do
    end subroutine gradVdm2D
end module vandermondeMod
