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
module vandermondeMod
    ! module to create the vandermonde matrix
    use constantsMod
    use simplexMod
    use triTrafoMod
    use jacobiMod
    use errorMessage

    implicit none

    contains

    subroutine vdm1d(v,r)
        ! initialize th 1D vandermonde matrix
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r
        real(kind=CUSTOM_REAL), dimension(size(r),NGLL), intent(out) :: v
        integer :: i

        do i=1,NGLL
            call jacobiP(v(:,i),r(:),0.0,0.0,i-1)
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
        real(kind=CUSTOM_REAL), dimension(2*Np) :: work
        integer, dimension(Np) :: ipvt
        integer :: ierr, iwork
        character(len=8) :: myname = "invVdm2D"
        character(len=50) :: errstr

        call addTrace(errmsg, myname)

        iwork=2*Np
        w=v
        !LU trafo
        call sgetrf(Np,Np,w,Np,ipvt,ierr)
        if (ierr/=0) then
            write(errstr,*) "Error LU vdm2d ",ierr
            call add(errmsg, 2, errstr, myname)
            call print(errmsg)
            stop
        end if

        ! ivert Pr
        call sgetri(Np,w,Np,ipvt,work,iwork,ierr)
        if (ierr/=0) then
            write(errstr,*) "Error invert vdm2d ",ierr
            call add(errmsg, 2, errstr, myname)
            call print(errmsg)
            stop
        end if
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
