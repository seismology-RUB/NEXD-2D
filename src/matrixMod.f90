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
module matrixMod
    ! some service routines for matrices
    use constantsMod
    use errorMessage
    implicit none

    contains

    subroutine invert(V, errmsg)
        !invert matrix
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: V
        real(kind=CUSTOM_REAL), dimension(2*size(v(:,1))) :: work
        integer, dimension(size(v(:,1))) :: ipvt
        integer :: ierr, iwork
        integer :: N
        character(len=6) :: myname = "invert"
        character(len=180) :: err_str
        call addTrace(errmsg, myname)

        N=size(v(:,1))
        iwork=2*N
        !LU trafo
        call sgetrf(N,N,V,N,ipvt,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error LU vdm2d ",ierr
            call add(errmsg, 2, trim(err_str), myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif

        ! ivert Pr
        call sgetri(N,V,N,ipvt,work,iwork,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error invert vdm2d ",ierr
            call add(errmsg, 2, trim(err_str), myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif
    end subroutine invert

    subroutine matrixLU(V, errmsg)
        !invert matrix
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: V
        integer, dimension(size(v(:,1))) :: ipvt
        integer :: ierr, iwork
        integer :: N
        character(len=6) :: myname = "matrixLU"
        character(len=180) :: err_str

        call addTrace(errmsg, myname)

        N=size(v(:,1))
        iwork=2*N
        !LU trafo
        call sgetrf(N,N,V,N,ipvt,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error LU vdm2d ",ierr
            call add(errmsg, 2, trim(err_str), myname)
            if (.level.errmsg == 2) then; call print(errmsg); stop; endif
        endif
    end subroutine matrixLU

    function mxv(M,v,N)
        ! matrix times vector
        integer :: N
        real(kind=CUSTOM_REAL), dimension(N,N) :: M
        real(kind=CUSTOM_REAL), dimension(N) :: v, mxv
        real(kind=CUSTOM_REAL) :: temp
        integer :: i,j

        do i=1,N
            temp=0
            do j=1,N
                temp=temp+M(i,j)*v(j)
            end do
            mxv(i)=temp
        end do
  end function mxv

end module matrixMod
