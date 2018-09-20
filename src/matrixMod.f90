!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2018 Andre Lamert (Ruhr-Universität Bochum, GER)
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
        if (CUSTOM_REAL==SIZE_REAL) then
            call sgetrf(N,N,V,N,ipvt,ierr)
        else
            call dgetrf(N,N,V,N,ipvt,ierr)
        endif
        if (ierr/=0) then
            write(err_str,*) "Error LU vdm2d ",ierr
            call add(errmsg, 2, trim(err_str), myname)
            if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif
        endif

        ! ivert Pr
        if (CUSTOM_REAL==SIZE_REAL) then
            call sgetri(N,V,N,ipvt,work,iwork,ierr)
        else
            call dgetri(N,V,N,ipvt,work,iwork,ierr)
        endif
        if (ierr/=0) then
            write(err_str,*) "Error invert vdm2d ",ierr
            call add(errmsg, 2, trim(err_str), myname)
            if (.level.errmsg == 2) then; call print(errmsg, .true.); stop; endif
        endif
    end subroutine invert

end module matrixMod
