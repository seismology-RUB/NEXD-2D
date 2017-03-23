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
module liftMod
    ! module to compute the lifting from surf to vol for dg
    use constantsMod
    use vandermondeMod
    use errorMessage
    implicit none

    contains

    subroutine lift2d(lift,fmask,v,r,s,errmsg)
        type(error_message) :: errmsg
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
        integer, dimension(NGLL,3), intent(in) :: fmask
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL), intent(out) :: lift
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL) :: temp
        real(kind=CUSTOM_REAL), dimension(Np,3*NGLL) :: emat
        real(kind=CUSTOM_REAL), dimension(NGLL) :: faceR, faceS
        real(kind=CUSTOM_REAL), dimension(NGLL,NGLL) :: v1d, masse1,masse2,masse3
        real(kind=CUSTOM_REAL), dimension(2*NGLL) :: work
        integer, dimension(NGLL) :: ipvt
        integer :: ierr, iwork
        integer :: i,j,k


        character(len=6) :: myname = "lift2d"
        character(len=180) :: err_str

        call addTrace(errmsg, myname)

        emat = 0.0
        lift = 0.0

        iwork=2*NGLL
        ! masse1
        faceR(:) = r(fmask(:,1))
        call vdm1d(v1d,faceR)

        masse1=matmul(v1d,transpose(v1d))
        ! LU
        call sgetrf(NGLL,NGLL,masse1,NGLL,ipvt,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error LU masse1 ",ierr
            call add(errmsg, 2, err_str, myname)
        end if

        ! ivert Pr
        call sgetri(NGLL,masse1,NGLL,ipvt,work,iwork,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error invert masse1 ",ierr
            call add(errmsg, 2, err_str, myname)
        endif

        do i=1,NGLL
            do j=1,NGLL
                emat(fmask(i,1),j)=masse1(i,j)
            end do
        end do

        ! masse2
        faceR(:) = r(fmask(:,2))
        call vdm1d(v1d,faceR)

        masse2=matmul(v1d,transpose(v1d))
        ! LU
        call sgetrf(NGLL,NGLL,masse2,NGLL,ipvt,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error LU masse2 ",ierr
            call add(errmsg, 2, err_str, myname)
        end if

        ! ivert Pr
        call sgetri(NGLL,masse2,NGLL,ipvt,work,iwork,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error invert masse2 ",ierr
            call add(errmsg, 2, err_str, myname)
        endif

        do i=1,NGLL
            k=1
            do j=NGLL+1,2*NGLL
                emat(fmask(i,2),j)=masse2(i,k)
                k=k+1
            end do
        end do

        ! masse3
        faceS(:) = s(fmask(:,3))
        call vdm1d(v1d,faceS)

        masse3=matmul(v1d,transpose(v1d))

        ! LU
        call sgetrf(NGLL,NGLL,masse3,NGLL,ipvt,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error LU masse3 ",ierr
            call add(errmsg, 2, err_str, myname)
        end if

        ! ivert Pr
        call sgetri(NGLL,masse3,NGLL,ipvt,work,iwork,ierr)
        if (ierr/=0) then
            write(err_str,*) "Error invert masse3 ",ierr
            call add(errmsg, 2, err_str, myname)
        endif

        do i=1,NGLL
            k=1
            do j=2*NGLL+1,3*NGLL
                emat(fmask(i,3),j)=masse3(i,k)
                k=k+1
            end do
        end do

        ! build lift matrix
        temp = matmul(transpose(v),emat)
        lift= matmul(v,temp)
    end subroutine lift2d
end module liftMod
