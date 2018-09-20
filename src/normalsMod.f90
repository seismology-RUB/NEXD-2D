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
module normalsMod
    ! module to compute the outwart pointing normals at element faces and the surface jacobians
    use constantsMod

    implicit none

    contains

    subroutine normals2d(nx,nz,sj,Dr,Ds,x,z,fmask)
        !input
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,z
        real(kind=CUSTOM_REAL), dimension(:,:):: Dr,Ds
        integer, dimension(:,:), intent(in) :: fmask
        !output
        real(kind=CUSTOM_REAL), dimension(3*NGLL), intent(out) :: nx,nz,sj
        !local
        real(kind=CUSTOM_REAL), dimension(Np) :: xr,zr,xs,zs
        real(kind=CUSTOM_REAL), dimension(3*NGLL) :: fxr,fzr,fxs,fzs
        integer, dimension(3*NGLL) :: fmaskv
        integer, dimension(NGLL) :: fid1,fid2,fid3
        integer :: i,k,l

        fid1 = 0
        fid2 = 0
        fid3 = 0
        k=1
        do l=1,3
            do i=1,NGLL
                fmaskv(k)=fmask(i,l)
                if (l==1) then
                    fid1(i)=k
                else if (l==2) then
                    fid2(i)=k
                else
                    fid3(i)=k
                end if
            k=k+1
            end do
        end do

        xr = matmul(Dr,x)
        zr = matmul(Dr,z)
        xs = matmul(Ds,x)
        zs = matmul(Ds,z)

        ! interpolate geometric factors to face nodes
        fxr = xr(fmaskv)
        fxs = xs(fmaskv)
        fzr = zr(fmaskv)
        fzs = zs(fmaskv)

        ! build normals
        ! face 1
        nx(fid1) = fzr(fid1)
        nz(fid1) = -fxr(fid1)
        ! face 2
        nx(fid2) = fzs(fid2)-fzr(fid2)
        nz(fid2) = -fxs(fid2)+fxr(fid2)
        ! face 3
        nx(fid3) = -fzs(fid3)
        nz(fid3) = fxs(fid3)

        ! normalize
        sj = sqrt(nx*nx+nz*nz)
        nx = nx/sJ
        nz = nz/sJ
    end subroutine normals2d
end module normalsMod
