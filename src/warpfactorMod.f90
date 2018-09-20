!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
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
module warpfactorMod
    ! module to calculate the warp factors for interpolation in the triangle with quasi gll points
    use constantsMod

    implicit none

    contains

    function warpfactor(xnodes,xout)
        integer :: p
        real(kind=CUSTOM_REAL), dimension(NGLL) :: xnodes,xeq
        real(kind=CUSTOM_REAL), dimension(Np) :: warpfactor,d, warp ,xout
        real(kind=CUSTOM_REAL) :: dx,prod,prod1,prod2,e,sf
        logical :: zerof
        integer :: i,j,k, zeroff

        p=NGLL-1

        !create equispaced vector in same ordering as xnodes
        dx=2./(NGLL-1)
        j=NGLL
        do i=1,NGLL
            xeq(i)=((j-1.0)*dx)-1.0
            j=j-1
        end do
        d(:)=0.
        do k=1,Np
            do i=1,NGLL
                e=(xnodes(i)-xeq(i))
                prod1=1.
                prod2=1.
                do j=1,NGLL
                    if (i/=j) then
                        prod1=prod1*(xout(k)-xeq(j))
                        prod2=prod2*(xeq(i)-xeq(j))
                    end if
                end do
                prod=e*(prod1/prod2)
                d(k)=d(k)+prod
            end do
            warp(k)=d(k)
            ! scale warp factor
            zerof = abs(xout(k)) < (1.0-EPS)
            if (zerof) then
                zeroff=1
            else
                zeroff=0
            end if
            sf = 1.0-(zeroff*xout(k))**2
            warpfactor(k)=warp(k)/sf + warp(k)*(zeroff-1.0)
        end do
    end function warpfactor
end module warpfactorMod
