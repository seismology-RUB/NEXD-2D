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
module externalModelMod
    implicit none

    contains

    subroutine define_external_model(x,y,iflag_element,rho,vp,vs,nx,nz,dx,dz,xmin,zmin,ipolg,xg,zg,rhog,vpg,vsg )

        !  include "constants.h"

        ! user can modify this routine to assign any different external Earth model (rho, vp, vs)
        ! based on the x and y coordinates of that grid point and the flag of the region it belongs to

        !  integer, intent(in) :: iflag_element,myrank
        integer :: iflag_element

        real, intent(in) :: x,y

        real, intent(out) :: rho,vp,vs
        integer :: cnt,is,ks,nx,nz,i,k
        integer, dimension(nx,nz) :: ipolg
        real, dimension(nx,nz) :: rhog,vpg,vsg
        real, dimension(nx) :: xg
        real, dimension(nz) :: zg
        real, dimension(4) :: xc,yc,xi
        integer, dimension(4) :: isinp,ksinp
        real :: xmin,zmin,dx,dz,bx,bz,d,xs,zs

        xs = x*1.e-3
        zs = -y*1.e-3

        ! interpolate values
        is = int((xs-xmin)/dx)+1
        ks = int((zs-zmin)/dz)+1

        if (iflag_element>7) then
            rho = rhog(is,ks)
            vp  = vpg(is,ks)
            vs  = vsg(is,ks)
        end if

        ! count how many of the four grid points have ipolg = iflag_element
        ! and store their indices
        cnt = 0
        do k = ks,min(ks+1,nz)
            do i = is,min(is+1,nx)
                if (ipolg(i,k) == iflag_element) then
                    cnt = cnt+1
                    isinp(cnt) = i; ksinp(cnt) = k
                end if
            enddo
        enddo

        !  take value at this single point
        if (cnt == 1) then
            xi(1) = 1.

        !  project point on line connecting grid points and do linear interpolation
        else if (cnt == 2) then
            xc(1:2) = xg(isinp(1:2)); yc(1:2) = zg(ksinp(1:2))
            d = sqrt((xc(2)-xc(1))**2+(yc(2)-yc(1))**2)          ! projected distance from P1
            xi(2) = (xs-xc(1))/d*(xc(2)-xc(1))/d+(zs-yc(1))/d*(yc(2)-yc(1))/d
            xi(1) = 1.-xi(2)
        !  do triangle interolation using barycentric coordinates
        else if (cnt == 3) then
            xc(1:3) = xg(isinp(1:3)); yc(1:3) = zg(ksinp(1:3))
            d = (yc(2)-yc(3))*(xc(1)-xc(3))+(xc(3)-xc(2))*(yc(1)-yc(3))
            xi(1) = ((yc(2)-yc(3))*(xs-xc(3))+(xc(3)-xc(2))*(zs-yc(3)))/d
            xi(2) = ((yc(3)-yc(1))*(xs-xc(3))+(xc(1)-xc(3))*(zs-yc(3)))/d
            xi(3) = 1.-xi(1)-xi(2)
        !  do bilinear interpolation
        else if (cnt == 4) then
            bx = (xs-xg(is))/dx; bz = (zs-zg(ks))/dz
            xi(1) = (1.-bx)*(1.-bz)
            xi(2) = bx*(1.-bz)
            xi(3) = (1.-bx)*bz
            xi(4) = bx*bz
        endif

        !  Interpolation
        if (cnt >0) then
            rho = 0
            vp = 0
            vs = 0
            do i = 1,cnt
                rho = rho +xi(i)*rhog(isinp(i),ksinp(i))
                vp  = vp + xi(i)*vpg(isinp(i),ksinp(i))
                vs  = vs + xi(i)*vsg(isinp(i),ksinp(i))
            enddo
        end if

        if (vp<=vs) write(*,*) "ERROR, vp smaller than vs", vp, vs, cnt, is, ks, isinp(1:cnt), ksinp(1:cnt), &
           (vpg(isinp(i),ksinp(i)),i=1,cnt),(vsg(isinp(i),ksinp(i)),i=1,cnt), xi(1:cnt)
    end subroutine define_external_model
end module externalModelMod
