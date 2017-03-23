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
module pmlMod
    use constantsMod

    implicit none

    contains

    subroutine dampingRegular(ddx,ddz,alphax,alphaz,kx,kz,vx,vz,vp,xm,zm,pml_delta,rc,amax,kmax,corner)
        logical :: corner
        real(kind=CUSTOM_REAL), dimension(Np) :: vx,vz
        real(kind=CUSTOM_REAL), dimension(Np) :: ddx,ddz,alphax,alphaz,kx,kz
        real(kind=CUSTOM_REAL), dimension(Np) :: ar,br
        real(kind=CUSTOM_REAL) :: pml_delta,vp,rc,amax,kmax, xm, zm
        integer :: j

        if (corner) then
            do j=1,Np
                if (abs(vz(j)-zm) > (abs(vx(j)-xm))) then
                    ar(j)=sqrt((xm+pml_delta)**2)
                    br(j)=sqrt(vx(j)**2)
                    br(j)=abs(ar(j)-br(j))
                else
                    ar(j)=sqrt((zm+pml_delta)**2)
                    br(j)=sqrt(vz(j)**2)
                    br(j)=abs(ar(j)-br(j))
                end if
            end do
        else
            ar=sqrt((xm+pml_delta)**2)
            br=sqrt(vx**2)
            br=abs(ar-br)
        end if

        ddx = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
        ddz = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
        alphax =amax*(1-br/pml_delta)
        alphaz =amax*(1-br/pml_delta)
        kx=1+(kmax-1)*(br/pml_delta)**2
        kz=1+(kmax-1)*(br/pml_delta)**2
    end subroutine dampingRegular

    subroutine dampingProfile(ddx,ddz,alphax,alphaz,kx,kz,vx,vz,vp,xmin,xmax,zmin,zmax,pml_delta,rc,amax,kmax,iv,ipml)
        real(kind=CUSTOM_REAL), dimension(:) :: vx,vz
        real(kind=CUSTOM_REAL), dimension(:) :: ddx,ddz,alphax,alphaz,kx,kz
        real(kind=CUSTOM_REAL), dimension(Np) :: ar,br
        real(kind=CUSTOM_REAL) :: pml_delta,vp,rc,amax,kmax, xmin, xmax, zmin, zmax
        integer :: j,ipml
        integer, dimension(NP) :: iv

        if ( ipml == 1) then !left
            ar=sqrt((xmin+pml_delta)**2)
            br=sqrt(vx(iv)**2)
            br=abs(ar-br)
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) =amax*(1-br/pml_delta)
            alphaz(iv) =amax*(1-br/pml_delta)
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( ipml == 2 ) then !bottom
            ar=sqrt((zmin+pml_delta)**2)
            br=sqrt(vz(iv)**2)
            br=abs(ar-br)
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            alphaz(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( ipml == 3 ) then !right
            ar=sqrt((xmax-pml_delta)**2)
            br=sqrt(vx(iv)**2)
            br=abs(ar-br)
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            alphaz(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( ipml == 4 ) then !top
            ar=sqrt((zmax-pml_delta)**2)
            br=sqrt(vz(iv)**2)
            br=abs(ar-br)
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            alphaz(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*( br/pml_delta )**2
            kz(iv)=1+(kmax-1)*( br/pml_delta )**2
        else if ( ipml == 5) then !bottom left
            do j=1,Np
                if (abs(vz(iv(j))-zmin) > (abs(vx(iv(j))-xmin))) then
                    ar(j)=sqrt((xmin+pml_delta)**2)
                    br(j)=sqrt(vx(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                else
                    ar(j)=sqrt((zmin+pml_delta)**2)
                    br(j)=sqrt(vz(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                end if
            end do
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) =amax*(1-br/pml_delta)
            alphaz(iv) =amax*(1-br/pml_delta)
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( ipml == 6 ) then !bottom right
            do j=1,Np
                if (abs(vz(iv(j))-zmin) > (abs(vx(iv(j))-xmax))) then
                    ar(j)=sqrt((xmax-pml_delta)**2)
                    br(j)=sqrt(vx(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                else
                    ar(j)=sqrt((zmin+pml_delta)**2)
                    br(j)=sqrt(vz(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                end if
            end do
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            alphaz(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( ipml == 7 ) then !top right
            do j=1,Np
                if (abs(vz(iv(j))-zmax) > (abs(vx(iv(j))-xmax))) then
                    ar(j)=sqrt((xmax-pml_delta)**2)
                    br(j)=sqrt(vx(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                else
                    ar(j)=sqrt((zmax-pml_delta)**2)
                    br(j)=sqrt(vz(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                end if
            end do
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            alphaz(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*(br/pml_delta)**2
            kz(iv)=1+(kmax-1)*(br/pml_delta)**2
        else if ( ipml == 8 ) then !top left
            do j=1,Np
                if (abs(vz(iv(j))-zmax) > (abs(vx(iv(j))-xmin))) then
                    ar(j)=sqrt((xmin+pml_delta)**2)
                    br(j)=sqrt(vx(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                else
                    ar(j)=sqrt((zmax-pml_delta)**2)
                    br(j)=sqrt(vz(iv(j))**2)
                    br(j)=abs(ar(j)-br(j))
                end if
            end do
            ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(rc)*( br/pml_delta )**2
            alphax(iv) = amax*(1- br/pml_delta )
            alphaz(iv) = amax*(1- br/pml_delta )
            kx(iv)=1+(kmax-1)*( br/pml_delta )**2
            kz(iv)=1+(kmax-1)*( br/pml_delta )**2
        end if
    end subroutine dampingProfile
end module pmlMod
