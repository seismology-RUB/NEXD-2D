!-----------------------------------------------------------------------
!   Copyright 2013 Wolfgang Friederich (Ruhr-Universität Bochum, GER)
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
!-----------------------------------------------------------------------
!  Fast Fourier transform routines
!     schnelle fouriertransformation der 2**n komplexen werte (x,y)
!     is negativ - in den frequenzbereich / positiv - in den zeitbereich
!     normierung - abs(is) =1 - ohne / 2 - mit 1/ng / 3 - mit 1/sqrt(ng)
!    4 - ohne normierung,ohne kontrollausdruck
!-----------------------------------------------------------------------
module fourierTransform
    use mathConstants
    implicit none
    interface fastFourierTransform
        module procedure realFastFourierTransform
        module procedure doubleFastFourierTransform
    end interface
!
contains
    !----------------------------------------------------------------------
    !  single precision version (ssft.f)
    !
    subroutine realFastFourierTransform(x,y,n,is)
        real, dimension(:) :: x,y
        integer :: is,n,m
        integer, dimension(21) :: zh
        integer :: l,ng,nar,lar,larh,jr,ja,nr,jb,j,js,k,ny,nny
        real :: gn,alpha,piz,beta,excos,exsin,zx,zy
        !
        piz=2.*mc_pi
        !
        !  tabelle der zweierpotenzen
        !
        zh(1)=1
        do l=1,n
            zh(l+1)=2*zh(l)
        enddo
        ng=zh(n+1)
        gn=1./float(ng)
        !
        !  kernprogramm, dreifache schleife ueber schritt/index/teilserie
        !
        do m=1,n
            nar=zh(m)
            lar=ng/nar
            larh=lar/2
            alpha =  piz/float(isign(lar,is))
            do jr=1,larh
                beta=alpha*float(jr-1)
                excos = cos(beta)
                exsin = sin(beta)
                ja=jr-lar
                do nr=1,nar
                    ja=ja+lar
                    jb=ja+larh
                    zx = x(ja)-x(jb)
                    zy = y(ja)-y(jb)
                    x(ja) = x(ja)+x(jb)
                    y(ja) = y(ja)+y(jb)
                    x(jb) = zx*excos-zy*exsin
                    y(jb) = zx*exsin+zy*excos
                enddo
            enddo
        enddo
        !
        !     normierung
        !
        if (iabs(is) == 3) gn=sqrt(gn)
        if (iabs(is) == 2 .or. iabs(is) == 3) then
            y = y*gn
            x = x*gn
        endif
        !
        !     umordnung nach "bitreversed" indizes
        !
        do j=1,ng
            js=j-1
            k=1
            nny=n+1
            do ny=1,n
                nny=nny-1
                if (js.lt.zh(nny)) cycle
                js=js-zh(nny)
                k=k+zh(ny)
            enddo
            if (j-k < 0) then
                zx = x(j)
                zy = y(j)
                x(j) = x(k)
                y(j) = y(k)
                x(k) = zx
                y(k) = zy
            else
                cycle
            endif
        enddo
    end subroutine realFastFourierTransform
    !----------------------------------------------------------------------
    !  double precision version (sft.f)
    !
    subroutine doubleFastFourierTransform(x,y,n,is)
        double precision, dimension(:) :: x,y
        integer :: is,n,m
        integer, dimension(21) :: zh
        integer :: l,ng,nar,lar,larh,jr,ja,nr,jb,j,js,k,ny,nny
        double precision :: gn,alpha,piz,beta,excos,exsin,zx,zy
        !
        piz=2.*mc_pid
        !
        !  tabelle der zweierpotenzen
        !
        zh(1)=1
        do l=1,n
            zh(l+1)=2*zh(l)
        enddo
        ng=zh(n+1)
        gn=1./dble(ng)
        !
        !  kernprogramm, dreifache schleife ueber schritt/index/teilserie
        !
        do m=1,n
            nar=zh(m)
            lar=ng/nar
            larh=lar/2
            alpha =  piz/dble(isign(lar,is))
            do jr=1,larh
                beta=alpha*dble(jr-1)
                excos = dcos(beta)
                exsin = dsin(beta)
                ja=jr-lar
                do nr=1,nar
                    ja=ja+lar
                    jb=ja+larh
                    zx = x(ja)-x(jb)
                    zy = y(ja)-y(jb)
                    x(ja) = x(ja)+x(jb)
                    y(ja) = y(ja)+y(jb)
                    x(jb) = zx*excos-zy*exsin
                    y(jb) = zx*exsin+zy*excos
                enddo
            enddo
        enddo
        !
        !     normierung
        !
        if (iabs(is) == 3) gn=dsqrt(gn)
        if (iabs(is) == 2 .or. iabs(is) == 3) then
            y = y*gn
            x = x*gn
        endif
        !
        !     umordnung nach "bitreversed" indizes
        !
        do j=1,ng
            js=j-1
            k=1
            nny=n+1
            do ny=1,n
                nny=nny-1
                if (js.lt.zh(nny)) cycle
                js=js-zh(nny)
                k=k+zh(ny)
            enddo
            if (j-k < 0) then
                zx = x(j)
                zy = y(j)
                x(j) = x(k)
                y(j) = y(k)
                x(k) = zx
                y(k) = zy
            else
                cycle
            endif
        enddo
    end subroutine doubleFastFourierTransform
!
end module fourierTransform
