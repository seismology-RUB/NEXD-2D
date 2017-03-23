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
module dtMod
    ! calculate scaling for dt in tri
    use constantsMod
    
    implicit none

    contains

    subroutine scaleDT(x,z,elem,nelem,sDT)
        !Input
        real(kind=CUSTOM_REAL), dimension(:) :: x,z
        integer, dimension(:,:) :: elem
        integer :: nelem
        !Output
        real(kind=CUSTOM_REAL), dimension(nelem) :: sDT
        !Local
        real(kind=CUSTOM_REAL) :: len1,len2,len3,sper,A
        integer :: ie

        do ie=1,nelem
           len1=sqrt( (x(elem(1,ie)) - x(elem(2,ie)))**2 + (z(elem(1,ie)) - z(elem(2,ie)))**2)
           len2=sqrt( (x(elem(2,ie)) - x(elem(3,ie)))**2 + (z(elem(2,ie)) - z(elem(3,ie)))**2)
           len3=sqrt( (x(elem(3,ie)) - x(elem(1,ie)))**2 + (z(elem(3,ie)) - z(elem(1,ie)))**2)
           sper = (len1+len2+len3)/2.0
           a= sqrt(sper*(sper-len1)*(sper-len2)*(sper-len3))
           sDT(ie) = a/sper
        end do
  end subroutine scaleDT
end module dtMod
    
