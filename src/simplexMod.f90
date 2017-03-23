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
module simplexMod
    ! Evaluate 2D orthonormal polynomial on simplex at (a,b) oft order (i,j)
    use constantsMod
    use jacobiMod

    implicit none

    contains

    subroutine simplex2DP(Pr,a,b,i,j)
        real(kind=CUSTOM_REAL), dimension(:), intent(out) :: Pr
        real(kind=CUSTOM_REAL), dimension(:), intent(in):: a,b
        real(kind=CUSTOM_REAL), dimension(size(a)) :: h1,h2
        integer, intent(in) :: i,j

        call jacobiP(h1,a,0.0,0.0,i)
        call jacobiP(h2,b,2.0*i+1.0,0.0,j)

        Pr(:) = sqrt(2.0)*h1*h2*(1.0-b)**i
    end subroutine simplex2DP

    subroutine gradSimplex2DP(dPdr,dPds,a,b,id,jd)
        ! computes the differentiation matrices Dr and Ds
        real(kind=CUSTOM_REAL), dimension(:), intent(in) ::a,b
        integer, intent(in) :: id,jd
        real(kind=CUSTOM_REAL), dimension(size(a)), intent(out) :: dPdr,dPds
        real(kind=CUSTOM_REAL), dimension(size(a)) :: fa,dfa
        real(kind=CUSTOM_REAL), dimension(size(b)) :: gb,dgb
        real(kind=CUSTOM_REAL), dimension(size(a)) :: temp

        call jacobiP(fa,a,0.0,0.0,id)
        call gradJacobiP(dfa,a,0.0,0.0,id)
        call jacobiP(gb,b,2.0*id+1.0,0.0,jd)
        call gradJacobiP(dgb,b,2.0*id+1.0,0.0,jd)

        ! r-deriverate
        !d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da

        dPdr(:)=dfa(:)*gb(:)
        if (id>0) then
            dPdr(:)=dPdr(:)*((0.5*(1.0-b(:)))**(id-1.0))
        end if

        ! s-deriverate
        !d/ds = ((1+a)/2)/((1-b)/s) d/da + d/db

        dPds(:) = dfa(:)*(gb(:)*(0.5*(1.0+a(:))))
        if (id>0) then
            dPds(:)=dPds(:)*((0.5*(1.0-b(:)))**(id-1.0))
        end if
    
        temp(:) = dgb(:)*((0.5*(1.0-b(:)))**id)
        if (id>0) then
            temp(:)=temp(:)-0.5*id*gb(:)*((0.5*(1.0-b(:)))**(id-1.0))
        end if
        dPds(:) = dPds(:)+fa(:)*temp(:)
    
        !normalize
        dPdr(:) = 2.0**(id+0.5)*dPdr(:)
        dPds(:) = 2.0**(id+0.5)*dPds(:)
    end subroutine gradSimplex2DP
end module simplexMod
