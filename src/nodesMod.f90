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
module nodesMod
    ! module to create warp and blend nodes in equi tri
    use constantsMod
    use gllMod
    use warpfactorMod
    implicit none

    contains

    subroutine nodes(alpha,x,y)
        real(kind=CUSTOM_REAL), intent(in) :: alpha
        real(kind=CUSTOM_REAL), dimension(:), intent(out) :: x,y
        real(kind=CUSTOM_REAL), dimension(Np) :: L1,L2,L3, blend1, blend2, blend3
        real(kind=CUSTOM_REAL), dimension(Np) :: warpfactor1,warpfactor2,warpfactor3
        real(kind=CUSTOM_REAL), dimension(Np) :: warp1,warp2,warp3
        real(kind=CUSTOM_REAL), dimension(NGLL) :: xgll
        integer :: p,sk,i,j

        p=NGLL-1

        ! get GLL nodes
        xgll=getGLL()

        ! create equdistributed nodes on equilateral triangle
        sk=1
        do i=1,NGLL
            do j=1,NGLL+1-i
                L1(sk) = (i-1.0) / p
                L3(sk) = (j-1.0) / p
                sk=sk+1
            end do
        end do
        L2(:) = 1.-L1(:)-L3(:)
        x(:) = -L2 + L3
        y(:) = (-L2(:)-L3(:)+2.*L1(:)) / sqrt(3.)

        ! compute blending
        blend1(:) = 4.0*L2(:)*L3(:)
        blend2(:) = 4.0*L1(:)*L3(:)
        blend3(:) = 4.0*L1(:)*L2(:)

        !get warpfactors
        warpfactor1 = warpfactor(xgll, L3(:)-L2(:))
        warpfactor2 = warpfactor(xgll, L1(:)-L3(:))
        warpfactor3 = warpfactor(xgll, L2(:)-L1(:))
    
        ! combine blend and warp
        warp1(:) = blend1(:)*warpfactor1(:)*(1.0+(alpha*L1(:))**2)
        warp2(:) = blend2(:)*warpfactor2(:)*(1.0+(alpha*L2(:))**2)
        warp3(:) = blend3(:)*warpfactor3(:)*(1.0+(alpha*L3(:))**2)
    
        ! accumulate deformations associated with each edge
        x(:) = x(:) + 1.0*warp1(:) + cos(2.0*PI/3.0)*warp2(:) + cos(4.0*PI/3.0)*warp3(:)
        y(:) = y(:) + 0.0*warp1(:) + sin(2.0*PI/3.0)*warp2(:) + sin(4.0*PI/3.0)*warp3(:)
    end subroutine nodes
end module nodesMod
