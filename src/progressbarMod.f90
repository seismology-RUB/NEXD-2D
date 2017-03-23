!--------------------------------------------------------------------------
!   Copyright 2014-2017 Marc S. Boxberg (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
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
module progressbarMod
    
    implicit none
    
    contains
    
    subroutine progress(j,total)
        !input
        integer :: j
        integer :: total
        !local
        integer :: k
        character(len=70) :: bar

        bar="Progress: [                                                  ] ???%   "

        write(unit=bar(64:66),fmt="(i3)") 100*j/total

        do k = 1, j*50/total
            bar(11+k:11+k)="#"
        enddo
        ! print the progress bar.
        write(unit=6,fmt="(a1,a70)",advance="no") char(13), bar
        if (j /= total) then
            flush(unit=6)
        else
            write(unit=6, fmt=*)
        endif
        return
    end subroutine progress

end module progressbarMod
