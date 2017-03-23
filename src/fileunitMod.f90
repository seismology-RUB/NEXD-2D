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
module fileunitMod
    use constantsMod

    implicit none

    contains

    function getFileUnit(fu_max)
        ! get_file_unit returns a unit number that is not in use
        integer :: getfileunit
        integer fu, fu_max, m, iostat
        logical opened

        m = fu_max ; if (m < 1) m = 97
        do fu = m,1,-1
            inquire (unit=fu, opened=opened, iostat=iostat)
            if (iostat.ne.0) cycle
            if (.not.opened) exit
        end do
        getFileUnit = fu
    end function getFileUnit
end module fileunitMod
