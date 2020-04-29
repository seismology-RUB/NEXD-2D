!-----------------------------------------------------------------------
!   Copyright 2020 Wolfgang Friederich (Ruhr-Universität Bochum, GER)
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
!--------------------------------------------------
!  Module with mathematical constants
!--------------------------------------------------
module mathConstants
    complex,parameter :: mc_ci = (0., 1.)                             ! complex i
    real,parameter :: mc_pi = 3.141592653589793                       ! pi single
    real,parameter :: mc_deg2rad = 3.141592653589793/180.             ! degree to radian
    !
    !  doubles
    !
    double precision,parameter :: mc_pid = 3.141592653589793          ! pi double precision
    double precision,parameter :: mc_ed  = 2.718281828459045          ! e double precision
!
end module mathConstants 
