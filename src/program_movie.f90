!--------------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2014-2017 Thomas Möller (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
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
program movie_pro
    use parameterMod
    use collectMovieMod
    use errorMessage

    type (parameterVar) :: par
    type (movie_parameter) :: movie
    type (error_message) :: errmsg
    logical :: log
    integer :: myrank = 0
    character(len=13) :: myname = "program_movie"

    !Create new error message for this program
    call new(errmsg, myname)

    ! read Parfile
    call readParfile(par, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    ! log?

    if (par%log) then
       write (*, "(a80)") "|             Generate single movie files for mpi version of dg2d              |"
       write (*, "(a80)") "|------------------------------------------------------------------------------|"
    end if

    call collectMovie(par, movie, errmsg)

end program movie_pro
