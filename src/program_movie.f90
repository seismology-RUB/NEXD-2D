!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2019 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
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
program movie_pro
    use parameterMod
    use collectMovieMod
    use errorMessage
    use slipInterfaceMod

    type (parameterVar) :: par
    type (lsi_parameter) :: lsipar
    type (movie_parameter) :: movie
    type (error_message) :: errmsg
    integer :: myrank = 0
    character(len=13) :: myname = "program_movie"

    !Create new error message for this program
    call new(errmsg, myname)

    ! read Parfile
    call readParfile(par, lsipar, movie, myrank, errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif

    if (par%log) then
       write (*, "(a80)") "|             Generate single movie files for mpi version of dg2d              |"
       write (*, "(a80)") "|------------------------------------------------------------------------------|"
    end if

    call collectMovie(par, movie, errmsg)

end program movie_pro
