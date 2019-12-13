!-----------------------------------------------------------------------
!   Copyright 2014-2019 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2019 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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
module timestampMod
    use parameterMod
    use calendar
    use convert_time

    implicit none

    type :: timestampVar
        real(kind=SIZE_DOUBLE) :: start
        real(kind=SIZE_DOUBLE) :: current
        integer :: year
        integer :: month
        integer :: day
        integer :: hours
        integer :: minutes
        integer :: timestamp
    end type

    contains

    subroutine initialTimestamp(timestamp)
        !in/out
        type(timestampVar) :: timestamp
        !local
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer, dimension(8) :: time_values

        call date_and_time(date,time,zone,time_values)
        timestamp%year    = time_values(1)
        timestamp%month   = time_values(2)
        timestamp%day     = time_values(3)
        timestamp%hours   = time_values(5)
        timestamp%minutes = time_values(6)
        call convtime(timestamp%timestamp,timestamp%year,timestamp%month,timestamp%day,timestamp%hours,timestamp%minutes)

        !convert the current time to seconds
        timestamp%start = timestamp%timestamp*dble(60) + time_values(7) + time_values(8)/dble(1000)
    end subroutine

    subroutine outputStamp(par, localtime, tstep, timestamp, normu, normv)
        !input
        type(parameterVar) :: par
        type(timestampVar) :: timestamp
        real(kind=custom_real), intent(in) :: normu
        real(kind=custom_real), intent(in) :: normv
        real(kind=custom_real), intent(in) :: localtime
        integer :: tstep
        !Local
        real(kind=SIZE_DOUBLE) :: cpu
        real(kind=SIZE_DOUBLE) :: remaining
        real(kind=SIZE_DOUBLE) :: total
        integer :: hours
        integer :: minutes
        integer :: seconds
        integer :: remaining_hours
        integer :: remaining_minutes
        integer :: remaining_seconds
        integer :: total_hours
        integer :: total_minutes
        integer :: total_seconds
        integer :: int_remaining
        integer :: int_total
        integer :: int_cpu
        !Variables to determine date and time at which the run will finish
        character(len=8)  :: date
        character(len=10) :: time
        character(len=5)  :: zone
        integer, dimension(8) :: time_values
        character(len=3), dimension(12) :: month_name
        character(len=3), dimension(0:6) :: weekday_name
        data month_name /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
        data weekday_name /'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/
        integer :: julian_day_number
        integer :: day_of_week

        !Create timestamp
        call date_and_time(date,time,zone,time_values)
        timestamp%year    = time_values(1)
        timestamp%month   = time_values(2)
        timestamp%day     = time_values(3)
        timestamp%hours   = time_values(5)
        timestamp%minutes = time_values(6)
        call convtime(timestamp%timestamp,timestamp%year,timestamp%month,timestamp%day,timestamp%hours,timestamp%minutes)

        !convert the current time to seconds
        timestamp%current = timestamp%timestamp*dble(60) + time_values(7) + time_values(8)/dble(1000)

        !Calculate elapsed time
        cpu     = timestamp%current - timestamp%start
        int_cpu = int(cpu)
        hours   = int_cpu/3600
        minutes = (int_cpu - 3600*hours)/60
        seconds = int_cpu - 3600*hours - 60*minutes

        !Calculate estimated remaining simulation time
        remaining         = (par%nt - tstep)*(cpu/dble(tstep))
        int_remaining     = int(remaining)
        remaining_hours   = int_remaining/3600
        remaining_minutes = (int_remaining - 3600*remaining_hours)/60
        remaining_seconds = int_remaining - 3600*remaining_hours-60*remaining_minutes

        !Calculate estimated total simulation time
        total         = remaining + cpu
        int_total     = int(total)
        total_hours   = int_total/3600
        total_minutes = (int_total - 3600*total_hours)/60
        total_seconds = int_total - 3600*total_hours - 60*total_minutes

        if (par%log) then
            write(*,"('|           Time step number ',i7,' (t = ',es8.1,' s) out of ',i7, '           |')") tstep,localtime,par%nt
            write(*,"('|           Elapsed time:                ',i4,' h ',i2.2,' m ',i2.2,' s', '                      |')") hours,minutes,seconds
            write(*,"('|           Estimated remaining time:    ',i4,' h ',i2.2,' m ',i2.2,' s', '                      |')")&
                remaining_hours, remaining_minutes, remaining_seconds
            write(*,"('|           Estimated total time:        ',i4,' h ',i2.2,' m ',i2.2,' s', '                      |')") &
                 total_hours,total_minutes,total_seconds
            write(*,"('|           Mean elapsed time per timestep:   ',es9.2,' s', '                      |')") cpu/dble(tstep)
            write (*,"(a43, es14.7, a23)") "|           Current maximum norm of U    : ", normu, "                      |"
            write (*,"(a43, es14.7, a23)") "|           Current maximum norm of V    : ", normv, "                      |"
        endif

        if(tstep <= par%nt) then
            ! compute date and time at which the run should finish (useful for long runs)
            ! add remaining minutes and get date and time of that future timestamp in minutes
            timestamp%timestamp = int((timestamp%current + remaining) / dble(60))
            call invtime(timestamp%timestamp,timestamp%year,timestamp%month,timestamp%day,timestamp%hours,timestamp%minutes)

            ! convert to Julian day to get day of the week
            call calndr(timestamp%day,timestamp%month,timestamp%year,julian_day_number)
            day_of_week = idaywk(julian_day_number)
            if (par%log) then
                write(*,"('|           This run will finish on: ',a3,' ',a3,' ',i2.2,', ',i4.4,' at ',i2.2,':',i2.2, '                 |')") &
                weekday_name(day_of_week),month_name(timestamp%month),&
                timestamp%day,timestamp%year,timestamp%hours,timestamp%minutes
            endif
        endif
        if (par%log) then
            write(*,'(a24, i3, a53)') "|                       ", 100*tstep/par%nt , "% of the simulation completed.                     |"
            write(*,"(a80)") "|------------------------------------------------------------------------------|"
        endif
    end subroutine



end module timestampMod
