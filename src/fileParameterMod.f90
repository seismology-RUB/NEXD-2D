!-----------------------------------------------------------------------
!   Copyright 2011-2016 Lasse Lambrecht (Ruhr-Universität Bochum, GER)
!   Copyright 2015-2018 Andre Lamert (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Thomas Möller (Ruhr-Universität Bochum, GER)
!   Copyright 2014-2018 Marc S. Boxberg (Ruhr-Universität Bochum, GER)
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
module fileParameterMod
    use constantsMod
    use errorMessage

    implicit none

    contains

    function getParameterName(line) result(par)
        !Function to call and compare the parameter name to the entry in the file -> find the loocation of a specific parameter!
        character (len=*), intent(in) :: line
        character (len=len(line)) :: tmp
        character (len=50) :: par

        integer, parameter :: tab = 9, space = 32, eq = 61
        integer :: i, ASCII_char, n

        n = 0
        par = " "
        tmp = line(1:index(line, "="))
        do i=1, len(tmp)
            ASCII_char = iachar(tmp(i:i))
            if (ASCII_char /= tab .and. ASCII_char /= space .and. ASCII_char /= eq) then
                n = n + 1
                par(n:n) = tmp(i:i)
            end if
        end do
    end function getParameterName

    function getParameterValue(name, filename, pos, errmsg) result(value)
        !Read the parameters value from the input string
        type(error_message) :: errmsg
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: pos

        character (len=80) :: tmp, par, value
        character (len=255) :: line
        character (len=17) :: myname = 'getParameterValue'


        integer :: i, j, k, ier=0
        integer, parameter :: seek_set = 0

        integer :: set_pos
        integer :: ftell

        do while (.not. is_iostat_end(ier))
            read(19, "(a255)", iostat=ier) line
            if (is_iostat_end(ier)) then
                call add(errmsg, 2, "Parameter, "//trim(name)//" , not found in parameter file...", myname, filename )
            endif
            line = trim(adjustl(line))
            if (line(1:1) == "#") cycle
            if (line(1:1) == "\n") cycle
            par = getParameterName(line)
            par = trim(par)
            j = scan(line, "=")
            if (par /= name) cycle
            if (par == name) then
                do k = j, len_trim(line)
                    if (line(k:k) /= "=" .and. line(k:k) /= " ") then
                        i=scan(line,'#')
                        if (i==0) i=len_trim(line)+1
                        tmp = trim(line(k:i-1))
                        read(tmp, "(a80)") value
                        value = trim(value)
                        exit
                    endif
                enddo
                exit
            endif
        enddo
        call fseek(19, pos, seek_set, ier)
        set_pos = int(ftell(19))
    end function getParameterValue

    subroutine readStringPar(par_to_read, name, filename, pos, errmsg)
        !Function to read Parametrs as text.
        type(error_message) :: errmsg
        character (len=*) :: name, par_to_read, filename
        character (len=80) :: value!, tmp
        integer :: pos, l

        value = getParameterValue(name, filename, pos, errmsg)
        read(value, "(a80)") par_to_read
        l = scan(par_to_read, "#")
        if (l /= 0) par_to_read = trim(par_to_read(:l-1))
    end subroutine readStringPar

    subroutine readIntPar(par_to_read, name, filename, pos, errmsg)
        !function to read interger values from the parameter-file
        type(error_message) :: errmsg
        character (len=*), intent(in) :: name, filename
        integer :: par_to_read, pos, ios
        character (len=80) :: value

        value = getParameterValue(name, filename, pos, errmsg)
        read(value, *, iostat=ios) par_to_read
        if (ios > 0) call add(errmsg, 2, "Type Error for parameter "//name//". Parameter is supposed to be an integer.", "readIntPar", filename )
    end subroutine readIntPar

    subroutine readFloatPar(par_to_read, name, filename, pos, errmsg)
        !Function to read floating point variables from the parameter-file
        type(error_message) :: errmsg
        character (len=*), intent(in) :: name, filename
        real(kind=CUSTOM_REAL) :: par_to_read
        character (len=80) :: value
        integer :: pos, ios

        value = getParameterValue(name, filename, pos, errmsg)
        read(value, *, iostat=ios) par_to_read
        if (ios > 0) call add(errmsg, 2, "Type Error for parameter "//name//". Parameter is supposed to be a float.", "readFloatPar", filename )
    end subroutine readFloatPar

    subroutine readLogicalPar(par_to_read, name, filename, pos, errmsg)
        !Function to read logical variables from the parameter-file
        type(error_message) :: errmsg
        character (len=*), intent(in) :: name, filename
        logical :: par_to_read
        character (len=7) :: value
        integer :: pos, ios

        write(value,'(a7)') getParameterValue(name, filename, pos, errmsg)
        read(value, *, iostat=ios) par_to_read
        if (ios > 0) call add(errmsg, 2, "Type Error for parameter "//name//". Parameter is supposed to be an bool.", "readLogicalPar", filename )
    end subroutine

    function setFilePosition(name, filename, nr, errmsg) result(seek_pos)
        !Function to set the starting position from which the file is read to a specific postion
        type(error_message) :: errmsg
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: nr

        character (len=80) :: tmp, par, value
        character (len=255) :: line, errstr
        character (len=15) :: myname = 'setFilePosition'
        integer :: seek_pos, ier, tmp_int, j, k
        integer :: ftell

        do
            seek_pos = -1
            read(19, "(a255)", iostat=ier) line
            if (is_iostat_end(ier)) then
                write(errstr, * ) "End of File reached. Parameter "//trim(name)//" with value", nr, "is not found."
                call add(errmsg, 2, errstr, myname, filename )
                call print(errmsg)
                stop
            endif
            line = trim(adjustl(line))
            if (line(1:1) == "#") cycle
            if (line(1:1) == "\n") cycle
            par = getParameterName(line)
            par = trim(par)
            j = scan(line, "=")
            if (par /= name) cycle
            if (par == name) then
                do k = j, len_trim(line)
                    if (line(k:k) /= "=" .and. line(k:k) /= " ") then
                        tmp = trim(line(k:len_trim(line)))
                        read(tmp, *) value
                        value = trim(value)
                        read(value,*) tmp_int
                        exit
                    endif
                enddo
            else
                call add(errmsg, 2, "Parameter, "//trim(name)//" , not found in parameter file...", myname, filename )
            endif
            if (tmp_int == nr) then
                seek_pos = int(ftell(19))
                !Exit the loop if a matching source has been found
                exit
            endif
        enddo
    end function setFilePosition
end module fileParameterMod
